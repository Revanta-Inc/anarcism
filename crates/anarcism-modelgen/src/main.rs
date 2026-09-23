//! Convert a versioned ANARCI-generated HMMER3/f database into the compact,
//! deterministic representation embedded by the engine.

use std::env;
use std::fmt::{self, Display};
use std::fs;
use std::path::PathBuf;

const MAGIC: &[u8; 8] = b"ANRCPRF3";
const FORMAT_VERSION: u16 = 3;
const SCORE_SCALE: f32 = 32_768.0;
const MSV_SCORE_SCALE: f32 = 3.0 / std::f32::consts::LN_2;
const MODEL_LENGTH: usize = 128;
const ALPHABET: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
const BACKGROUND: [f32; 20] = [
    0.078_794_5,
    0.015_160_0,
    0.053_522_2,
    0.066_829_8,
    0.039_706_2,
    0.069_507_1,
    0.022_919_8,
    0.059_009_2,
    0.059_442_2,
    0.096_372_8,
    0.023_771_8,
    0.041_438_6,
    0.048_290_4,
    0.039_563_9,
    0.054_097_8,
    0.068_336_4,
    0.054_068_7,
    0.067_341_7,
    0.011_413_5,
    0.030_413_3,
];

#[derive(Debug)]
struct Error(String);

impl Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Error {}

type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
struct RawProfile {
    name: String,
    length: usize,
    checksum: u32,
    msv_mu: f32,
    msv_lambda: f32,
    forward_tau: f32,
    forward_lambda: f32,
    consensus: Vec<u8>,
    match_nlog: Vec<[f32; 20]>,
    transitions_nlog: Vec<[f32; 7]>,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("anarcism-modelgen: {error}");
        std::process::exit(2);
    }
}

fn run() -> Result<()> {
    let mut args = env::args_os().skip(1);
    let input = args
        .next()
        .map(PathBuf::from)
        .ok_or_else(|| Error("usage: anarcism-modelgen <ANARCI ALL.hmm> <profiles.bin>".into()))?;
    let output = args
        .next()
        .map(PathBuf::from)
        .ok_or_else(|| Error("usage: anarcism-modelgen <ANARCI ALL.hmm> <profiles.bin>".into()))?;
    if args.next().is_some() {
        return Err(Error(
            "usage: anarcism-modelgen <ANARCI ALL.hmm> <profiles.bin>".into(),
        ));
    }

    let source = fs::read_to_string(&input)
        .map_err(|error| Error(format!("could not read {}: {error}", input.display())))?;
    let profiles = parse_hmmer3(&source)?;
    let encoded = encode(&profiles)?;
    fs::write(&output, &encoded)
        .map_err(|error| Error(format!("could not write {}: {error}", output.display())))?;

    println!(
        "encoded {} profiles ({} bytes) from {}",
        profiles.len(),
        encoded.len(),
        input.display()
    );
    Ok(())
}

fn parse_hmmer3(source: &str) -> Result<Vec<RawProfile>> {
    let lines: Vec<&str> = source.lines().collect();
    let mut cursor = 0;
    let mut profiles = Vec::new();

    while cursor < lines.len() {
        if !lines[cursor].starts_with("HMMER3/f") {
            cursor += 1;
            continue;
        }
        cursor += 1;

        let mut name = None;
        let mut length = None;
        let mut checksum = 0;
        let mut msv = None;
        let mut forward = None;

        while cursor < lines.len() && !lines[cursor].starts_with("HMM ") {
            let line = lines[cursor];
            cursor += 1;
            if let Some(value) = line.strip_prefix("NAME  ") {
                name = Some(value.trim().to_owned());
            } else if let Some(value) = line.strip_prefix("LENG  ") {
                length = Some(parse_usize(value, "LENG")?);
            } else if let Some(value) = line.strip_prefix("CKSUM ") {
                checksum = value
                    .trim()
                    .parse()
                    .map_err(|error| Error(format!("invalid CKSUM {value:?}: {error}")))?;
            } else if let Some(value) = line.strip_prefix("STATS LOCAL MSV") {
                let fields: Vec<&str> = value.split_whitespace().collect();
                if fields.len() != 2 {
                    return Err(Error(format!("invalid MSV stats line: {line}")));
                }
                msv = Some((
                    parse_f32(fields[0], "MSV mu")?,
                    parse_f32(fields[1], "MSV lambda")?,
                ));
            } else if let Some(value) = line.strip_prefix("STATS LOCAL FORWARD") {
                let fields: Vec<&str> = value.split_whitespace().collect();
                if fields.len() != 2 {
                    return Err(Error(format!("invalid FORWARD stats line: {line}")));
                }
                forward = Some((
                    parse_f32(fields[0], "FORWARD tau")?,
                    parse_f32(fields[1], "FORWARD lambda")?,
                ));
            }
        }

        let name = name.ok_or_else(|| Error("profile is missing NAME".into()))?;
        let length = length.ok_or_else(|| Error(format!("profile {name} is missing LENG")))?;
        if length != MODEL_LENGTH {
            return Err(Error(format!(
                "profile {name} has length {length}; expected {MODEL_LENGTH}"
            )));
        }
        let (msv_mu, msv_lambda) =
            msv.ok_or_else(|| Error(format!("profile {name} is missing MSV stats")))?;
        let (forward_tau, forward_lambda) =
            forward.ok_or_else(|| Error(format!("profile {name} is missing FORWARD stats")))?;
        if !msv_mu.is_finite() || !msv_lambda.is_finite() || msv_lambda <= 0.0 {
            return Err(Error(format!("profile {name} has invalid MSV stats")));
        }
        if !forward_tau.is_finite() || !forward_lambda.is_finite() || forward_lambda <= 0.0 {
            return Err(Error(format!("profile {name} has invalid FORWARD stats")));
        }

        let alphabet_line = lines
            .get(cursor)
            .ok_or_else(|| Error(format!("profile {name} is missing HMM header")))?;
        let observed_alphabet: Vec<u8> = alphabet_line
            .split_whitespace()
            .skip(1)
            .filter(|field| field.len() == 1)
            .map(|field| field.as_bytes()[0])
            .collect();
        if observed_alphabet != ALPHABET {
            return Err(Error(format!("profile {name} uses an unexpected alphabet")));
        }
        cursor += 1;

        cursor += 1;
        require_line(&lines, cursor, &name, "COMPO")?;
        cursor += 1;
        require_line(&lines, cursor, &name, "node-0 insert emissions")?;
        cursor += 1;
        let mut transitions_nlog = Vec::with_capacity(length + 1);
        transitions_nlog.push(parse_scores::<7>(
            require_line(&lines, cursor, &name, "node-0 transitions")?,
            "node-0 transitions",
        )?);
        cursor += 1;

        let mut consensus = Vec::with_capacity(length);
        let mut match_nlog = Vec::with_capacity(length);
        for expected_node in 1..=length {
            let node_line = require_line(&lines, cursor, &name, "match emissions")?;
            cursor += 1;
            let fields: Vec<&str> = node_line.split_whitespace().collect();
            if fields.len() < 22 {
                return Err(Error(format!("truncated node line in {name}: {node_line}")));
            }
            let observed_node = parse_usize(fields[0], "node index")?;
            if observed_node != expected_node {
                return Err(Error(format!(
                    "profile {name}: expected node {expected_node}, found {observed_node}"
                )));
            }
            match_nlog.push(parse_scores::<20>(
                &fields[1..21].join(" "),
                "match emissions",
            )?);
            consensus.push(fields[21].as_bytes()[0].to_ascii_uppercase());

            // Protein insert emissions have zero configured log odds.
            require_line(&lines, cursor, &name, "insert emissions")?;
            cursor += 1;
            transitions_nlog.push(parse_scores::<7>(
                require_line(&lines, cursor, &name, "transitions")?,
                "transitions",
            )?);
            cursor += 1;
        }

        while cursor < lines.len() && lines[cursor].trim().is_empty() {
            cursor += 1;
        }
        if lines.get(cursor).map(|line| line.trim()) != Some("//") {
            return Err(Error(format!("profile {name} is missing // terminator")));
        }
        cursor += 1;

        profiles.push(RawProfile {
            name,
            length,
            checksum,
            msv_mu,
            msv_lambda,
            forward_tau,
            forward_lambda,
            consensus,
            match_nlog,
            transitions_nlog,
        });
    }

    if profiles.is_empty() {
        return Err(Error("no HMMER3/f profiles found".into()));
    }
    Ok(profiles)
}

fn require_line<'a>(lines: &'a [&str], cursor: usize, name: &str, field: &str) -> Result<&'a str> {
    lines
        .get(cursor)
        .copied()
        .ok_or_else(|| Error(format!("profile {name} is missing {field}")))
}

fn parse_scores<const N: usize>(line: &str, context: &str) -> Result<[f32; N]> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < N {
        return Err(Error(format!(
            "{context}: expected {N} scores, found {}",
            fields.len()
        )));
    }
    let mut values = [0.0; N];
    for (index, field) in fields.iter().take(N).enumerate() {
        values[index] = if *field == "*" {
            f32::INFINITY
        } else {
            parse_f32(field, context)?
        };
    }
    Ok(values)
}

fn parse_f32(value: &str, context: &str) -> Result<f32> {
    value
        .parse()
        .map_err(|error| Error(format!("invalid {context} value {value:?}: {error}")))
}

fn parse_usize(value: &str, context: &str) -> Result<usize> {
    value
        .trim()
        .parse()
        .map_err(|error| Error(format!("invalid {context} value {value:?}: {error}")))
}

fn encode(profiles: &[RawProfile]) -> Result<Vec<u8>> {
    let mut output = Vec::with_capacity(16 + profiles.len() * 7_400);
    output.extend_from_slice(MAGIC);
    push_u16(&mut output, FORMAT_VERSION);
    push_u16(&mut output, SCORE_SCALE as u16);
    push_u16(
        &mut output,
        profiles
            .len()
            .try_into()
            .map_err(|_| Error("too many profiles".into()))?,
    );
    push_u16(&mut output, MODEL_LENGTH as u16);

    for profile in profiles {
        let (species, chain) = profile.name.rsplit_once('_').ok_or_else(|| {
            Error(format!(
                "profile name {:?} is not species_chain",
                profile.name
            ))
        })?;
        if chain.len() != 1
            || !matches!(
                chain.as_bytes()[0],
                b'H' | b'K' | b'L' | b'A' | b'B' | b'G' | b'D'
            )
        {
            return Err(Error(format!(
                "unsupported chain in profile {}",
                profile.name
            )));
        }
        push_short_string(&mut output, &profile.name)?;
        push_short_string(&mut output, species)?;
        output.push(chain.as_bytes()[0]);
        output.push(u8::from(matches!(
            chain.as_bytes()[0],
            b'A' | b'B' | b'G' | b'D'
        )));
        push_u32(&mut output, profile.checksum);
        output.extend_from_slice(&profile.msv_mu.to_le_bytes());
        output.extend_from_slice(&profile.msv_lambda.to_le_bytes());
        output.extend_from_slice(&profile.forward_tau.to_le_bytes());
        output.extend_from_slice(&profile.forward_lambda.to_le_bytes());
        output.extend_from_slice(&profile.consensus);

        let match_scores: Vec<[f32; 20]> = profile
            .match_nlog
            .iter()
            .map(|node| std::array::from_fn(|residue| -node[residue] - BACKGROUND[residue].ln()))
            .collect();
        let (msv_bias, msv_match_costs) = msv_match_costs(&match_scores);
        output.push(msv_bias);
        output.extend_from_slice(&msv_match_costs);

        let entries = local_entry_scores(&profile.transitions_nlog, profile.length);
        for score in entries {
            push_i24(&mut output, quantize(score));
        }

        for node in &match_scores {
            for score in node {
                push_i24(&mut output, quantize(*score));
            }
        }

        // DP transitions exist only for nodes 1..M-1.
        for node in profile.transitions_nlog.iter().take(profile.length).skip(1) {
            for negative_log_probability in node {
                push_i24(&mut output, quantize(-*negative_log_probability));
            }
        }
    }
    Ok(output)
}

/// Reproduce HMMER 3.4's one-third-bit `mf_conversion()`.
fn msv_match_costs(match_scores: &[[f32; 20]]) -> (u8, Vec<u8>) {
    let maximum = match_scores
        .iter()
        .flatten()
        .copied()
        .fold(0.0_f32, f32::max);
    let bias = msv_unbiased_byte_cost(-maximum);
    // Residue-major storage avoids gathers in the MSV recurrence.
    let costs = (0..ALPHABET.len())
        .flat_map(|residue| {
            match_scores
                .iter()
                .map(move |node| msv_biased_byte_cost(node[residue], bias))
        })
        .collect();
    (bias, costs)
}

fn msv_unbiased_byte_cost(score: f32) -> u8 {
    let cost = -(MSV_SCORE_SCALE * score).round();
    cost.clamp(0.0, 255.0) as u8
}

fn msv_biased_byte_cost(score: f32, bias: u8) -> u8 {
    if !score.is_finite() {
        return u8::MAX;
    }
    let cost = -(MSV_SCORE_SCALE * score).round() as i32 + i32::from(bias);
    cost.clamp(0, i32::from(u8::MAX)) as u8
}

fn local_entry_scores(transitions: &[[f32; 7]], length: usize) -> Vec<f32> {
    const MM: usize = 0;
    const MI: usize = 1;
    const IM: usize = 3;
    const DM: usize = 5;

    let probability = |negative_log: f32| {
        if negative_log.is_infinite() {
            0.0
        } else {
            (-negative_log).exp()
        }
    };

    let mut occupancy = vec![0.0_f32; length + 1];
    occupancy[1] = probability(transitions[0][MI]) + probability(transitions[0][MM]);
    for node in 2..=length {
        let previous = occupancy[node - 1];
        occupancy[node] = previous
            * (probability(transitions[node - 1][MM]) + probability(transitions[node - 1][MI]))
            + (1.0 - previous) * probability(transitions[node - 1][DM]);
    }

    // Insert occupancy uses I->M; local match entry does not.
    let _ = IM;
    let normalizer: f32 = (1..=length)
        .map(|node| occupancy[node] * (length - node + 1) as f32)
        .sum();
    (1..=length)
        .map(|node| (occupancy[node] / normalizer).ln())
        .collect()
}

fn quantize(score: f32) -> i32 {
    if !score.is_finite() {
        -(1 << 23)
    } else {
        (score * SCORE_SCALE)
            .round()
            .clamp((-(1 << 23) + 1) as f32, ((1 << 23) - 1) as f32) as i32
    }
}

fn push_short_string(output: &mut Vec<u8>, value: &str) -> Result<()> {
    let length: u8 = value
        .len()
        .try_into()
        .map_err(|_| Error(format!("string is too long: {value:?}")))?;
    output.push(length);
    output.extend_from_slice(value.as_bytes());
    Ok(())
}

fn push_u16(output: &mut Vec<u8>, value: u16) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn push_i24(output: &mut Vec<u8>, value: i32) {
    output.extend_from_slice(&value.to_le_bytes()[..3]);
}

fn push_u32(output: &mut Vec<u8>, value: u32) {
    output.extend_from_slice(&value.to_le_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantization_reserves_minimum_for_impossible_scores() {
        assert_eq!(quantize(f32::NEG_INFINITY), -(1 << 23));
        assert_ne!(quantize(-32.0), -(1 << 23));
    }

    #[test]
    fn expected_hmmer_alphabet_is_stable() {
        assert_eq!(
            std::str::from_utf8(ALPHABET).unwrap(),
            "ACDEFGHIKLMNPQRSTVWY"
        );
    }

    #[test]
    fn msv_byte_costs_use_hmmer_third_bit_rounding_and_bias() {
        let scores = [[2.0_f32; 20], [-1.0_f32; 20]];
        let (bias, costs) = msv_match_costs(&scores);
        assert_eq!(bias, (MSV_SCORE_SCALE * 2.0).round() as u8);
        assert_eq!(costs[0], 0);
        assert_eq!(costs[1], bias + MSV_SCORE_SCALE.round() as u8);
        assert_eq!(costs[2], 0);
        assert_eq!(msv_biased_byte_cost(f32::NEG_INFINITY, bias), u8::MAX);
    }
}
