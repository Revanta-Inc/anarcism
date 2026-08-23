//! Convert a versioned ANARCI germline export into a compact deterministic
//! binary representation. Aligned residues use a five-bit alphabet.

use std::collections::BTreeMap;
use std::env;
use std::fmt::{self, Display};
use std::fs;
use std::path::PathBuf;

const MAGIC: &[u8; 8] = b"ANRCGER1";
const FORMAT_VERSION: u16 = 1;
const ALIGNMENT_LENGTH: usize = 128;
const PACKED_LENGTH: usize = ALIGNMENT_LENGTH * 5 / 8;
const ALPHABET: &[u8; 21] = b"-ACDEFGHIKLMNPQRSTVWY";

#[derive(Debug)]
struct Error(String);

impl Display for Error {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl std::error::Error for Error {}

type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
struct Germline<'a> {
    segment: u8,
    chain: u8,
    species: &'a str,
    gene: &'a str,
    sequence: &'a [u8],
}

fn main() {
    if let Err(error) = run() {
        eprintln!("anarcism-germlinegen: {error}");
        std::process::exit(2);
    }
}

fn run() -> Result<()> {
    let mut args = env::args_os().skip(1);
    let input = args.next().map(PathBuf::from).ok_or_else(usage)?;
    let output = args.next().map(PathBuf::from).ok_or_else(usage)?;
    if args.next().is_some() {
        return Err(usage());
    }

    let source = fs::read_to_string(&input)
        .map_err(|error| Error(format!("could not read {}: {error}", input.display())))?;
    let germlines = parse(&source)?;
    let encoded = encode(&germlines)?;
    fs::write(&output, &encoded)
        .map_err(|error| Error(format!("could not write {}: {error}", output.display())))?;
    println!(
        "encoded {} germlines ({} bytes) from {}",
        germlines.len(),
        encoded.len(),
        input.display()
    );
    Ok(())
}

fn usage() -> Error {
    Error("usage: germlinegen <germlines.tsv> <germlines.bin>".into())
}

fn parse(source: &str) -> Result<Vec<Germline<'_>>> {
    let mut result = Vec::new();
    for (line_index, line) in source.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 5 {
            return Err(Error(format!(
                "line {} has {} fields; expected 5",
                line_index + 1,
                fields.len()
            )));
        }
        let segment = one_of(fields[0], b"VJ", "segment", line_index)?;
        let chain = one_of(fields[1], b"HKLABGD", "chain", line_index)?;
        let sequence = fields[4].as_bytes();
        if sequence.len() != ALIGNMENT_LENGTH {
            return Err(Error(format!(
                "line {} alignment has length {}; expected {ALIGNMENT_LENGTH}",
                line_index + 1,
                sequence.len()
            )));
        }
        if sequence.iter().any(|residue| !ALPHABET.contains(residue)) {
            return Err(Error(format!(
                "line {} alignment contains a noncanonical residue",
                line_index + 1
            )));
        }
        if fields[2].is_empty() || fields[3].is_empty() {
            return Err(Error(format!(
                "line {} has empty species or gene metadata",
                line_index + 1
            )));
        }
        result.push(Germline {
            segment,
            chain,
            species: fields[2],
            gene: fields[3],
            sequence,
        });
    }
    if result.is_empty() {
        return Err(Error("germline export is empty".into()));
    }
    Ok(result)
}

fn one_of(value: &str, allowed: &[u8], field: &str, line_index: usize) -> Result<u8> {
    let bytes = value.as_bytes();
    if bytes.len() == 1 && allowed.contains(&bytes[0]) {
        Ok(bytes[0])
    } else {
        Err(Error(format!(
            "line {} has invalid {field} {value:?}",
            line_index + 1
        )))
    }
}

fn encode(germlines: &[Germline<'_>]) -> Result<Vec<u8>> {
    let mut species_indices = BTreeMap::new();
    let mut species_names = Vec::new();
    for germline in germlines {
        if !species_indices.contains_key(germline.species) {
            let index: u8 = species_indices
                .len()
                .try_into()
                .map_err(|_| Error("more than 255 germline species".into()))?;
            species_indices.insert(germline.species, index);
            species_names.push(germline.species);
        }
    }

    let mut output = Vec::with_capacity(16 + germlines.len() * (PACKED_LENGTH + 20));
    output.extend_from_slice(MAGIC);
    push_u16(&mut output, FORMAT_VERSION);
    output.push(
        species_indices
            .len()
            .try_into()
            .map_err(|_| Error("more than 255 germline species".into()))?,
    );
    push_u16(
        &mut output,
        germlines
            .len()
            .try_into()
            .map_err(|_| Error("more than 65535 germlines".into()))?,
    );
    for species in species_names {
        push_short_string(&mut output, species)?;
    }

    for germline in germlines {
        output.push(germline.segment);
        output.push(germline.chain);
        output.push(species_indices[germline.species]);
        push_short_string(&mut output, germline.gene)?;
        output.extend_from_slice(&pack_sequence(germline.sequence)?);
    }
    Ok(output)
}

fn pack_sequence(sequence: &[u8]) -> Result<[u8; PACKED_LENGTH]> {
    let mut packed = [0_u8; PACKED_LENGTH];
    let mut bit_offset = 0;
    for residue in sequence {
        let value = ALPHABET
            .iter()
            .position(|candidate| candidate == residue)
            .ok_or_else(|| Error("alignment contains a noncanonical residue".into()))?
            as u16;
        let byte_offset = bit_offset / 8;
        let shift = bit_offset % 8;
        packed[byte_offset] |= (value << shift) as u8;
        if shift > 3 {
            packed[byte_offset + 1] |= (value >> (8 - shift)) as u8;
        }
        bit_offset += 5;
    }
    Ok(packed)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn five_bit_sequence_round_trip() {
        let sequence: Vec<u8> = ALPHABET.iter().copied().cycle().take(128).collect();
        let packed = pack_sequence(&sequence).unwrap();
        for (position, expected) in sequence.into_iter().enumerate() {
            let bit_offset = position * 5;
            let byte_offset = bit_offset / 8;
            let shift = bit_offset % 8;
            let mut value = u16::from(packed[byte_offset]) >> shift;
            if shift > 3 {
                value |= u16::from(packed[byte_offset + 1]) << (8 - shift);
            }
            assert_eq!(ALPHABET[usize::from(value & 0x1f)], expected);
        }
    }
}
