use anarcism_core::{
    DomainResult, GermlineAssignment, NumberingOptions, ProfileHit, number_sequence,
};
#[cfg(not(target_family = "wasm"))]
use anarcism_core::{SequenceInput, number_sequences, number_sequences_parallel};
use serde::Deserialize;
use std::collections::HashMap;
#[cfg(not(target_family = "wasm"))]
use std::num::NonZeroUsize;

#[derive(Deserialize)]
struct InputCase {
    id: String,
    seq: String,
    category: String,
}

#[derive(Deserialize)]
struct ReferenceCase {
    id: String,
    category: String,
    domains: Vec<ReferenceDomain>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct ReferenceDomain {
    profile: String,
    start: usize,
    end: usize,
    bit_score: f32,
    e_value: f64,
    bias: f32,
    alternative_hits: Vec<ReferenceHit>,
    germline: Option<ReferenceGermline>,
    numbering: String,
    alignment: String,
}

/// `[profile, bitScore, eValue, bias, queryStart, queryEnd]`
type HitRecord = (String, f32, f64, f32, usize, usize);

/// `[species, vGene, vIdentity, jGene, jIdentity]`
type GermlineRecord = (
    String,
    Option<String>,
    Option<f32>,
    Option<String>,
    Option<f32>,
);

#[derive(Deserialize)]
#[serde(from = "HitRecord")]
struct ReferenceHit {
    profile: String,
    bit_score: f32,
    e_value: f64,
    bias: f32,
    query_start: usize,
    query_end: usize,
}

impl From<HitRecord> for ReferenceHit {
    fn from((profile, bit_score, e_value, bias, query_start, query_end): HitRecord) -> Self {
        Self {
            profile,
            bit_score,
            e_value,
            bias,
            query_start,
            query_end,
        }
    }
}

#[derive(Deserialize)]
#[serde(from = "GermlineRecord")]
struct ReferenceGermline {
    species: String,
    v_gene: Option<String>,
    v_identity: Option<f32>,
    j_gene: Option<String>,
    j_identity: Option<f32>,
}

impl From<GermlineRecord> for ReferenceGermline {
    fn from((species, v_gene, v_identity, j_gene, j_identity): GermlineRecord) -> Self {
        Self {
            species,
            v_gene,
            v_identity,
            j_gene,
            j_identity,
        }
    }
}

/// Reads `corpus_reference.jsonl` by case id, skipping its header record.
fn reference_cases() -> HashMap<String, ReferenceCase> {
    include_str!("../../../tests/golden/corpus_reference.jsonl")
        .lines()
        .skip(1)
        .map(|line| {
            let case: ReferenceCase =
                serde_json::from_str(line).expect("reference record is valid JSON");
            (case.id.clone(), case)
        })
        .collect()
}

fn input_cases() -> Vec<InputCase> {
    serde_json::from_str(include_str!("../../../tests/golden/corpus.json"))
        .expect("checked-in corpus is valid JSON")
}

#[test]
fn corpus_matches_versioned_anarci_reference() {
    let inputs = input_cases();
    let reference = reference_cases();
    assert_eq!(inputs.len(), reference.len());

    let worker_count = std::thread::available_parallelism()
        .map_or(1, std::num::NonZero::get)
        .min(inputs.len().max(1));
    std::thread::scope(|scope| {
        let mut workers = Vec::with_capacity(worker_count);
        for shard_index in 0..worker_count {
            let inputs = &inputs;
            let reference_cases = &reference;
            workers.push(scope.spawn(move || {
                let mut failures = Vec::new();
                for case_index in (shard_index..inputs.len()).step_by(worker_count) {
                    if std::panic::catch_unwind(|| {
                        let input = &inputs[case_index];
                        let expected = reference_cases
                            .get(&input.id)
                            .unwrap_or_else(|| panic!("{} has no reference record", input.id));
                        assert_case(input, expected);
                    })
                    .is_err()
                    {
                        failures.push(inputs[case_index].id.clone());
                    }
                }
                failures
            }));
        }
        let mut failures = Vec::new();
        for worker in workers {
            failures.extend(worker.join().expect("expanded parity worker panicked"));
        }
        failures.sort_unstable();
        assert!(
            failures.is_empty(),
            "{} corpus cases failed parity: {}",
            failures.len(),
            failures.join(", ")
        );
    });
}

#[test]
#[cfg(not(target_family = "wasm"))]
fn full_native_batch_scheduler_matches_the_serial_batch() {
    let inputs: Vec<_> = input_cases()
        .into_iter()
        .map(|case| SequenceInput {
            id: case.id,
            sequence: case.seq,
        })
        .collect();
    let options = NumberingOptions::default();
    let serial = number_sequences(&inputs, &options).expect("serial corpus numbering succeeds");
    let parallel = number_sequences_parallel(
        &inputs,
        &options,
        NonZeroUsize::new(8).expect("worker count is nonzero"),
    )
    .expect("parallel corpus numbering succeeds");
    assert_eq!(parallel, serial);
}

fn assert_case(input: &InputCase, expected_case: &ReferenceCase) {
    assert_eq!(input.id, expected_case.id);
    assert_eq!(input.category, expected_case.category);
    let result = number_sequence(
        &input.seq,
        &NumberingOptions {
            assign_germline: true,
            alternative_hit_count: if strict_alternative_reference(&input.id) {
                28
            } else {
                NumberingOptions::default().alternative_hit_count
            },
            ..NumberingOptions::default()
        },
    )
    .unwrap_or_else(|error| panic!("{} failed: {error}", input.id));
    assert_eq!(
        result.domains.len(),
        expected_case.domains.len(),
        "{} domain count ({})",
        input.id,
        input.category
    );
    for (observed, expected) in result.domains.iter().zip(&expected_case.domains) {
        assert_eq!(
            observed_profile(observed),
            expected.profile,
            "{} profile",
            input.id
        );
        assert_eq!(observed.start, expected.start, "{} start", input.id);
        assert_eq!(observed.end, expected.end, "{} end", input.id);
        assert!(
            (observed.bit_score - expected.bit_score).abs() <= 0.2,
            "{} score: Rust {}, ANARCI {}",
            input.id,
            observed.bit_score,
            expected.bit_score
        );
        assert_e_value_parity(&input.id, observed.e_value, expected.e_value);
        assert!(
            (observed.bias - expected.bias).abs() <= 0.2,
            "{} bias: Rust {}, ANARCI {}",
            input.id,
            observed.bias,
            expected.bias
        );
        let matched_alternatives = expected
            .alternative_hits
            .iter()
            .filter(|expected_hit| {
                observed
                    .alternative_hits
                    .iter()
                    .any(|hit| hit.profile == expected_hit.profile)
            })
            .count();
        assert!(
            matched_alternatives >= expected.alternative_hits.len().min(1),
            "{} only matched {matched_alternatives}/{} pinned top alternative profiles",
            input.id,
            expected.alternative_hits.len()
        );
        if strict_alternative_reference(&input.id) {
            for (hit_index, (actual_hit, expected_hit)) in observed
                .alternative_hits
                .iter()
                .zip(&expected.alternative_hits)
                .enumerate()
            {
                assert_hit_parity(&input.id, hit_index, actual_hit, expected_hit);
            }
        }
        assert_germline_parity(
            &input.id,
            observed.germline.as_ref(),
            expected.germline.as_ref(),
        );
        assert_eq!(
            encode_numbering(input, observed),
            expected.numbering,
            "{} residue parity",
            input.id
        );
        assert_eq!(
            observed.padded_imgt_alignment, expected.alignment,
            "{} padded alignment",
            input.id
        );
    }
}

fn observed_profile(domain: &DomainResult) -> String {
    format!("{}_{:?}", domain.species, domain.chain_type)
}

/// Encodes numbering like the reference: IMGT labels in residue order, with
/// runs of consecutive plain positions written as `a-b`.
fn encode_numbering(input: &InputCase, domain: &DomainResult) -> String {
    let mut tokens: Vec<String> = Vec::new();
    let mut run: Option<(u16, u16)> = None;
    for (offset, residue) in domain.numbering.iter().enumerate() {
        assert_eq!(
            residue.sequence_index,
            domain.start + offset,
            "{} residue {offset} index",
            input.id
        );
        assert_eq!(
            input.seq[residue.sequence_index..].chars().next(),
            Some(residue.amino_acid),
            "{} residue {offset} amino acid",
            input.id
        );
        if residue.insertion_code.is_empty() {
            match run {
                Some((first, last)) if residue.position == last + 1 => {
                    run = Some((first, residue.position));
                    continue;
                }
                _ => {
                    push_run(&mut tokens, run.take());
                    run = Some((residue.position, residue.position));
                }
            }
        } else {
            push_run(&mut tokens, run.take());
            tokens.push(format!("{}{}", residue.position, residue.insertion_code));
        }
    }
    push_run(&mut tokens, run);
    tokens.join(" ")
}

fn push_run(tokens: &mut Vec<String>, run: Option<(u16, u16)>) {
    match run {
        Some((first, last)) if first == last => tokens.push(first.to_string()),
        Some((first, last)) => tokens.push(format!("{first}-{last}")),
        None => {}
    }
}

fn assert_e_value_parity(id: &str, observed: f64, expected: f64) {
    assert!(observed.is_finite() && observed > 0.0, "{id} Rust E-value");
    assert!(
        expected.is_finite() && expected > 0.0,
        "{id} ANARCI E-value"
    );
    let relative_difference = (observed - expected).abs() / expected;
    assert!(
        relative_difference <= 0.08,
        "{id} E-value: Rust {observed:e}, ANARCI {expected:e}, relative difference {relative_difference}"
    );
}

fn assert_hit_parity(id: &str, hit_index: usize, observed: &ProfileHit, expected: &ReferenceHit) {
    let label = format!("{id} alternative hit {hit_index}");
    assert_eq!(observed.profile, expected.profile, "{label} profile");
    assert_eq!(
        format!("{}_{:?}", observed.species, observed.chain_type),
        observed.profile,
        "{label} profile components"
    );
    assert!(
        (observed.bit_score - expected.bit_score).abs() <= 0.2,
        "{label} score: Rust {}, ANARCI {}",
        observed.bit_score,
        expected.bit_score
    );
    let relative_e_value_difference =
        (observed.e_value - expected.e_value).abs() / expected.e_value;
    assert!(
        relative_e_value_difference <= 0.08,
        "{label} E-value: Rust {:e}, ANARCI {:e}, relative difference {relative_e_value_difference}",
        observed.e_value,
        expected.e_value
    );
    assert!(
        (observed.bias - expected.bias).abs() <= 0.2,
        "{label} bias: Rust {}, ANARCI {}",
        observed.bias,
        expected.bias
    );
    assert_eq!(observed.query_start, expected.query_start, "{label} start");
    assert_eq!(observed.query_end, expected.query_end, "{label} end");
}

fn strict_alternative_reference(id: &str) -> bool {
    matches!(
        id,
        "trastuzumab_vh"
            | "trastuzumab_vl"
            | "human_trav12_2_traj33"
            | "human_trbv19_trbj2_7"
            | "human_trgv9_trgjp"
            | "human_trdv2_trdj1"
    )
}

fn assert_germline_parity(
    id: &str,
    observed: Option<&GermlineAssignment>,
    expected: Option<&ReferenceGermline>,
) {
    match (observed, expected) {
        (None, None) => {}
        (Some(observed), Some(expected)) => {
            assert_eq!(observed.species, expected.species, "{id} germline species");
            assert_eq!(observed.v_gene, expected.v_gene, "{id} V gene");
            assert_eq!(observed.j_gene, expected.j_gene, "{id} J gene");
            assert_optional_identity(id, "V", observed.v_identity, expected.v_identity);
            assert_optional_identity(id, "J", observed.j_identity, expected.j_identity);
        }
        _ => panic!("{id} germline presence differs"),
    }
}

fn assert_optional_identity(id: &str, segment: &str, observed: Option<f32>, expected: Option<f32>) {
    match (observed, expected) {
        (None, None) => {}
        (Some(observed), Some(expected)) => assert!(
            (observed - expected).abs() <= 0.08,
            "{id} {segment} identity: Rust {observed}, ANARCI {expected}"
        ),
        _ => panic!("{id} {segment} identity presence differs"),
    }
}
