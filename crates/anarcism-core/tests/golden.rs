use anarcism_core::{
    DomainResult, GermlineAssignment, NumberingOptions, ProfileHit, number_sequence,
};
#[cfg(not(target_family = "wasm"))]
use anarcism_core::{SequenceInput, number_sequences, number_sequences_parallel};
use serde::Deserialize;
use std::fmt::Write;
#[cfg(not(target_family = "wasm"))]
use std::num::NonZeroUsize;

#[derive(Deserialize)]
struct V2InputCase {
    id: String,
    seq: String,
    category: String,
}

#[derive(Deserialize)]
struct V2Reference {
    cases: Vec<V2ReferenceCase>,
}

#[derive(Deserialize)]
struct V2ReferenceCase {
    id: String,
    category: String,
    domains: Vec<V2ReferenceDomain>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct V2ReferenceDomain {
    chain_type: String,
    species: String,
    start: usize,
    end: usize,
    bit_score: f32,
    e_value: f64,
    e_value_display_significant_digits: u32,
    score_components: V2ReferenceScoreComponents,
    alternative_hits: Vec<V2ReferenceHit>,
    germline: Option<V2ReferenceGermline>,
    numbering_length: usize,
    numbering_fnv1a64: String,
    padded_alignment_fnv1a64: String,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct V2ReferenceHit {
    profile: String,
    chain_type: String,
    species: String,
    bit_score: f32,
    e_value: f64,
    e_value_display_significant_digits: u32,
    bias: f32,
    query_start: usize,
    query_end: usize,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct V2ReferenceGermline {
    species: String,
    v_gene: Option<String>,
    v_identity: Option<f32>,
    j_gene: Option<String>,
    j_identity: Option<f32>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct V2ReferenceScoreComponents {
    envelope_start: usize,
    envelope_end: usize,
    envelope_forward_nats: f32,
    outside_envelope_nats: f32,
    null_one_nats: f32,
    null2_bias_bits: f32,
    display_precision_bits: f32,
}

#[test]
fn corpus_v2_matches_versioned_anarci_reference() {
    let inputs: Vec<V2InputCase> =
        serde_json::from_str(include_str!("../../../tests/golden/corpus_v2.json"))
            .expect("checked-in corpus_v2 input is valid JSON");
    let reference: V2Reference = serde_json::from_str(include_str!(
        "../../../tests/golden/corpus_v2_reference.json"
    ))
    .expect("checked-in corpus_v2 reference is valid JSON");
    assert_eq!(inputs.len(), reference.cases.len());

    let worker_count = std::thread::available_parallelism()
        .map_or(1, std::num::NonZero::get)
        .min(inputs.len().max(1));
    std::thread::scope(|scope| {
        let mut workers = Vec::with_capacity(worker_count);
        for shard_index in 0..worker_count {
            let inputs = &inputs;
            let reference_cases = &reference.cases;
            workers.push(scope.spawn(move || {
                let mut failures = Vec::new();
                for case_index in (shard_index..inputs.len()).step_by(worker_count) {
                    if std::panic::catch_unwind(|| {
                        assert_v2_case(&inputs[case_index], &reference_cases[case_index]);
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
            "{} corpus_v2 cases failed parity: {}",
            failures.len(),
            failures.join(", ")
        );
    });
}

#[test]
#[cfg(not(target_family = "wasm"))]
fn full_native_batch_scheduler_matches_the_serial_batch() {
    let cases: Vec<V2InputCase> =
        serde_json::from_str(include_str!("../../../tests/golden/corpus_v2.json"))
            .expect("checked-in corpus_v2 input is valid JSON");
    let inputs: Vec<_> = cases
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

fn assert_v2_case(input: &V2InputCase, expected_case: &V2ReferenceCase) {
    assert_eq!(input.id, expected_case.id);
    assert_eq!(input.category, expected_case.category);
    let result = number_sequence(
        &input.seq,
        &NumberingOptions {
            assign_germline: true,
            alternative_hit_count: 28,
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
        assert_reference_score_components(input, expected);
        assert_eq!(
            format!("{:?}", observed.chain_type),
            expected.chain_type,
            "{} chain",
            input.id
        );
        assert_eq!(observed.species, expected.species, "{} species", input.id);
        assert_eq!(observed.start, expected.start, "{} start", input.id);
        assert_eq!(observed.end, expected.end, "{} end", input.id);
        assert!(
            (observed.bit_score - expected.bit_score).abs() <= 0.2,
            "{} score: Rust {}, ANARCI {}",
            input.id,
            observed.bit_score,
            expected.bit_score
        );
        assert_e_value_parity(
            &input.id,
            observed.e_value,
            expected.e_value,
            expected.e_value_display_significant_digits,
        );
        assert!(
            (observed.bias - expected.score_components.null2_bias_bits).abs() <= 0.2,
            "{} bias: Rust {}, ANARCI {}",
            input.id,
            observed.bias,
            expected.score_components.null2_bias_bits
        );
        let mut matched_alternatives = 0;
        for (hit_index, expected_hit) in expected.alternative_hits.iter().enumerate() {
            if let Some(actual_hit) = observed
                .alternative_hits
                .iter()
                .find(|hit| hit.profile == expected_hit.profile)
            {
                matched_alternatives += 1;
                assert_eq!(
                    format!("{:?}", actual_hit.chain_type),
                    expected_hit.chain_type,
                    "{} alternative hit {hit_index} chain",
                    input.id
                );
                assert_eq!(
                    actual_hit.species, expected_hit.species,
                    "{} alternative hit {hit_index} species",
                    input.id
                );
            }
        }
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
            observed.numbering.len(),
            expected.numbering_length,
            "{} residue count",
            input.id
        );
        assert_eq!(
            numbering_fingerprint(observed),
            expected.numbering_fnv1a64,
            "{} residue parity",
            input.id
        );
        assert_eq!(
            fnv1a64(observed.padded_imgt_alignment.as_bytes()),
            expected.padded_alignment_fnv1a64,
            "{} padded alignment",
            input.id
        );
    }
}

fn assert_e_value_parity(id: &str, observed: f64, expected: f64, significant_digits: u32) {
    assert_eq!(significant_digits, 2, "{id} reference E-value precision");
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

fn assert_hit_parity(id: &str, hit_index: usize, observed: &ProfileHit, expected: &V2ReferenceHit) {
    let label = format!("{id} alternative hit {hit_index}");
    assert_eq!(observed.profile, expected.profile, "{label} profile");
    assert_eq!(
        format!("{:?}", observed.chain_type),
        expected.chain_type,
        "{label} chain"
    );
    assert_eq!(observed.species, expected.species, "{label} species");
    assert!(
        (observed.bit_score - expected.bit_score).abs() <= 0.2,
        "{label} score: Rust {}, ANARCI {}",
        observed.bit_score,
        expected.bit_score
    );
    assert_eq!(
        expected.e_value_display_significant_digits, 2,
        "{label} reference E-value precision"
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
    expected: Option<&V2ReferenceGermline>,
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

fn assert_reference_score_components(input: &V2InputCase, expected: &V2ReferenceDomain) {
    let components = &expected.score_components;
    assert!(
        components.envelope_start < components.envelope_end,
        "{} reference envelope",
        input.id
    );
    assert!(
        components.envelope_end <= input.seq.len(),
        "{} reference envelope end",
        input.id
    );
    assert_eq!(components.display_precision_bits, 0.1);
    let recomposed = (components.envelope_forward_nats + components.outside_envelope_nats
        - components.null_one_nats)
        / std::f32::consts::LN_2
        - components.null2_bias_bits;
    assert!(
        (recomposed - expected.bit_score).abs() <= 1.0e-4,
        "{} reference score components recompose to {recomposed}, expected {}",
        input.id,
        expected.bit_score
    );
}

fn numbering_fingerprint(domain: &DomainResult) -> String {
    let mut canonical = String::new();
    for residue in &domain.numbering {
        writeln!(
            canonical,
            "{}|{}|{}|{}",
            residue.sequence_index, residue.amino_acid, residue.position, residue.insertion_code
        )
        .expect("writing to a String is infallible");
    }
    fnv1a64(canonical.as_bytes())
}

fn fnv1a64(bytes: &[u8]) -> String {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    format!("{hash:016x}")
}
