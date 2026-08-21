use anarcism_core::{
    DomainResult, NumberingOptions, PairValidationOptions, number_sequence, validate_antibody_pair,
};
use serde::Deserialize;
use std::fmt::Write;

#[derive(Deserialize)]
struct Corpus {
    cases: Vec<Case>,
    pairs: Vec<Pair>,
}

#[derive(Deserialize)]
struct Case {
    id: String,
    sequence: String,
    #[serde(rename = "referenceDomains")]
    reference_domains: Vec<ReferenceDomain>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct ReferenceDomain {
    chain_type: String,
    species: String,
    start: usize,
    end: usize,
    bit_score: f32,
    numbering: Vec<ReferenceResidue>,
    padded_imgt_alignment: String,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct ReferenceResidue {
    sequence_index: usize,
    amino_acid: char,
    position: u16,
    insertion_code: String,
}

#[derive(Deserialize)]
struct Pair {
    id: String,
    vh: String,
    vl: String,
    ok: bool,
}

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
    numbering_length: usize,
    numbering_fnv1a64: String,
    padded_alignment_fnv1a64: String,
}

fn corpus() -> Corpus {
    serde_json::from_str(include_str!("../../../tests/golden/corpus.json"))
        .expect("checked-in golden corpus is valid JSON")
}

#[test]
fn detection_boundaries_and_imgt_labels_match_pinned_anarci() {
    for case in corpus().cases {
        let result = number_sequence(&case.sequence, &NumberingOptions::default())
            .unwrap_or_else(|error| panic!("{} failed: {error}", case.id));
        assert_eq!(
            result.domains.len(),
            case.reference_domains.len(),
            "{} domain count",
            case.id
        );
        for (observed, expected) in result.domains.iter().zip(&case.reference_domains) {
            assert_eq!(
                format!("{:?}", observed.chain_type),
                expected.chain_type,
                "{} chain",
                case.id
            );
            assert_eq!(observed.species, expected.species, "{} species", case.id);
            assert_eq!(observed.start, expected.start, "{} start", case.id);
            assert_eq!(observed.end, expected.end, "{} end", case.id);
            assert!(
                (observed.bit_score - expected.bit_score).abs() <= 4.0,
                "{} score: Rust {}, ANARCI {}",
                case.id,
                observed.bit_score,
                expected.bit_score
            );
            assert_eq!(
                observed.padded_imgt_alignment, expected.padded_imgt_alignment,
                "{} padded alignment",
                case.id
            );
            assert_eq!(
                observed.numbering.len(),
                expected.numbering.len(),
                "{} residue count",
                case.id
            );
            for (observed_residue, expected_residue) in
                observed.numbering.iter().zip(&expected.numbering)
            {
                assert_eq!(
                    (
                        observed_residue.sequence_index,
                        observed_residue.amino_acid,
                        observed_residue.position,
                        observed_residue.insertion_code.as_str(),
                    ),
                    (
                        expected_residue.sequence_index,
                        expected_residue.amino_acid,
                        expected_residue.position,
                        expected_residue.insertion_code.as_str(),
                    ),
                    "{} residue parity",
                    case.id
                );
            }
        }
    }
}

#[test]
fn antibody_pair_outcomes_match_the_golden_cases() {
    for pair in corpus().pairs {
        let observed =
            validate_antibody_pair(&pair.vh, &pair.vl, &PairValidationOptions::default())
                .unwrap_or_else(|error| panic!("{} failed: {error}", pair.id));
        assert_eq!(observed.ok, pair.ok, "{} pair outcome", pair.id);
    }
}

#[test]
fn expanded_corpus_matches_pinned_anarci() {
    let inputs: Vec<V2InputCase> =
        serde_json::from_str(include_str!("../../../tests/golden/corpus_v2.json"))
            .expect("checked-in expanded corpus is valid JSON");
    let reference: V2Reference = serde_json::from_str(include_str!(
        "../../../tests/golden/corpus_v2_reference.json"
    ))
    .expect("checked-in expanded reference is valid JSON");
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
                for case_index in (shard_index..inputs.len()).step_by(worker_count) {
                    assert_v2_case(&inputs[case_index], &reference_cases[case_index]);
                }
            }));
        }
        for worker in workers {
            worker.join().expect("expanded parity worker panicked");
        }
    });
}

fn assert_v2_case(input: &V2InputCase, expected_case: &V2ReferenceCase) {
    assert_eq!(input.id, expected_case.id);
    assert_eq!(input.category, expected_case.category);
    let result = number_sequence(&input.seq, &NumberingOptions::default())
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
            format!("{:?}", observed.chain_type),
            expected.chain_type,
            "{} chain",
            input.id
        );
        assert_eq!(observed.species, expected.species, "{} species", input.id);
        assert_eq!(observed.start, expected.start, "{} start", input.id);
        assert_eq!(observed.end, expected.end, "{} end", input.id);
        assert!(
            (observed.bit_score - expected.bit_score).abs() <= 4.0,
            "{} score: Rust {}, ANARCI {}",
            input.id,
            observed.bit_score,
            expected.bit_score
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
