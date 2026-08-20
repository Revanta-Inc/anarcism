use anarcism_core::{
    NumberingOptions, PairValidationOptions, number_sequence, validate_antibody_pair,
};
use serde::Deserialize;

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
