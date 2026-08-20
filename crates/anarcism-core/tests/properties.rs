use anarcism_core::{
    GermlineDatabase, NumberingOptions, ProfileDatabase, number_fasta, number_sequence,
};
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn compact_decoders_reject_or_parse_arbitrary_bytes_without_panicking(
        bytes in prop::collection::vec(any::<u8>(), 0..2_048),
    ) {
        let _ = ProfileDatabase::from_bytes(&bytes);
        let _ = GermlineDatabase::from_bytes(&bytes);
    }

    #[test]
    fn untrusted_text_returns_a_result_without_panicking(
        bytes in prop::collection::vec(any::<u8>(), 0..512),
    ) {
        let text = String::from_utf8_lossy(&bytes);
        let _ = number_sequence(&text, &NumberingOptions::default());
        let _ = number_fasta(&text, &NumberingOptions::default());
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    #[test]
    fn canonical_sequence_results_are_deterministic(
        residues in prop::collection::vec(
            prop::sample::select("ACDEFGHIKLMNPQRSTVWY".chars().collect::<Vec<_>>()),
            1..180,
        ),
    ) {
        let sequence: String = residues.into_iter().collect();
        let first = number_sequence(&sequence, &NumberingOptions::default());
        let second = number_sequence(&sequence, &NumberingOptions::default());
        prop_assert_eq!(first, second);
    }
}
