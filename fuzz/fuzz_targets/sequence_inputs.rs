#![no_main]

use anarcism_core::{NumberingOptions, number_fasta, number_sequence};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(text) = std::str::from_utf8(data) {
        let options = NumberingOptions::default();
        let _ = number_sequence(text, &options);
        let _ = number_fasta(text, &options);
    }
});
