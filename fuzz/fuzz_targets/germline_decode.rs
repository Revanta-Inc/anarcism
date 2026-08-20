#![no_main]

use anarcism_core::GermlineDatabase;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = GermlineDatabase::from_bytes(data);
});
