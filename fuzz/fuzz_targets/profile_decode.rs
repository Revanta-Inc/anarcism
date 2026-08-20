#![no_main]

use anarcism_core::ProfileDatabase;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = ProfileDatabase::from_bytes(data);
});
