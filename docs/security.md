# Security and privacy

## Local-only execution

The JavaScript loader fetches only the package's static `anarcism.wasm` file when a URL is used. Profile and germline bytes are inside that file. Sequence, FASTA, result, and error payloads cross only the JavaScript/WASM memory boundary in the same page. The library contains no analytics, logging, `fetch`, XHR, socket, filesystem, or backend analysis path.

The browser test records network requests and verifies that numbering causes only fixture, JavaScript, and WASM GETs—not a sequence-bearing request. Applications can still transmit or log the returned `normalizedSequence`; that is under the host application's control.

## Input and allocation limits

Default limits are:

| Input | Limit |
|---|---:|
| Normalized sequence | 10,000 residues |
| Batch | 1,000 sequences |
| FASTA document | 10 MiB |
| Raw JSON ABI request | 12 MiB |

The FASTA byte limit is checked before parsing. Batch count and individual normalized length are checked before HMM allocation. Profile and germline decoders cap format dimensions and use checked field-size arithmetic. Viterbi memory is bounded by the sequence limit and fixed 128-state models.

Only canonical amino acids are accepted. Structured errors never include a sequence prefix, suffix, or full payload. The native-safe Rust core forbids unsafe code. The WASM crate confines unsafe code to allocation/free and request-slice conversion behind the internal C ABI; the public JavaScript wrapper returns every allocation exactly once.

## Failure behavior

Core operations return `Result` and do not panic on normal invalid input. The bridge returns JSON envelopes and converts failures to `AnarcismError`. Model corruption, invalid options, oversized inputs, invalid FASTA, and invalid residues have distinct codes.

The raw C ABI is not a public untrusted interface. Arbitrary direct calls with forged pointers can trap inside the WASM instance; the JavaScript wrapper passes only pointers returned by the module allocator.

## Verification

- Rust property tests feed arbitrary bytes to both compact decoders and arbitrary text to sequence/FASTA entry points.
- libFuzzer targets cover profile decoding, germline decoding, and sequence/FASTA text.
- Golden and browser tests cover malformed residues and biological edge cases.
- Clippy runs with warnings denied; the core has `#![forbid(unsafe_code)]`.

Run fuzzing with a nightly toolchain:

```sh
cd fuzz
cargo +nightly fuzz run profile_decode
cargo +nightly fuzz run germline_decode
cargo +nightly fuzz run sequence_inputs
```

## Threat-model exclusions

The engine is a sequence classifier/numberer, not a clinical diagnostic. It does not authenticate model provenance at runtime because assets are part of the signed/served WASM resource; distributors should apply their normal package integrity and content-security controls. Side-channel resistance against hostile code running in the same JavaScript realm is not claimed.
