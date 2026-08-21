# anarcism

`anarcism` is a browser-native Rust/WebAssembly implementation of the ANARCI behavior needed to recognize and IMGT-number immunoglobulin and T-cell-receptor variable domains. It runs locally, bundles every H/K/L/A/B/G/D profile in ANARCI 2026.2.13.2, supports multiple domains and optional V/J assignment, and makes no runtime network request for models or sequence analysis.

The complete counted browser package is 308,064 bytes with gzip and 260,884 bytes with Brotli, below the 500,000-byte limits. See [size and performance measurements](docs/benchmarks.md) for the exact reproducible report.

## Browser use

Build the package, then initialize the WASM module before using the synchronous numbering calls:

```sh
tools/build-browser.sh
```

The release artifact requires WebAssembly SIMD128 but does not require WASM threads or shared memory.

```js
import init, {
  numberSequence,
  numberSequences,
  numberFasta,
  validateAntibodyPair,
} from "@revanta/anarcism";

const engine = await init();

engine.version;   // package version
engine.chains();  // chain types present in the embedded profiles
engine.species(); // species present in the embedded profiles

const result = numberSequence(vh, {
  allowedChains: ["H", "K", "L"],
  allowedSpecies: ["human", "mouse"],
  minBitScore: 80,
  alternativeHitCount: 3,
  assignGermline: true,
});

const batch = numberSequences([
  { id: "heavy", sequence: vh },
  { id: "light", sequence: vl },
]);

const fastaResults = numberFasta(`>heavy\n${vh}\n`);
const pair = validateAntibodyPair(vh, vl, { startMax: 10, endMin: 100 });
```

`initSync(bytesOrModule)` is also available for applications that already have WASM bytes or a compiled `WebAssembly.Module`. Once initialized, all analysis methods are synchronous. No method sends a sequence anywhere.

Errors are thrown as `AnarcismError` with `code`, `message`, and optional `inputId` fields. TypeScript declarations are included in the package.

Native Rust callers can opt into multicore batches with `number_sequences_parallel(inputs, options, worker_count)`. The explicit `NonZeroUsize` worker count avoids oversubscribing applications that already manage their own thread pools. Results remain in input order; the browser API stays single-threaded and requires no WASM threads or shared memory.

## Output conventions

- `start` and `end` are zero-based, half-open input coordinates: `[start, end)`.
- `numbering` contains residues only; HMM deletion states are omitted.
- `paddedImgtAlignment` is local to one detected domain. It contains `-` for missing IMGT positions, includes insertion residues beside their anchor, can exceed 128 characters when insertions occur, and may end before 128 for a C-terminal truncation. Input flanks are not included.
- Region annotations use FR1 1–26, CDR1 27–38, FR2 39–55, CDR2 56–65, FR3 66–104, CDR3 105–117, and FR4 118–128.
- `eValue` is absent in this release because fixed-point profile scores and scalar arithmetic are not bit-identical to optimized HMMER, and its complete final E-value accounting is outside the supported subset. A misleading approximate E-value is not returned.
- There is no confidence field. Alternative profile hits expose the underlying score separation directly.

## Scope and limits

The embedded profile inventory covers human, mouse, rat, rabbit, rhesus, pig, alpaca, and cow wherever the pinned ANARCI release supplies a profile. Input limits are 10,000 residues per sequence, 1,000 records per batch, and 10 MiB per FASTA document. Only the 20 canonical amino-acid symbols are accepted after ASCII whitespace removal and uppercase normalization.

The default bit-score threshold is 80. `allowedSpecies` is deliberately strict: unavailable or nonmatching species do not silently fall back to all profiles as some ANARCI paths do.

## Validation status

The checked-in golden corpora contain 1,412 cases and 1,414 accepted domains across IG H/K/L and TR A/B/G/D, including truncation, CDR-length ladders, insertion/deletion, scFv, multidomain, boundary, swapped-pair, constant-domain, and negative controls. Against pinned ANARCI/HMMER they currently have:

- 1,414/1,414 exact domain calls, chain/species classifications, boundaries, and padded alignments
- 154,968/154,968 exact residue indices, IMGT positions, and insertion codes
- maximum absolute bit-score difference of 0.8449 bits

The residual score difference comes from high-precision profile quantization and scalar arithmetic; see [compatibility](docs/compatibility.md). Chromium and WebKit pass locally. Firefox is configured in the Playwright/CI matrix, but the downloaded Firefox 153 runner could not launch in the current macOS sandbox; even a blank-page launch stalled before loading this project.

## Repository map

- `crates/anarcism-core`: safe Rust sequence, HMM, numbering, FASTA, pair, and germline engine
- `crates/anarcism-modelgen`: deterministic compact profile and germline generators
- `crates/anarcism-wasm`: minimal JSON-over-C-ABI WASM bridge
- `browser`: JavaScript package, declarations, Node and Playwright tests, benchmarks
- `assets`: embedded model binaries and provenance manifest
- `tests/golden`: pinned reference corpus
- `fuzz`: libFuzzer targets for both decoders and untrusted text
- `demo`: local browser demo

## Development

```sh
cargo test --workspace --locked
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check
tools/build-browser.sh
cd browser
npm ci
npm test
npm run test:types
npm run parity
npm run test:browsers
npm run benchmark:compare
```

The comparison benchmark requires `uv sync --locked --no-install-project --no-default-groups --group reference` and HMMER 3.4 on `PATH`; ordinary Rust and browser builds do not.

Run `node tools/serve-demo.mjs` from the repository root and open <http://127.0.0.1:8080/> for the demo. Analysis runs in a dedicated module Web Worker so synchronous WASM calls do not block the page.

Further documentation:

- [Feasibility gate](docs/feasibility.md)
- [Algorithm and data formats](docs/algorithm.md)
- [Compatibility and intentional differences](docs/compatibility.md)
- [Model provenance and licensing](docs/model-provenance.md)
- [Security and privacy](docs/security.md)
- [Build and release](docs/build-release.md)
- [Size and benchmark results](docs/benchmarks.md)

Third-party attributions and license text are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The repository's own licensing remains `UNLICENSED`/`LicenseRef-Proprietary`; no additional project license is implied.
