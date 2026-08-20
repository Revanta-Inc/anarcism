# anarcism

`anarcism` is a browser-native Rust/WebAssembly implementation of the ANARCI behavior needed to recognize and IMGT-number immunoglobulin and T-cell-receptor variable domains. It runs locally, bundles every H/K/L/A/B/G/D profile in ANARCI 2026.2.13.2, supports multiple domains and optional V/J assignment, and makes no runtime network request for models or sequence analysis.

The complete counted browser package is below 250 kB with gzip and below 210 kB with Brotli. See [size and performance measurements](docs/benchmarks.md) for the exact reproducible report.

## Browser use

Build the package, then initialize the WASM module before using the synchronous numbering calls:

```sh
tools/build-browser.sh
```

```js
import init, {
  numberSequence,
  numberSequences,
  numberFasta,
  validateAntibodyPair,
} from "@revanta/anarcism";

await init();

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

## Output conventions

- `start` and `end` are zero-based, half-open input coordinates: `[start, end)`.
- `numbering` contains residues only; HMM deletion states are omitted.
- `paddedImgtAlignment` is local to one detected domain. It contains `-` for missing IMGT positions, includes insertion residues beside their anchor, can exceed 128 characters when insertions occur, and may end before 128 for a C-terminal truncation. Input flanks are not included.
- Region annotations use FR1 1–26, CDR1 27–38, FR2 39–55, CDR2 56–65, FR3 66–104, CDR3 105–117, and FR4 118–128.
- `eValue` is absent in this release because exact HMMER posterior null2 correction is not implemented. A misleading approximate E-value is not returned.
- There is no confidence field. Alternative profile hits expose the underlying score separation directly.

## Scope and limits

The embedded profile inventory covers human, mouse, rat, rabbit, rhesus, pig, alpaca, and cow wherever the pinned ANARCI release supplies a profile. Input limits are 10,000 residues per sequence, 1,000 records per batch, and 10 MiB per FASTA document. Only the 20 canonical amino-acid symbols are accepted after ASCII whitespace removal and uppercase normalization.

The default bit-score threshold is 80. `allowedSpecies` is deliberately strict: unavailable or nonmatching species do not silently fall back to all profiles as some ANARCI paths do.

## Validation status

The checked-in golden corpus has 15 accepted domains across IG H/K/L and TR A/B/G/D, plus truncation, insertion, deletion, scFv, boundary, swapped-pair, and negative controls. Against pinned ANARCI/HMMER it currently has:

- 15/15 exact domain calls, chain/species classifications, boundaries, and padded alignments
- 1,719/1,719 exact residue indices, IMGT positions, and insertion codes
- maximum absolute bit-score difference of 3.3592 bits

The score difference is caused by compact trace-based null2 correction; see [compatibility](docs/compatibility.md). Chromium and WebKit pass locally. Firefox is configured in the Playwright/CI matrix, but the downloaded Firefox 153 runner could not launch in the current macOS sandbox; even a blank-page launch stalled before loading this project.

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
```

Run `node tools/serve-demo.mjs` from the repository root and open <http://127.0.0.1:8080/> for the demo.

Further documentation:

- [Feasibility gate](docs/feasibility.md)
- [Algorithm and data formats](docs/algorithm.md)
- [Compatibility and intentional differences](docs/compatibility.md)
- [Model provenance and licensing](docs/model-provenance.md)
- [Security and privacy](docs/security.md)
- [Build and release](docs/build-release.md)
- [Size and benchmark results](docs/benchmarks.md)

Third-party attributions and license text are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The repository's own licensing remains `UNLICENSED`/`LicenseRef-Proprietary`; no additional project license is implied.
