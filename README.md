<p align="center">
  <img src="demo/public/anarcism.png" alt="anarcism logo" width="300">
</p>

# anarcism

`anarcism` is a self-contained Rust implementation of the ANARCI behavior needed to recognize and IMGT-number immunoglobulin and T-cell-receptor variable domains, with WebAssembly/browser and Python interfaces. It runs locally, bundles every H/K/L/A/B/G/D profile in ANARCI 2026.2.13.2, supports multiple domains and optional V/J assignment, and makes no runtime network request for models or sequence analysis.

The complete counted browser package is 314,492 bytes with gzip and 265,996 bytes with Brotli, below the 500,000-byte limits. See [size and performance measurements](docs/benchmarks.md) for the exact reproducible report.

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

For lower browser batch latency, the optional asynchronous worker-pool entry point creates one independently initialized WASM engine per worker. It balances records by sequence length and restores input order before resolving:

```js
import { AnarcismWorkerPool } from "@revanta/anarcism/worker-pool";

const pool = new AnarcismWorkerPool();
await pool.initialize(4);
const batch = await pool.numberSequences(inputs, { minBitScore: 80 });
await pool.resize(8);
const fastaResults = await pool.numberFasta(fasta);
pool.terminate();
```

The synchronous browser API remains single-threaded and requires no WASM threads or shared memory. Native Rust callers can opt into multicore batches with `number_sequences_parallel(inputs, options, worker_count)`. Both parallel APIs take an explicit worker count to avoid oversubscribing applications that already manage parallel work, and both preserve input order.

## Python use

The Python wheel contains the same embedded Rust engine and has no runtime
dependency on ANARCI, HMMER, or model files:

```python
import anarcism

result = anarcism.number_sequence(vh, assign_germline=True)
batch = anarcism.number_sequences(
    [("heavy", vh), ("light", vl)],
    workers=4,
)
```

It installs the `ANARCI` compatibility command plus the `anarcism`
package-native spelling; `python -m anarcism` is equivalent. The command uses
ANARCI's vertical output by default and supports its IMGT-capable options:

```sh
ANARCI -i EVQLVESGGGLVQPGGSLRLSC...
ANARCI -i sequences.fasta -o numbered.anarci -p 8
ANARCI -i sequences.fasta -o numbered --csv --assign_germline
ANARCI -i sequences.fasta -r ig --use_species human
```

`-i/--sequence`, `-o/--outfile`, `-s/--scheme` (`imgt` or `i`),
`-r/--restrict`, `--csv`, `-p/--ncpu`, `--assign_germline`,
`--use_species`, and `--bit_score_threshold` match ANARCI's flag names and
output layout. Scores are rendered to one decimal place and E-values to two
significant digits, matching the precision of ANARCI's HMMER text output.
Other numbering schemes are intentionally rejected. `--hmmerpath` is
inapplicable because the backend is embedded, and `--outfile_hits` is not yet
exposed because the public result does not contain HMMER's full hit-table bias
and coordinate fields.

## Output conventions

- `start` and `end` are zero-based, half-open input coordinates: `[start, end)`.
- `numbering` contains residues only; HMM deletion states are omitted.
- `paddedImgtAlignment` is local to one detected domain. It contains `-` for missing IMGT positions, includes insertion residues beside their anchor, can exceed 128 characters when insertions occur, and may end before 128 for a C-terminal truncation. Input flanks are not included.
- Region annotations use FR1 1–26, CDR1 27–38, FR2 39–55, CDR2 56–65, FR3 66–104, CDR3 105–117, and FR4 118–128.
- `eValue` is HMMER's independent domain E-value (`i-Evalue`), calculated from the profile's Forward `tau`/`lambda`, the null2-corrected domain bit score, and ANARCI's complete 29-profile search space. It is present for the selected domain and every alternative hit. Conditional domain E-values and whole-sequence E-values are not returned.
- There is no confidence field. Alternative profile hits expose the underlying score separation directly.

## Scope and limits

The embedded profile inventory covers human, mouse, rat, rabbit, rhesus, pig, alpaca, and cow wherever the pinned ANARCI release supplies a profile. Input limits are 10,000 residues per sequence, 1,000 records per batch, and 10 MiB per FASTA document. Only the 20 canonical amino-acid symbols are accepted after ASCII whitespace removal and uppercase normalization.

The default bit-score threshold is 80. `allowedSpecies` is deliberately strict: unavailable or nonmatching species do not silently fall back to all profiles as some ANARCI paths do.

## Validation status

The checked-in `corpus_v2` golden corpus contains 1,397 cases and 1,399 accepted domains across IG H/K/L and TR A/B/G/D, including truncation, CDR-length ladders, insertion/deletion, scFv, multidomain, boundary, constant-domain, and negative controls. Against pinned ANARCI/HMMER it currently has:

- 1,399/1,399 exact domain calls, chain/species classifications, boundaries, and padded alignments
- 153,249/153,249 exact residue indices, IMGT positions, and insertion codes
- maximum absolute bit-score difference of 0.1432 bits
- 1,129/1,399 independent E-values within HMMER's two-significant-digit display precision, with 1.16% mean and 6.79% maximum relative difference

The residual score difference comes from high-precision profile quantization and scalar arithmetic; see [compatibility](docs/compatibility.md). Chromium and WebKit pass locally. Firefox is configured in the Playwright/CI matrix, but the downloaded Firefox 153 runner could not launch in the current macOS sandbox; even a blank-page launch stalled before loading this project.

## Repository map

- `crates/anarcism-core`: safe Rust sequence, HMM, numbering, FASTA, pair, and germline engine
- `crates/anarcism-modelgen`: deterministic compact profile and germline generators
- `crates/anarcism-wasm`: minimal JSON-over-C-ABI WASM bridge
- `browser`: JavaScript package, declarations, Node and Playwright tests, benchmarks
- `assets`: embedded model binaries and provenance manifest
- `tests/golden`: pinned reference corpus
- `fuzz`: libFuzzer targets for both decoders and untrusted text
- `demo`: Preact, Vite, and Tailwind demo application, deployed to GitHub Pages

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

### Demo

The demo is a Preact application built with Vite and Tailwind. It imports `browser/dist` through the `@anarcism` alias, so it exercises the published package the way an npm dependent would; `tools/build-browser.sh` has to run first.

```sh
tools/build-browser.sh
cd demo
npm ci
npm run dev          # http://localhost:5173/
npm run build        # type-checks, then writes demo/dist
```

Analysis runs in a configurable pool of module Web Workers so synchronous WASM calls do not block the page. The demo starts one worker for low single-sequence latency, lazily grows the pool for FASTA batches, balances records across those workers, and returns results in input order.

`.github/workflows/pages.yml` builds the WASM package and the demo on every push to `main` and publishes `demo/dist` to GitHub Pages, which requires Pages to be set to the GitHub Actions source once in the repository settings. The build uses a relative `base`, so the same output works under any Pages path or a custom domain.

Further documentation:

- [Feasibility gate](docs/feasibility.md)
- [Algorithm and data formats](docs/algorithm.md)
- [Compatibility and intentional differences](docs/compatibility.md)
- [Model provenance and licensing](docs/model-provenance.md)
- [Security and privacy](docs/security.md)
- [Build and release](docs/build-release.md)
- [Size and benchmark results](docs/benchmarks.md)

Third-party attributions and license text are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The repository's own licensing remains `UNLICENSED`/`LicenseRef-Proprietary`; no additional project license is implied.
