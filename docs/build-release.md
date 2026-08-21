# Build and release

## Toolchain

The workspace uses Rust edition 2024 and pins Rust 1.94.0 in `rust-toolchain.toml`. The measured release used:

- rustc 1.94.0
- Node.js 24.13.1 (package minimum: Node 20 for development tools)
- Binaryen/`wasm-opt` 132
- Playwright 1.62.1
- TypeScript 7.0.2

Runtime users need a modern browser with WebAssembly SIMD128, `i64`/BigInt integration, `TextEncoder`, and `TextDecoder`.

## Reproducible browser build

From a clean checkout:

```sh
cargo test --workspace --locked
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check

cd browser
npm ci
cd ..
tools/build-browser.sh
```

The script builds the `wasm32-unknown-unknown` release with `simd128`, `opt-level="z"`, LTO, one codegen unit, aborting panics, and stripped symbols. It then runs `wasm-opt -Oz` with SIMD and the other modern browser features emitted by Rust, copies the JavaScript/declarations/third-party notice, and fails if either gzip or Brotli reaches 500,000 bytes.

The output directory is `browser/dist`. It contains no source profile sidecar and makes no model download; `profiles.bin` and `germlines.bin` are linked into the WASM data segment.

Cargo and npm lockfiles pin all build/test dependencies. Source order, floating-point quantization, binary encoding, JSON generation, and profile tie breaks are deterministic. The final hash can differ across compiler/Binaryen versions, which is why release tooling versions are recorded above and in CI.

## Regenerating model assets

This is necessary only when deliberately changing the ANARCI pin. The checked-in runtime build does not need Python or HMMER.

```sh
uv sync --locked --no-install-project --no-default-groups --group reference

cargo run --locked -p anarcism-modelgen --bin anarcism-modelgen -- \
  .venv/lib/python3.14/site-packages/anarci/dat/HMMs/ALL.hmm \
  assets/profiles.bin

.venv/bin/python tools/export_germlines.py /tmp/anarcism-germlines.tsv
cargo run --locked -p anarcism-modelgen --bin germlinegen -- \
  /tmp/anarcism-germlines.tsv assets/germlines.bin
```

Then update and verify `assets/MANIFEST.toml`, inspect profile/germline inventory changes, rerun every test, and redo the license review. Both generators fail on an unexpected version, alphabet, model length, symbol, or metadata shape.

## Regenerating golden data

Install HMMER 3.4 on `PATH`, use the locked Python environment, and run:

```sh
.venv/bin/python tools/generate_golden.py tests/golden/corpus.json
cargo test --locked -p anarcism-core --test golden
tools/build-browser.sh
node tools/parity-report.mjs
```

Golden JSON is reviewed and committed; it is not silently rewritten in normal CI.

## Test matrix

```sh
cargo test --workspace --locked
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check

cd browser
npm test
npm run test:types
npm run parity
npm run test:browsers
npm run benchmark:browser
```

The Playwright configuration runs Chromium, Firefox, and WebKit. CI installs all three official runtimes. On the development macOS sandbox used for the recorded measurements, Chromium and WebKit pass; Playwright's downloaded Firefox 153 binary stalls during a minimal blank-page `firefox.launch()` before any repository page or WASM is loaded. That environmental exception is not converted into a skipped CI project.

Scheduled CI fuzzes all three targets with libFuzzer. Local fuzz runs should be longer before a release that changes parsing or DP indexing.

## Packaging

After a clean successful build:

```sh
cd browser
npm pack --dry-run
npm pack
```

Confirm that the tarball contains `dist/anarcism.wasm`, `dist/index.js`, `dist/index.d.ts`, and `dist/THIRD_PARTY_NOTICES.md`. The npm package is marked `UNLICENSED` because the repository owner has not selected a project license; upstream notices do not license original project code.

### Python wheels

The root `pyproject.toml` builds the `anarcism` wheel from the `crates/anarcism-python` bindings crate. The build is driven from the repository root because the sdist has to carry the whole Cargo workspace.

```sh
uvx maturin build --release   # honors [tool.maturin] profile = "python"
uvx maturin sdist
```

Three settings are load-bearing and are verified by the release checklist rather than by a test:

- `profile = "python"` selects `[profile.python]` from the workspace `Cargo.toml`. The shipping `[profile.release]` sets `panic = "abort"` for the browser byte budget, and under that profile any Rust panic aborts the host interpreter with `SIGABRT` instead of raising. `[profile.python]` overrides `panic = "unwind"` so PyO3 converts a panic into `PanicException`, and restores `opt-level = 3` because a wheel has no size budget.
- `include` names `assets/**/*`, `rust-toolchain.toml`, and `THIRD_PARTY_NOTICES.md` for the sdist. `cargo package --list` only reports files inside the selected crate, so without these entries the sdist builds a wheel that fails on the `include_bytes!` of `assets/profiles.bin`.
- abi3 produces one wheel per platform covering every supported CPython, rather than one per interpreter version. abi3 does not cover free-threaded builds; those need their own non-abi3 wheels.

Verify the sdist in a clean environment, since a wheel built in the source tree can succeed while the sdist cannot:

```sh
uvx maturin sdist
tar tzf target/wheels/anarcism-*.tar.gz   # expect assets/ and the workspace Cargo.toml
pip install --no-binary :all: target/wheels/anarcism-*.tar.gz
```

The wheel is self-contained: both model binaries are linked into the extension, so there is no runtime HMMER, MUSCLE, ANARCI, or data-file dependency. `THIRD_PARTY_NOTICES.md` is declared through `license-files` and lands in `anarcism-*.dist-info/licenses/`.

The release checklist is:

1. All Rust, property, golden, Node, type, and three-browser CI jobs pass.
2. Fuzz workflow is green.
3. `npm audit --omit=dev` reports no runtime findings.
4. Parity counts and maximum score difference are unchanged or reviewed.
5. Both compressed totals are below 500,000 bytes.
6. Model hashes match the manifest.
7. Third-party notices are present in the tarball.
8. Benchmarks and this documentation reflect the release artifact.
9. The Python sdist installs with `--no-binary :all:` in a clean environment.
10. Wheels are abi3 and one per target platform, and a Rust panic raises `PanicException` rather than aborting the interpreter.
