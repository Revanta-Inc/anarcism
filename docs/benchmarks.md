# Size and performance results

Measurements below are from the current release artifact built with rustc 1.94.0 and Binaryen 132. Compression is deterministic Node zlib gzip level 9 and Brotli quality 11. The hard threshold is strictly less than 500,000 bytes for both totals.

## Complete distribution

The notice-inclusive final report is produced by `node tools/size-report.mjs`:

| File | Raw | gzip | Brotli |
|---|---:|---:|---:|
| `anarcism.wasm` | 752,231 | 303,339 | 257,032 |
| `index.js` | 4,420 | 1,434 | 1,237 |
| `index.d.ts` | 2,853 | 842 | 740 |
| `THIRD_PARTY_NOTICES.md` | 5,288 | 2,005 | 1,491 |
| `package.json` | 958 | 444 | 384 |
| **Complete distribution** | **765,750** | **308,064** | **260,884** |

This is 61.6% of the gzip limit and 52.2% of the Brotli limit. Run the report after any source change and treat its output as authoritative.

The report counts:

- optimized `anarcism.wasm`, including both embedded assets
- JavaScript glue
- TypeScript declarations
- required third-party notice
- npm package metadata

It does not add standalone asset sizes to the package total because that would double-count bytes already inside WASM.

Standalone asset attribution from the same report:

| Embedded input | Raw | gzip | Brotli |
|---|---:|---:|---:|
| `profiles.bin` | 315,721 | 152,430 | 132,231 |
| `germlines.bin` | 226,705 | 44,108 | 36,195 |

Before Binaryen, the Rust release WASM was 824,987 bytes raw / 294,882 gzip / 249,949 Brotli. `wasm-opt -Oz` reduced raw WASM to 752,231 bytes but increased gzip/Brotli to 303,339/257,032 bytes. The optimized artifact is retained for the release because `-Oz` materially reduces raw transfer/storage and is the specified release pipeline; both compressed totals remain below budget.

## Browser performance

Measured in Playwright Headless Chromium 151.0.7922.34 on an Apple-Silicon arm64 laptop running macOS 27.0, with Node 24.13.1 driving the fixture. The benchmark reads local cached files and reports no network decompression time. Cold initialization includes local fetch, array-buffer creation, compilation, and instantiation; cached initialization receives an already compiled `WebAssembly.Module`.

| Measure | Actual | Target | Status |
|---|---:|---:|---|
| Cold initialization | 6.9 ms | <100 ms cached | Pass even cold |
| Cached initialization | 0.2 ms | <100 ms | Pass |
| Single VH median / p95 | 55.2 / 55.3 ms | <50 ms | Miss |
| Single VL median / p95 | 49.4 / 49.5 ms | <50 ms | Pass |
| 100 VH/VL pairs (200 sequences) | 10,472.0 ms | <1,000 ms | Miss |

The single-sequence and batch throughput targets remain unmet, although replacing per-cell transcendental logsum operations with HMMER's table-driven approximation cut all three measurements by roughly 68%. Analysis is deterministic and single-threaded. The release enables WASM SIMD128 and LLVM emits vector instructions for suitable operations, but the Plan7 recurrences are not hand-vectorized and batching does not process sequences in parallel. Meeting the targets would require explicitly SIMD/batched DP kernels or parallel workers; the shortfall is not hidden with result caching or excluded work.

The Node/V8 cross-check produced 16.23 ms cold initialization, 0.118 ms cached initialization, 56.58/56.77 ms median/p95 VH, 50.85/51.46 ms median/p95 VL, and 10,803.3 ms for 100 pairs. A fresh-process first VH call took 69.5 ms because it also initializes the 64 KiB logsum table; subsequent calls reuse it.

## Cross-runtime comparison

The same VH and VL inputs were measured on August 20, 2026 on an Apple M1 Max with 10 logical CPUs. Each runtime used default numbering options, an IMGT scheme, an 80-bit threshold, no germline assignment, and one worker/thread. After three VH/VL warmups, the harness collected 25 hot single-sequence samples and timed one call containing 100 repeated VH/VL pairs (200 sequences).

| Runtime | VH median / p95 | VL median / p95 | 200-sequence batch | ms/sequence | sequences/s | Batch time vs Rust |
|---|---:|---:|---:|---:|---:|---:|
| Native Rust (`-O3`) | 25.7 / 25.9 ms | 23.0 / 23.1 ms | 4,881.8 ms | 24.409 | 40.968 | 1.00x |
| Browser WASM | 55.2 / 55.3 ms | 49.4 / 49.5 ms | 10,472.0 ms | 52.360 | 19.099 | 2.15x |
| ANARCI/HMMER | 26.1 / 26.3 ms | 22.0 / 22.2 ms | 3,602.5 ms | 18.013 | 55.517 | 0.74x |

Native Rust uses the dedicated `benchmark` Cargo profile with `opt-level = 3`. Browser WASM uses the shipping size-optimized `-Oz` artifact with SIMD128 enabled and calls the package API directly in Playwright Headless Chromium 151, excluding Web Worker messaging, DOM rendering, and initialization. ANARCI 2026.2.13.2 runs under Python 3.14.5 with HMMER 3.4 and `ncpu=1`. Process startup, imports, and WASM initialization are outside the hot timings; each timed ANARCI API call still includes its inherent temporary-file I/O and `hmmscan` subprocess.

The batch result is a public-API throughput comparison, not an equivalent DP-kernel microbenchmark. ANARCI submits all 200 queries to one HMMER scan, while the serial Rust comparison and WASM `numberSequences` implementation loop over 200 independent sequence calls. That amortization makes ANARCI 1.36x faster than serial native Rust for this repeated batch even though their hot single-sequence medians are close. The browser result measures direct single-threaded WASM computation and therefore does not include the demo's responsiveness benefit from moving work to a Web Worker.

### Native batch scaling

Native callers can use the opt-in `number_sequences_parallel` API with an explicit nonzero worker count. It uses scoped standard-library threads, retains input order and exact serial results, adds no dependency, and is omitted from WASM builds. The existing `number_sequences` API remains serial so library callers are never surprised by automatic CPU saturation.

The same 200-sequence workload, measured as the median of three complete batches under the native `-O3` profile, scales as follows:

| Workers | Batch total | ms/sequence | sequences/s | Speedup |
|---:|---:|---:|---:|---:|
| 1 | 4,892.9 ms | 24.465 | 40.875 | 1.00x |
| 2 | 2,456.7 ms | 12.283 | 81.412 | 1.99x |
| 4 | 1,231.0 ms | 6.155 | 162.467 | 3.97x |
| 8 | 662.9 ms | 3.315 | 301.689 | 7.38x |
| 10 | 644.9 ms | 3.225 | 310.102 | 7.59x |

Eight workers bring the native batch below the 1,000 ms target. The small gain from 8 to 10 workers is consistent with this machine's eight performance cores plus two efficiency cores; callers should benchmark an appropriate cap rather than assuming every logical CPU improves latency.

### Optimization attribution

An immediate same-harness Node control separated the effect of enabling `simd128` from the algorithmic change:

| Build | VH median | VL median | 100 pairs |
|---|---:|---:|---:|
| Previous scalar release | 183.46 ms | 161.93 ms | 35,198.3 ms |
| SIMD128 only | 184.15 ms | 162.48 ms | 34,744.1 ms |
| SIMD128 + HMMER logsum table + inverse score scale | 56.58 ms | 50.85 ms | 10,803.3 ms |

SIMD alone was effectively neutral for single sequences; LLVM emitted vector instructions, but the dependent Plan7 recurrences remain scalar. The approximately 68% reduction comes primarily from replacing per-cell `exp`/`ln1p` operations with HMMER's table lookup. The SIMD-only optimized WASM was 748,016 bytes versus 752,616 bytes for the scalar control.

## Running measurements

```sh
uv sync --locked --no-install-project --no-default-groups --group reference
tools/build-browser.sh

cargo run --profile benchmark -p anarcism-core --example native_batch_benchmark

cd browser
npm run size
npm run benchmark
npm run benchmark:browser
npm run benchmark:compare
```

The cross-runtime comparison additionally requires HMMER 3.4 on `PATH`; the runner rejects a different ANARCI or HMMER version. The browser benchmark starts a loopback-only static server, loads the same release files as the demo, warms each path three times, measures 25 single-domain iterations, and includes one complete 200-sequence batch call. `benchmark:compare` runs native Rust, browser WASM, and the pinned `.venv` ANARCI installation sequentially so they do not compete for CPU. Set `ANARCISM_BENCH_WARMUP`, `ANARCISM_BENCH_SINGLE_ITERATIONS`, or `ANARCISM_BENCH_PAIR_COUNT` to positive integers to override the shared workload. Both commands print only runtime metadata and timings, never molecular payloads.

## Optimization priorities

If batch latency becomes a hard acceptance gate, the highest-value work is:

1. explicit WASM SIMD profile-major DP across several profiles or sequences;
2. a gapped linear-memory filter accurate enough to reduce full Viterbi/Forward candidates further;
3. batched profile-major DP to reuse decoded scores and improve cache locality;
4. optional browser worker-pool parallelism for large batches.

These changes need expanded marginal/truncated golden data because overly aggressive filtering can preserve chain calls while losing the correct species, boundary, or alternative profile.
