# Size and performance results

Measurements below are from the current release artifact built with rustc 1.94.0 and Binaryen 132. Compression is deterministic Node zlib gzip level 9 and Brotli quality 11. The hard threshold is strictly less than 500,000 bytes for both totals.

## Complete distribution

The notice-inclusive final report is produced by `node tools/size-report.mjs`:

| File | Raw | gzip | Brotli |
|---|---:|---:|---:|
| `anarcism.wasm` | 757,038 | 305,404 | 258,404 |
| `index.js` | 4,420 | 1,434 | 1,237 |
| `index.d.ts` | 2,853 | 842 | 740 |
| `worker-pool.js` | 11,842 | 3,243 | 2,818 |
| `worker-pool.d.ts` | 1,368 | 516 | 435 |
| `worker.js` | 1,666 | 586 | 470 |
| `THIRD_PARTY_NOTICES.md` | 5,288 | 2,005 | 1,491 |
| `package.json` | 1,070 | 462 | 401 |
| **Complete distribution** | **785,545** | **314,492** | **265,996** |

This is 62.9% of the gzip limit and 53.2% of the Brotli limit. Run the report after any source change and treat its output as authoritative.

The report counts:

- optimized `anarcism.wasm`, including both embedded assets
- JavaScript glue
- asynchronous worker-pool runtime and worker entry point
- TypeScript declarations
- required third-party notice
- npm package metadata

It does not add standalone asset sizes to the package total because that would double-count bytes already inside WASM.

Standalone asset attribution from the same report:

| Embedded input | Raw | gzip | Brotli |
|---|---:|---:|---:|
| `profiles.bin` | 315,721 | 152,430 | 132,231 |
| `germlines.bin` | 226,705 | 44,108 | 36,195 |

Before Binaryen, the Rust release WASM was 830,883 bytes raw / 296,960 gzip / 251,705 Brotli. `wasm-opt -Oz` reduced raw WASM to 757,038 bytes but increased gzip/Brotli to 305,404/258,404 bytes. The optimized artifact is retained for the release because `-Oz` materially reduces raw transfer/storage and is the specified release pipeline; both compressed totals remain below budget.

## Browser performance

Measured in Playwright Headless Chromium 151.0.7922.34 on an Apple-Silicon arm64 laptop running macOS 27.0, with Node 24.13.1 driving the fixture. The benchmark reads local cached files and reports no network decompression time. Cold initialization includes local fetch, array-buffer creation, compilation, and instantiation; cached initialization receives an already compiled `WebAssembly.Module`.

| Measure | Actual | Target | Status |
|---|---:|---:|---|
| Cold initialization | 8.0 ms | <100 ms cached | Pass even cold |
| Cached initialization | 0.5 ms | <100 ms | Pass |
| Single VH median / p95 | 33.6 / 33.8 ms | <50 ms | Pass |
| Single VL median / p95 | 29.8 / 30.1 ms | <50 ms | Pass |
| 100 VH/VL pairs, synchronous (200 sequences) | 6,339.0 ms | <1,000 ms | Miss |
| 100 VH/VL pairs, 8-worker pool (200 sequences) | 798.4 ms | <1,000 ms | Pass |

Both single-sequence targets pass. The synchronous API remains above the batch-latency target, while the opt-in eight-worker pool brings the same workload below it without WASM threads or shared memory. Replacing per-cell transcendental logsum operations with HMMER's table-driven approximation cut the original single-worker measurements by roughly 68%, and expanding packed profile scores once on database initialization removed another repeated inner-loop operation. The release enables WASM SIMD128 and LLVM emits vector instructions for suitable operations, but the Plan7 recurrences are not hand-vectorized.

Across three complete Node/V8 cross-checks, the median run produced 21.15 ms cold initialization, 0.475 ms cached initialization, 33.09/33.37 ms median/p95 VH, 29.55/29.80 ms median/p95 VL, and 6,280.5 ms for 100 pairs. The first numbering call initializes the 64 KiB logsum table and expands about 405 KiB of cached profile scores; subsequent calls reuse both.

### Browser worker-pool scaling

The asynchronous browser entry point creates an independent WASM engine in each module worker, balances inputs by sequence length, and restores results to input order. After warming two sequences per worker, the 200-sequence workload scaled as follows in Chromium:

| Workers | Pool initialization | Batch total | ms/sequence | sequences/s | Speedup vs synchronous |
|---:|---:|---:|---:|---:|---:|
| 1 | 10.1 ms | 6,392.5 ms | 31.963 | 31.287 | 0.99x |
| 4 | 15.8 ms | 1,618.8 ms | 8.094 | 123.548 | 3.92x |
| 8 | 24.7 ms | 798.4 ms | 3.992 | 250.501 | 7.92x |

Pool initialization is outside the hot batch totals. The one-worker control shows that messaging and structured cloning add about 0.6% on this workload. Eagerly creating a large pool can nevertheless regress a cold single-sequence interaction because it pays for engines that cannot help that request. The demo therefore initializes one worker and lazily grows to `min(configured workers, FASTA records)` only when a batch arrives.

## Cross-runtime comparison

The same VH and VL inputs were measured on August 20, 2026 on an Apple M1 Max with 10 logical CPUs. Each runtime used default numbering options, an IMGT scheme, an 80-bit threshold, no germline assignment, and one worker/thread. After three VH/VL warmups, the harness collected 25 hot single-sequence samples and timed one call containing 100 repeated VH/VL pairs (200 sequences).

| Runtime | VH median / p95 | VL median / p95 | 200-sequence batch | ms/sequence | sequences/s | Batch time vs Rust |
|---|---:|---:|---:|---:|---:|---:|
| Native Rust (`-O3`) | 19.8 / 19.9 ms | 17.6 / 17.8 ms | 3,745.5 ms | 18.728 | 53.397 | 1.00x |
| Browser WASM | 33.7 / 33.9 ms | 30.0 / 30.2 ms | 6,384.9 ms | 31.925 | 31.324 | 1.70x |
| ANARCI/HMMER | 28.5 / 29.1 ms | 24.4 / 24.6 ms | 3,774.7 ms | 18.874 | 52.984 | 1.01x |

Native Rust uses the dedicated `benchmark` Cargo profile with `opt-level = 3`. Browser WASM uses the shipping size-optimized `-Oz` artifact with SIMD128 enabled and calls the package API directly in Playwright Headless Chromium 151, excluding Web Worker messaging, DOM rendering, and initialization. ANARCI 2026.2.13.2 runs under Python 3.14.5 with HMMER 3.4 and `ncpu=1`. Process startup, imports, and WASM initialization are outside the hot timings; each timed ANARCI API call still includes its inherent temporary-file I/O and `hmmscan` subprocess.

The batch result is a public-API throughput comparison, not an equivalent DP-kernel microbenchmark. ANARCI submits all 200 queries to one HMMER scan, while the serial Rust comparison and synchronous WASM `numberSequences` implementation loop over 200 independent sequence calls. In this run their batch totals were effectively tied, with serial native Rust 0.8% faster. The browser row deliberately remains the one-worker baseline; optional worker-pool scaling is reported separately above.

### Native batch scaling

Native callers can use the opt-in `number_sequences_parallel` API with an explicit nonzero worker count. It uses scoped standard-library threads, retains input order and exact serial results, adds no dependency, and is omitted from WASM builds. The existing `number_sequences` API remains serial so library callers are never surprised by automatic CPU saturation.

The same 200-sequence workload, measured as the median of three complete batches under the native `-O3` profile, scales as follows:

| Workers | Batch total | ms/sequence | sequences/s | Speedup |
|---:|---:|---:|---:|---:|
| 1 | 3,735.6 ms | 18.678 | 53.539 | 1.00x |
| 2 | 1,860.0 ms | 9.300 | 107.526 | 2.01x |
| 4 | 932.8 ms | 4.664 | 214.409 | 4.00x |
| 8 | 473.9 ms | 2.369 | 422.046 | 7.88x |
| 10 | 464.0 ms | 2.320 | 431.051 | 8.05x |

Four workers bring the native batch below the 1,000 ms target. The smaller gain from 8 to 10 workers is consistent with this machine's eight performance cores plus two efficiency cores; callers should benchmark an appropriate cap rather than assuming every logical CPU improves latency.

### Optimization attribution

An immediate same-harness Node control separated the effect of enabling `simd128` from the algorithmic change:

| Build | VH median | VL median | 100 pairs |
|---|---:|---:|---:|
| Previous scalar release | 183.46 ms | 161.93 ms | 35,198.3 ms |
| SIMD128 only | 184.15 ms | 162.48 ms | 34,744.1 ms |
| SIMD128 + HMMER logsum table + inverse score scale | 56.58 ms | 50.85 ms | 10,803.3 ms |
| + cached `f32` profile scores + fill-local logsum table reference | 41.79 ms | 37.34 ms | 7,914.1 ms |
| + one-sided `-∞` sentinel + fixed-stride DP/profile views | 33.09 ms | 29.55 ms | 6,280.5 ms |

SIMD alone was effectively neutral for single sequences; LLVM emitted vector instructions, but the dependent Plan7 recurrences remain scalar. The approximately 68% initial reduction comes primarily from replacing per-cell `exp`/`ln1p` operations with HMMER's table lookup. Expanding signed-24-bit profile scores once and acquiring the shared logsum table once per Forward/Backward fill then reduced the native 200-sequence batch from 4,971.4 to 3,731.5 ms (24.9%), while the packed asset itself remained unchanged. The SIMD-only optimized WASM was 748,016 bytes versus 752,616 bytes for the scalar control.

The newest pass exposes fixed-stride match/transition rows and treats each DP state cell as a three-float view inside Forward/Backward fills. Logsum explicitly handles only a left-hand `-∞`; a right-hand sentinel naturally takes the existing difference-cutoff path and returns the finite left operand. These changes are neutral in native Rust (3,740.5 ms before versus 3,745.5 ms in the final cross-runtime run), but reduce the median Node/V8 batch from 7,914.1 to 6,280.5 ms (20.6%) and the Chromium batch from 9,230.3 to 6,390.1 ms (30.8%).

Apple Time Profiler corroborates the attribution. Before cached score expansion and fill-local lookup reuse, packed-score decoding accounted for 36.4% of exclusive samples and logsum lookup plus `OnceLock` access for 15.2%. In the follow-up trace, decoding appeared in one first-use sample (0.02%) and no per-cell `OnceLock` access remained. Its next leading costs were scalar logsum's sentinel/absolute-difference operations and direct transition-score/DP-matrix accesses, which motivated the newest pass.

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

For further single-worker throughput gains, the highest-value work is:

1. explicit WASM SIMD profile-major DP across several profiles or sequences;
2. a gapped linear-memory filter accurate enough to reduce full Viterbi/Forward candidates further;
3. batched profile-major DP to improve cache locality.

These changes need expanded marginal/truncated golden data because overly aggressive filtering can preserve chain calls while losing the correct species, boundary, or alternative profile.
