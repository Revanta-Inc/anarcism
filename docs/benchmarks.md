# Size and performance results

Measurements below are from the current release artifact built with rustc 1.94.0 and Binaryen 132. Compression is deterministic Node zlib gzip level 9 and Brotli quality 11. The hard threshold is strictly less than 500,000 bytes for both totals.

## Complete distribution

The notice-inclusive final report is produced by `node tools/size-report.mjs`:

| File | Raw | gzip | Brotli |
|---|---:|---:|---:|
| `anarcism.wasm` | 611,167 | 239,494 | 202,560 |
| `index.js` | 4,043 | 1,333 | 1,153 |
| `index.d.ts` | 2,777 | 824 | 735 |
| `THIRD_PARTY_NOTICES.md` | 5,263 | 1,990 | 1,482 |
| `package.json` | 891 | 425 | 373 |
| **Complete distribution** | **624,141** | **244,066** | **206,303** |

This is 48.8% of the gzip limit and 41.3% of the Brotli limit. Run the report after any source change and treat its output as authoritative.

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
| `profiles.bin` | 211,988 | 103,243 | 90,376 |
| `germlines.bin` | 226,705 | 44,108 | 36,195 |

Before Binaryen, the Rust release WASM was 676,691 bytes raw / 230,591 gzip / 195,793 Brotli. `wasm-opt -Oz` reduced raw WASM to 611,167 bytes but increased gzip/Brotli to 239,494/202,560 bytes. The optimized artifact is retained for the release because `-Oz` materially reduces raw transfer/storage and is the specified release pipeline; both compressed totals remain far below budget.

## Browser performance

Measured in Playwright Headless Chromium 151.0.7922.34 on an Apple-Silicon arm64 laptop running macOS 27.0, with Node 24.13.1 driving the fixture. The benchmark reads local cached files and reports no network decompression time. Cold initialization includes local fetch, array-buffer creation, compilation, and instantiation; cached initialization receives an already compiled `WebAssembly.Module`.

| Measure | Actual | Target | Status |
|---|---:|---:|---|
| Cold initialization | 16.0 ms | <100 ms cached | Pass even cold |
| Cached initialization | 0.1 ms | <100 ms | Pass |
| Single VH median / p95 | 40.3 / 40.6 ms | <50 ms | Pass |
| Single VL median / p95 | 35.9 / 36.1 ms | <50 ms | Pass |
| 100 VH/VL pairs (200 sequences) | 7,648.2 ms | <1,000 ms | Miss |

The batch target is not met. Analysis is deterministic and single-threaded; batching currently removes JSON call overhead but does not vectorize DP across sequences. Meeting one second would require a materially different SIMD/batched filter implementation. The baseline intentionally does not require WASM SIMD or threads, and the shortfall is not hidden with result caching or excluded work.

The Node/V8 cross-check produced 17.29 ms cold initialization, 0.118 ms cached initialization, 41.54 ms median VH, 37.04 ms median VL, and 7,916.0 ms for 100 pairs.

## Running measurements

```sh
tools/build-browser.sh

cd browser
npm run size
npm run benchmark
npm run benchmark:browser
```

The browser benchmark starts a loopback-only static server, loads the same release files as the demo, warms each path three times, measures 25 single-domain iterations, and includes one complete 200-sequence batch call. It prints only runtime metadata and timings, never molecular payloads.

## Optimization priorities

If batch latency becomes a hard acceptance gate, the highest-value work is:

1. a fixed-width integer/SIMD ungapped filter across several profiles or sequences;
2. a gapped linear-memory filter accurate enough to reduce full Viterbi/Forward candidates further;
3. batched profile-major DP to reuse decoded scores and improve cache locality;
4. optional WASM SIMD feature detection with the current scalar build retained as fallback.

These changes need expanded marginal/truncated golden data because overly aggressive filtering can preserve chain calls while losing the correct species, boundary, or alternative profile.
