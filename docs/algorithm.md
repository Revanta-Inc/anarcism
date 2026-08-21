# Algorithm and compact formats

## Search pipeline

For each normalized input, the engine performs these bounded, deterministic stages:

1. Decode the embedded profile header and borrow score slices directly from the WASM data segment.
2. Apply chain/species filters.
3. Compute an ungapped local emission score for every eligible profile. Based on sequence length, retain enough profiles for the requested alternatives plus reserve candidates. This is a size-focused filter, not an attempt to reproduce HMMER's byte/striped MSV implementation.
4. Run generic Plan7 local, multihit Viterbi on retained profiles. The DP includes M/I/D states and N/B/E/J/C special states and stores the full bounded matrix for deterministic traceback.
5. Cluster paths that overlap at least half of the shorter domain. Keep one path per profile in each physical-domain cluster, then use Viterbi trace scores to remove reserve profiles. The winner and every requested alternative continue to full scoring.
6. Run multihit Forward/Backward decoding for each finalist to locate posterior regions. A simple region becomes one envelope. When the posterior supports multiple domains, draw 200 deterministic seeded stochastic traces, single-link their domain segments, choose supported endpoints, and compute the trace-ensemble null2 correction.
7. Rescore each isolated envelope in local/unihit mode with full Forward and Backward matrices. Decode posterior state occupancy, compute expectation-based null2 composition bias where applicable, subtract null1, and apply HMMER's full-target outside-envelope length correction.
8. Reject candidates below `minBitScore`, rank them by calibrated significance, and return physical domains in input order. A 0.01 multidomain significance uncertainty band uses raw domain score and then profile name as deterministic tie breaks.
9. Build a posterior optimal-accuracy M/I/D display for the winning profile, then convert it to IMGT positions and regions. Framework deletions become `-` in the padded alignment. CDR residues are redistributed symmetrically around the IMGT center positions, matching ANARCI insertion ordering.
10. For a single domain that ends before the J region while leaving a long suffix, perform ANARCI's permissive second profile scan after IMGT 104 and splice the recovered J trace around the reconstructed CDR3.
11. If requested, project match states into a 128-character string and choose the highest-identity V germline across allowed species, then the highest-identity J germline from the V-assigned species. Ties keep source order, as pinned Python `max()` does.

The release build requires WASM SIMD128, but uses no WASM threads, shared memory, filesystem, server, native executable, or runtime model fetch.

## HMMER-compatible subset

Implemented:

- protein Plan7 match emission log odds
- M→M/I/D, I→M/I, and D→M/D transitions
- configured local entry probabilities
- multihit local N/B/E/J/C behavior
- Viterbi traceback with multiple domains
- multihit Forward/Backward posterior region detection
- fixed-seed stochastic traceback, segment clustering, and supported endpoint selection
- unihit Forward/Backward posterior decoding
- trace-ensemble and expectation-based null2 plus null1 bit scoring
- HMMER's 0.001-nat table-driven floating-point logsum approximation
- posterior optimal-accuracy alignment and traceback
- stored Forward `tau`/`lambda` calibration values

Not implemented:

- HMMER file/database APIs at runtime
- exact HMMER MSV or Viterbi filter thresholds
- HMMER's hand-striped SIMD matrix kernels (LLVM may vectorize suitable operations because the release enables `simd128`; the seeded traceback preserves optimized striped choice order)
- accurate final E-values
- complete HMMER diagnostics or domain-envelope reporting

E-values are intentionally omitted. Scalar log-space arithmetic and fixed-point profile scores are not bit-identical to HMMER's optimized probability-space pipeline, and calibration parameters alone do not reproduce its complete final E-value accounting.

## IMGT mapping

Match state `k` maps initially to IMGT position `k`. Insert and delete states are retained through traceback. Before returning residues:

- CDR1 uses positions 27–38 with nominal length 12.
- CDR2 uses positions 56–65 with nominal length 10.
- CDR3 uses positions 105–117 with nominal length 13.
- Short CDRs place residues from alternating left/right ends.
- Long CDRs add insertion codes around the central anchors using ANARCI's symmetric order; insertion codes continue `A…Z, AA, BB…`.
- A J-less isolated truncation after the conserved position-104 cysteine uses ANARCI's low-threshold J recovery and short-tail behavior.

`NumberedResidue.sequenceIndex` always refers to the normalized input. `start` is the minimum numbered index and `end` is one past the maximum.

## `ANRCPRF2`

The profile generator parses HMMER3/f text at build time and validates a 20-residue alphabet and 128 match states. The little-endian asset contains:

- magic/version, score scale, profile count, model length
- profile and species strings, chain/receptor byte, source checksum
- Forward `tau` and `lambda`
- 128-byte consensus
- 128 quantized local-entry scores
- 128 × 20 quantized match log-odds scores
- 127 × 7 quantized transition log-probabilities

Each score occupies three little-endian bytes. Impossible scores reserve signed 24-bit value `-2^23`; finite values use 1/32768-natural-log quantization. Parsing validates dimensions, UTF-8, calibration values, chain/receptor consistency, truncation, invalid codes, and trailing bytes.

## `ANRCGER1`

The germline asset contains an ordered species table followed by V/J records. Each record stores segment, chain, species index, gene name, and exactly 80 bytes for a 128-symbol five-bit alignment over `-ACDEFGHIKLMNPQRSTVWY`. Parsing validates all dimensions, indices, symbols, truncation, and trailing bytes.

The model database borrows embedded byte slices; it does not inflate a second floating-point copy of all profiles. Individual scores are decoded to `f32` on access.

## Determinism and complexity

Tie breaks use profile-name order after score and preserve germline source order. Stochastic domain definition uses HMMER's fixed seed of 42; no system entropy or hash-map iteration affects public results.

For selected profiles, time is `O(P × L × 128)`. Viterbi, Forward/Backward posterior decoding, and optimal-accuracy alignment use bounded `O(L × 128 × 3)` matrices, with at most two full matrices live in a stage and `L ≤ 10,000`. All user-controlled collection sizes are checked before allocation.
