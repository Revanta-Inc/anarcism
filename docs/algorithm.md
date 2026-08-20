# Algorithm and compact formats

## Search pipeline

For each normalized input, the engine performs these bounded, deterministic stages:

1. Decode the embedded profile header and borrow score slices directly from the WASM data segment.
2. Apply chain/species filters.
3. Compute an ungapped local emission score for every eligible profile. Based on sequence length, retain enough profiles for the requested alternatives plus ordering slack. This is a size-focused filter, not an attempt to reproduce HMMER's byte/striped MSV implementation.
4. Run generic Plan7 local, multihit Viterbi on retained profiles. The DP includes M/I/D states and N/B/E/J/C special states and stores the full bounded matrix for deterministic traceback.
5. Cluster paths that overlap at least half of the shorter domain. Keep one path per profile in each physical-domain cluster.
6. Compute the HMMER-style local/multihit Forward score over each candidate domain envelope, subtract the length-dependent null1 score, and apply a deterministic trace-based null2 composition correction.
7. Reject candidates below `minBitScore`, order remaining hits by score then profile name, and return physical domains in input order.
8. Convert the winning M/I/D path to IMGT positions and regions. Framework deletions become `-` in the padded alignment. CDR residues are redistributed symmetrically around the IMGT center positions, matching ANARCI insertion ordering.
9. If requested, project match states into a 128-character string and choose the highest-identity V germline across allowed species, then the highest-identity J germline from the V-assigned species. Ties keep source order, as pinned Python `max()` does.

The baseline uses no WASM threads, filesystem, SIMD requirement, server, native executable, or runtime model fetch.

## HMMER-compatible subset

Implemented:

- protein Plan7 match emission log odds
- M→M/I/D, I→M/I, and D→M/D transitions
- configured local entry probabilities
- multihit local N/B/E/J/C behavior
- Viterbi traceback with multiple domains
- Forward/null1 bit scoring
- stored Forward `tau`/`lambda` calibration values

Not implemented:

- HMMER file/database APIs at runtime
- exact HMMER MSV or Viterbi filter thresholds
- optimized SIMD stripes
- posterior Backward decoding and posterior null2
- accurate final E-values
- complete HMMER diagnostics or domain-envelope reporting

E-values are intentionally omitted. Applying `exp(-lambda * (score - tau))` to the approximate final score would produce a precise-looking but non-HMMER-equivalent number.

## IMGT mapping

Match state `k` maps initially to IMGT position `k`. Insert and delete states are retained through traceback. Before returning residues:

- CDR1 uses positions 27–38 with nominal length 12.
- CDR2 uses positions 56–65 with nominal length 10.
- CDR3 uses positions 105–117 with nominal length 13.
- Short CDRs place residues from alternating left/right ends.
- Long CDRs add insertion codes around the central anchors using ANARCI's symmetric order; insertion codes continue `A…Z, AA…` without an unbounded table.
- A J-less isolated truncation after the conserved position-104 cysteine uses ANARCI-compatible short-tail handling.

`NumberedResidue.sequenceIndex` always refers to the normalized input. `start` is the minimum numbered index and `end` is one past the maximum.

## `ANRCPRF1`

The profile generator parses HMMER3/f text at build time and validates a 20-residue alphabet and 128 match states. The little-endian asset contains:

- magic/version, score scale, profile count, model length
- profile and species strings, chain/receptor byte, source checksum
- Forward `tau` and `lambda`
- 128-byte consensus
- 128 quantized local-entry scores
- 128 × 20 quantized match log-odds scores
- 127 × 7 quantized transition log-probabilities

Impossible scores reserve `i16::MIN`; finite values use 1/1024-natural-log quantization. Parsing validates dimensions, UTF-8, calibration values, chain/receptor consistency, truncation, invalid codes, and trailing bytes.

## `ANRCGER1`

The germline asset contains an ordered species table followed by V/J records. Each record stores segment, chain, species index, gene name, and exactly 80 bytes for a 128-symbol five-bit alignment over `-ACDEFGHIKLMNPQRSTVWY`. Parsing validates all dimensions, indices, symbols, truncation, and trailing bytes.

The model database borrows embedded byte slices; it does not inflate a second floating-point copy of all profiles. Individual scores are decoded to `f32` on access.

## Determinism and complexity

Tie breaks use profile name after score and preserve germline source order. There is no randomized algorithm or hash-map iteration in public results.

For selected profiles, time is `O(P × L × 128)` and the largest Viterbi allocation for one profile is `O(L × 128 × 3)`, with `L ≤ 10,000`. Forward uses two DP rows. All user-controlled collection sizes are checked before allocation.
