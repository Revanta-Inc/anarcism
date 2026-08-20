# Feasibility gate

The gate was completed before expanding the proof of concept into the browser API. Its conclusion was that the full requested profile inventory and germline table can legally and technically fit below 500,000 compressed bytes.

## Required behavior

The pinned reference is ANARCI 2026.2.13.2 at commit `edcc29a08c40ac5acd49ce09f60a5ebfb7ccdd0c`, executed with HMMER 3.4 and the IMGT scheme. The required subset is:

1. Normalize canonical protein input with bounded allocation.
2. Search the ANARCI H/K/L/A/B/G/D Plan7 profiles in local, multihit mode.
3. Recover gapped match/insert/delete state paths and non-overlapping domains.
4. Apply ANARCI's IMGT CDR redistribution and insertion-code rules.
5. Rank profiles with Forward/null1 bit scores and return requested alternatives.
6. Optionally compare the 128 match-state projection with aligned V/J germlines.

The proof of concept showed that a full HMMER pipeline was unnecessary:

- A cheap ungapped local filter is useful for profile ordering, but byte-accurate HMMER MSV output is not a public requirement.
- Generic Plan7 Viterbi with traceback is required for boundaries and numbering.
- Forward is required for close bit-score parity.
- HMMER's full acceleration pipeline, sequence database layer, file parsers, SIMD code, threading, and C/Easel runtime are not required.
- Calibration parameters alone are insufficient for accurate E-values because the final score also requires posterior null2 composition correction. The public optional field therefore remains absent.

## Inventory and initial budget

The pinned `ALL.hmm` contains 29 profiles of length 128. The exact inventory is recorded in `assets/MANIFEST.toml`. Original source sizes were:

| Source | Raw bytes | Purpose |
|---|---:|---|
| `ALL.hmm` | 1,755,428 | 29 H/K/L/A/B/G/D HMMER profiles |
| `germlines.py` | 351,633 | 2,389 aligned V/J records |
| `schemes.py` | 86,509 | Numbering behavior reference only; not embedded |

The byte-level gate allocated the 500,000-byte compressed ceiling as follows:

| Component | Gate budget (Brotli) | Actual standalone Brotli |
|---|---:|---:|
| Quantized profiles | 120,000 | 90,376 |
| Packed germlines | 50,000 | 36,195 |
| Rust engine + JSON ABI | 250,000 | Not independently additive; WASM including both assets is 202,560 |
| JavaScript, declarations, package metadata/notices | 30,000 | Under 6,000 |
| Contingency | 50,000 | More than 290,000 remains against the final package |

The compact profile file is 211,988 bytes raw. It removes insert emissions whose configured protein log-odds score is zero, shares the alphabet, precomputes local entries, and quantizes all retained scores to signed 16-bit values at 1/1024 natural-log units. The germline file is 226,705 bytes raw and encodes 128 aligned symbols with five bits each.

## Proof-of-concept result

The initial VH and VL sequences produced the same chain, species, half-open numbered boundary, every IMGT label, and padded alignment as ANARCI. Scores were:

| Sequence | Pinned ANARCI | Rust proof of concept | Difference |
|---|---:|---:|---:|
| Mouse VH | 176.2 | 176.378 | +0.178 |
| Mouse K | 173.3 | 173.049 | -0.251 |

The initial compact profile asset was 103,243 bytes gzip / 90,376 bytes Brotli. That evidence supported continuing to the complete implementation.

## Legal gate

ANARCI and HMMER use permissive BSD three-clause terms, subject to preservation of their notices. IMGT currently provides its data and metadata under CC BY 4.0. The release includes attribution, change notices, and the relevant BSD license text in `THIRD_PARTY_NOTICES.md`. `assets/MANIFEST.toml` traces each transformed asset to the exact pinned input and SHA-256 digest.

This is an engineering audit, not legal advice. It establishes a documented redistribution basis and preserves the stated conditions; downstream distributors should perform their own review for their use and jurisdiction.

## Gate decision

Proceed. The implemented, notice-inclusive distribution is roughly 41% of the 500,000-byte Brotli ceiling and 49% of the gzip ceiling. No reduced-scope tier is needed.
