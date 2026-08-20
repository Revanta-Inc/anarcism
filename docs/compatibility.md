# Compatibility with pinned ANARCI

## Reference and corpus

The reference is ANARCI 2026.2.13.2 (`edcc29a08c40ac5acd49ce09f60a5ebfb7ccdd0c`) with HMMER 3.4 and `scheme="imgt"`. `tools/generate_golden.py` records the native output in `tests/golden/corpus.json`; ordinary tests consume that file and do not invoke Python or HMMER.

The 15-case corpus contains:

- mouse VH and kappa VL feasibility sequences
- synthetic human lambda and cow VH domains built deterministically from pinned V/J tables
- human TCR alpha, beta, gamma, and delta domains
- N- and C-terminal truncations
- flanking residues around a complete domain
- a CDR3 insertion and CDR1 deletion
- a two-domain VL–VH scFv
- lysozyme as a negative control
- valid and swapped VH/VL pair cases

Run the machine-readable report with:

```sh
tools/build-browser.sh
node tools/parity-report.mjs
```

Current result:

| Measure | Result |
|---|---:|
| Cases | 15 |
| Exact domain calls/class/species/boundaries/padded alignments | 15 / 15 |
| Exact residue index/position/insertion-code tuples | 1,719 / 1,719 |
| Pair outcomes | 2 / 2 |
| Maximum absolute bit-score difference | 3.3592 bits |

This is full numbering parity for the checked-in corpus, not merely chain-classification parity.

## Intentional differences

### Coordinates

Public `start`/`end` are zero-based and half-open. ANARCI's numbered tuple exposes an inclusive end even though other internal fields are inconsistent for some terminal profiles. Golden generation converts the numbered tuple to `[start, end + 1)`.

### Strict species filtering

`allowedSpecies` is an exact case-insensitive filter over bundled profiles and germlines. An empty list is invalid. If the selected species has no matching profile, no domain is returned. ANARCI has paths that fall back or make unavailable germline species nonfatal; this project deliberately does not broaden a caller's filter.

### Scores and E-values

Forward/null1 is implemented, but null2 uses state occupancy from the deterministic Viterbi trace rather than HMMER's Forward/Backward posterior expectations. Strong single domains are usually within a fraction of a bit; the multidomain corpus currently reaches 3.3592 bits.

Because final scores are not exactly HMMER-equivalent, `eValue` is omitted rather than applying stored calibration to an approximate score. The profile format retains `tau` and `lambda` for a future posterior-null2 implementation.

### Profile acceleration

The initial ungapped filter is used only to select profiles for expensive gapped scoring. It is not HMMER's optimized MSV filter. The retained count scales with requested alternatives and estimated domain count. This has exact profile ranking on the corpus, but a marginal exotic profile outside the retained set could be absent from alternatives. Raising `alternativeHitCount` causes more profiles to receive full scoring; 28 scores the entire inventory.

### Input validation

Only the 20 canonical amino acids are accepted. ANARCI can pass ambiguous symbols through some native paths. ASCII whitespace and lowercase are normalized with explicit warnings; other symbols produce `INVALID_SEQUENCE`.

### Terminal handling

ANARCI's displayed alignments, numbered tuple, and HMMER `query_end` can disagree by one residue for unusual J records. This implementation defines its domain boundary from returned numbered residues. It reproduces corpus truncations and the pinned K/L/B terminal behavior rather than exposing conflicting fields.

## Germline parity

Assignment compares only germline non-gap match positions, matching `get_identity` in pinned ANARCI. On the feasibility sequences:

| Sequence | V assignment / identity | J assignment / identity |
|---|---|---|
| Mouse VH | `IGHV14-4*02`, 0.9489796 | `IGHJ3*01`, 0.9285714 |
| Mouse K | `IGKV6-13*01`, 0.9891304 | `IGKJ5*01`, 1.0 |

The V-assigned species constrains J assignment. If strict filters leave no eligible V record, `germline` is absent.

## Adding cases

Add a sequence and category to `tools/generate_golden.py`, regenerate with the exact locked Python environment and HMMER 3.4, inspect the JSON diff for molecular-data changes, then run both the Rust golden test and browser parity report. Do not regenerate reference data silently during ordinary CI.
