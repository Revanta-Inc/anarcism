# Compatibility with pinned ANARCI

## Reference and corpus

The reference is ANARCI 2026.2.13.2 (`edcc29a08c40ac5acd49ce09f60a5ebfb7ccdd0c`) with HMMER 3.4 and `scheme="imgt"`. `tools/generate_golden.py` records native output in the checked-in fixtures; ordinary tests consume those files and do not invoke Python or HMMER.

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

The expanded 1,397-case corpus adds published and germline VH/VL sequences, all bundled species, TCRs, systematic N/C truncations, CDR-length ladders, framework indels, edge thresholds, scFvs and four-domain constructs, long/ultralong CDR3s, constant domains, and negative proteins. Its full ANARCI output is stored as compact FNV-1a fingerprints in `corpus_v2_reference.json`; the unhashed chain, species, boundary, count, and score fields remain directly inspectable.

Run the machine-readable report with:

```sh
npm --prefix browser run build
npm --prefix browser run parity
```

The parity runner shards the expanded corpus across `os.availableParallelism()` worker threads. Set `PARITY_WORKERS=N` to override that count.

Current result:

| Measure | Result |
|---|---:|
| Cases | 1,412 |
| Exact domain calls/class/species/boundaries/padded alignments | 1,414 / 1,414 |
| Exact residue index/position/insertion-code tuples | 154,968 / 154,968 |
| Pair outcomes | 2 / 2 |
| Maximum absolute bit-score difference | 0.8449 bits |

This is full numbering parity for the checked-in corpus, not merely chain-classification parity.

## Intentional differences

### Coordinates

Public `start`/`end` are zero-based and half-open. ANARCI's numbered tuple exposes an inclusive end even though other internal fields are inconsistent for some terminal profiles. Golden generation converts the numbered tuple to `[start, end + 1)`.

### Strict species filtering

`allowedSpecies` is an exact case-insensitive filter over bundled profiles and germlines. An empty list is invalid. If the selected species has no matching profile, no domain is returned. ANARCI has paths that fall back or make unavailable germline species nonfatal; this project deliberately does not broaden a caller's filter.

### Scores and E-values

Domain definition uses multihit Forward/Backward posterior regions, HMMER-compatible fixed-seed stochastic trace clustering for ambiguous regions, and unihit posterior decoding with null2 correction for isolated envelopes. Remaining score drift comes from 1/32768-nat profile quantization and scalar log-space arithmetic differing slightly from HMMER's optimized probability-space implementation.

Because final scores are not bit-identical to HMMER, `eValue` is omitted rather than applying stored calibration to an approximate score. The profile format retains `tau` and `lambda` for deterministic significance ordering and future exact E-values.

### Profile acceleration

The initial ungapped filter is used only to select profiles for expensive gapped scoring. It is not HMMER's optimized MSV filter. The retained count scales with requested alternatives and estimated domain count, with at least one reserve finalist. Viterbi trace scores remove that reserve per physical-domain cluster before posterior rescoring; the winner and every returned alternative still receive full Forward/Backward scoring. This has exact profile ranking on the corpus, but a marginal exotic profile outside the retained set could be absent. Raising `alternativeHitCount` expands both stages; 28 scores the entire inventory.

### Input validation

Only the 20 canonical amino acids are accepted. ANARCI can pass ambiguous symbols through some native paths. ASCII whitespace and lowercase are normalized with explicit warnings; other symbols produce `INVALID_SEQUENCE`.

### Terminal handling

ANARCI's displayed alignments, numbered tuple, and HMMER `query_end` can disagree by one residue for unusual J records. This implementation defines its domain boundary from returned numbered residues. It reproduces ANARCI's terminal junction buffering, CDR padding, and permissive second J-region scan for long CDR3s rather than exposing conflicting fields.

## Germline parity

Assignment compares only germline non-gap match positions, matching `get_identity` in pinned ANARCI. On the feasibility sequences:

| Sequence | V assignment / identity | J assignment / identity |
|---|---|---|
| Mouse VH | `IGHV14-4*02`, 0.9489796 | `IGHJ3*01`, 0.9285714 |
| Mouse K | `IGKV6-13*01`, 0.9891304 | `IGKJ5*01`, 1.0 |

The V-assigned species constrains J assignment. If strict filters leave no eligible V record, `germline` is absent.

## Adding cases

Add built-in cases to `tools/generate_golden.py`, or add `{id, seq, category, note}` records to `corpus_v2.json` and regenerate its compact reference with:

```sh
.venv/bin/python tools/generate_golden.py \
  tests/golden/corpus_v2_reference.json \
  --cases-file tests/golden/corpus_v2.json --compact
```

Use the exact locked Python environment and HMMER 3.4, inspect molecular-data diffs, then run both the Rust golden test and browser parity report. Reference data is never regenerated during ordinary CI.
