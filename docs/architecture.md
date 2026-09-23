# Architecture and provenance

## Search and numbering

Inputs are normalized and checked before model allocation. Each eligible 128-state Plan7 profile then passes through the following pipeline:

1. A byte-quantized MSV recurrence rejects profiles above HMMER's default `F1 = 0.02` P-value threshold.
2. Large batches rank candidates with an eight-lane, rolling-row Viterbi pass. Native batches distribute profile work and exact finalist scoring across an explicit worker count.
3. Generic local/multihit Viterbi produces traceable M/I/D and N/B/E/J/C paths.
4. Overlapping paths are clustered into physical domains. Multihit Forward/Backward posterior decoding and fixed-seed stochastic traces define ambiguous envelopes. Compatible batch candidates share an eight-lane Forward/Backward pass.
5. Each envelope is rescored in local/unihit mode. Null1, null2 composition bias, and the outside-envelope length correction produce the domain score. Equal-length envelopes use eight-lane Forward/Backward and posterior probability decoding.
6. Posterior optimal-accuracy alignment is converted to IMGT positions. Framework deletions become gaps and CDR residues are distributed around the IMGT center positions using ANARCI's insertion ordering.
7. A long J-less suffix after position 104 receives ANARCI's permissive J recovery pass.
8. Optional germline assignment compares the 128-state projection with aligned V records, followed by J records from the assigned species.

The engine implements only the HMMER behavior needed by this pipeline. It does not include HMMER database/file APIs, the composition-bias, F2, or F3 filters, whole-sequence E-values, conditional-domain E-values, or native HMMER hit-table diagnostics. Single-sequence matrices and all tracebacks use the scalar reference path; the MSV gate, batch candidate ranking, and compatible Forward/Backward/posterior work use SIMD. Sparse tails, stochastic traces, and null2 reduction also remain scalar.

The public E-value is HMMER's independent-domain value:

```text
29 × exp(-lambda × (bitScore - tau))
```

The search-space size remains 29 when callers filter chains or species because ANARCI searches the complete profile database before applying those filters.

Match state `k` initially maps to IMGT position `k`. CDR1, CDR2, and CDR3 use positions 27–38, 56–65, and 105–117. Short CDRs alternate residues from the left and right edges; long CDRs add insertion codes around the central anchors. Insertion codes continue `A…Z, AA, BB…`.

## Embedded data

`assets/profiles.bin` uses the `ANRCPRF3` format. It stores profile metadata, MSV and Forward calibration, consensus residues, local-entry scores, match emissions, and transitions. Finite values use signed 24-bit fixed point at 1/32768 natural-log units; `-2^23` represents an impossible score. Insert emissions are omitted because configured protein insert log odds are zero.

`assets/germlines.bin` uses `ANRCGER1`. Each ordered V/J record contains its chain, species, gene name, and a 128-symbol alignment packed into 80 bytes with a five-bit alphabet.

Both decoders validate their magic, dimensions, strings, indices, symbols, calibration, arithmetic, truncation, and trailing bytes. Scores are expanded once and cached; model metadata continues to borrow the embedded bytes.

The model source is pinned in `assets/MANIFEST.toml`:

- official [oxpig/ANARCI](https://github.com/oxpig/ANARCI), commit `79f6c575056dedef86cb8f405ebb039197923eec`
- HMMER 3.4
- IMGT/GENE-DB 3.1.43, snapshot 2026-09-23, release `202638-7`
- 29 profiles and 2,450 aligned V/J records

The ANARCI pipeline was rerun against the named IMGT snapshot. Its generated `ALL.hmm` was converted to `ANRCPRF3`; its ordered `all_germlines` table was converted to `ANRCGER1`. The manifest records source and output SHA-256 hashes, the official IMGT amino-acid export hash, every profile name, and record counts.

The weekly database workflow checks the current IMGT release, program version, and reference-export checksum. A change fails the workflow until a maintainer reviews and deliberately refreshes the assets.

## Determinism, concurrency, and safety

Profile-name order resolves score ties and source order resolves germline ties. Stochastic domain definition uses HMMER's seed 42. Scheduling cannot affect result order or scoring order within a sequence.

Native parallel batches normalize once, pack same-profile and equal-length SIMD groups across the request, assign profile work by estimated DP cells, and dynamically dispatch bounded finalist batches. Incomplete groups use the scalar reference path. Browser concurrency uses independent WASM instances in module workers. Both APIs require an explicit worker count.

The JavaScript runtime contains no analytics, sequence logging, filesystem access, or analysis network path. Models are part of the WASM module. Length limits and checked arithmetic bound all user-controlled allocation, and public failures are structured errors. The safe Rust core forbids unsafe code; the WASM bridge uses unsafe code only at its internal allocation and request-slice boundary.

Property tests and libFuzzer cover model decoders, sequences, and FASTA. Run the fuzz targets with:

```sh
cd fuzz
cargo +nightly fuzz run profile_decode
cargo +nightly fuzz run germline_decode
cargo +nightly fuzz run sequence_inputs
```

## Refreshing models and golden data

Model updates are deliberate compatibility changes. Install HMMER 3.4 and the locked Python build dependency, then run the pinned official ANARCI pipeline in a disposable checkout. `pyproject.toml` and `uv.lock` identify the official Git source and exact revision. Upstream generates its model files with a legacy installation command that does not produce an installable wheel with current build frontends, so these steps install only its locked Python dependency and use the checkout itself as the reference:

```sh
uv sync --locked --no-install-project --no-default-groups --group reference \
  --no-install-package anarci
source .venv/bin/activate

git clone https://github.com/oxpig/ANARCI /tmp/anarci-assets
git -C /tmp/anarci-assets checkout 79f6c575056dedef86cb8f405ebb039197923eec
(
  cd /tmp/anarci-assets/build_pipeline
  PATH="/tmp/anarci-assets/bin:$PATH" bash RUN_pipeline.sh
)

cp /tmp/anarci-assets/build_pipeline/curated_alignments/germlines.py \
  /tmp/anarci-assets/lib/python/anarci/germlines.py
mkdir -p /tmp/anarci-assets/lib/python/anarci/dat/HMMs
cp /tmp/anarci-assets/build_pipeline/HMMs/ALL.hmm* \
  /tmp/anarci-assets/lib/python/anarci/dat/HMMs/

cargo run --locked -p anarcism-modelgen --bin anarcism-modelgen -- \
  /tmp/anarci-assets/lib/python/anarci/dat/HMMs/ALL.hmm assets/profiles.bin

PYTHONPATH=/tmp/anarci-assets/lib/python \
  python tools/export_germlines.py /tmp/anarcism-germlines.tsv
cargo run --locked -p anarcism-modelgen --bin germlinegen -- \
  /tmp/anarcism-germlines.tsv assets/germlines.bin
```

`RUN_pipeline.sh` replaces generated directories in its checkout. Record the retrieval time, IMGT release and program version, reference and generated-file hashes, and inventory counts in `assets/MANIFEST.toml`. Update the Git revision in `pyproject.toml` and refresh `uv.lock` when the upstream revision changes, then regenerate the reference:

```sh
PYTHONPATH=/tmp/anarci-assets/lib/python .venv/bin/python tools/generate_golden.py

cargo test --locked -p anarcism-core --test golden
tools/build-js.sh
node tools/parity-report.mjs
```

This reads `tests/golden/corpus.json` and writes `tests/golden/corpus_reference.jsonl`: a header record with provenance and field layout, then one record per case holding native ANARCI's domains, scores, top alternative hits, germlines, range-encoded IMGT numbering, and padded alignment.

Review inventory changes, parity, package size, and upstream terms before committing generated files. License text and attribution are retained in `THIRD_PARTY_NOTICES.md`.
