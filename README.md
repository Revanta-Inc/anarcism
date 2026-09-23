<p align="center">
  <img src="demo/public/anarcism.png" alt="anarcism logo" width="300">
</p>

# anarcism

`anarcism` is an independent reimplementation of [ANARCI](https://github.com/oxpig/ANARCI) by James Dunbar and Charlotte M. Deane. It recognizes and IMGT-numbers antibody and T-cell-receptor variable domains. The engine is written in Rust and ships as a self-contained Python extension and WebAssembly package: analysis runs locally, with no HMMER binary, model download, or backend request.

It supports H/K/L/A/B/G/D chains, multidomain inputs, FASTA, configurable profile filters, alternative hits, V/J germline assignment, and VH/VL pair validation.

## Features

Compared with the official ANARCI implementation:

- Self-contained installation without legacy hooks or HMMER
- Bounded-memory FASTA streaming with lower end-to-end time on large inputs
- Efficient multicore scaling
- WebAssembly package for browsers and Node.js

## Python and CLI

```sh
pip install anarcism
```

```python
import anarcism

result = anarcism.number_sequence(vh, assign_germline=True)
batch = anarcism.number_sequences(
    [("heavy", vh), ("light", vl)],
    workers=4,
)

# Iterable input and output keep large datasets bounded in memory.
with open("sequences.fasta") as fasta:
    for result in anarcism.iter_number_fasta(fasta, workers=4):
        consume(result)
```

The `anarcism` command and `python -m anarcism` provide ANARCI-compatible IMGT vertical, CSV, and hit-table output:

```sh
anarcism -i sequences.fasta -o numbered.anarci -p 8
anarcism -i sequences.fasta -o numbered --csv --assign_germline
anarcism -i sequences.fasta -o numbered.anarci -ht hits.txt
gunzip -c sequences.fasta.gz | anarcism -i - -o numbered.anarci
```

Run `anarcism --help` for the supported ANARCI flags. Only IMGT numbering is implemented; `--hmmerpath` is unnecessary because the backend is embedded.

`anarcism` is a drop-in replacement for the official ANARCI executable. To avoid changing scripts that call `ANARCI`, create a symlink:

```sh
ln -s "$(command -v anarcism)" "$(dirname "$(command -v anarcism)")/ANARCI"
```

## JavaScript

```sh
npm install @revanta/anarcism
```

The package runs in modern browsers and Node.js 20.16+. It requires WebAssembly SIMD128 but not WASM threads or shared memory.

```js
import { Anarcism } from "@revanta/anarcism";

const anarcism = await Anarcism.create();

const result = anarcism.numberSequence(vh, {
	allowedChains: ["H", "K", "L"],
	minBitScore: 80,
	alternativeHitCount: 3,
	assignGermline: true,
});

const batch = anarcism.numberSequences([
	{ id: "heavy", sequence: vh },
	{ id: "light", sequence: vl },
]);

const fastaResults = anarcism.numberFasta(`>heavy\n${vh}\n`);
const pair = anarcism.validateAntibodyPair(vh, vl);
```

`Anarcism.create()` loads the bundled `anarcism.wasm`; pass `{ source }` with a URL, `Response`, bytes, or a compiled `WebAssembly.Module` to load it from elsewhere. `Anarcism.createSync(bytesOrModule)` instantiates without awaiting. Each instance owns an independent WebAssembly engine.

The synchronous API is suitable for interactive calls. For large batches in the browser, the worker-pool entry point runs independent WASM engines in Web Workers without blocking the page:

```js
import { AnarcismWorkerPool } from "@revanta/anarcism/worker-pool";

const pool = await AnarcismWorkerPool.create({ workers: 4 });
const results = await pool.numberSequences(inputs);
pool.terminate();
```

Every asynchronous call accepts an `AbortSignal`. An aborted call rejects with `signal.reason`; the pool replaces the workers still computing it, so it stays usable at the same size:

```js
const controller = new AbortController();
const pending = pool.numberSequences(inputs, { signal: controller.signal });
controller.abort();

const anarcism = await Anarcism.create({ signal: AbortSignal.timeout(5_000) });
```

TypeScript declarations are included. Invalid input raises `AnarcismError` with a stable `code`, message, and optional input ID.

## Compatibility and size

The reference is [ANARCI](https://github.com/oxpig/ANARCI) (OPIG) with HMMER 3.4, pinned in [`assets/MANIFEST.toml`](assets/MANIFEST.toml). The parity corpus holds 1,397 sequences (1,398 domains). It combines natural antibody and TCR sequences, IMGT germlines, and published PDB chains with deliberately adversarial stress tests: synthetic CDR-length ladders, framework indels, truncations, scFvs, multidomain constructs, constant domains, and non-antibody decoys.

| criterion                        | result                          |
| -------------------------------- | ------------------------------- |
| domain detection                 | 1,397/1,397 sequences (100%)    |
| chain classification             | 1,398/1,398 domains (100%)      |
| species assignment               | 1,398/1,398 domains (100%)      |
| domain boundaries                | 1,398/1,398 domains (100%)      |
| IMGT numbering, exact per domain | 1,398/1,398 domains (100%)      |
| IMGT numbering, per residue      | 153,213/153,213 residues (100%) |
| germline V and J genes           | 1,398/1,398 domains (100%)      |

A domain counts as exactly numbered only if every residue's IMGT position and insertion code, and its padded IMGT alignment, match ANARCI. Bit scores stay within 0.15 bits of HMMER's (mean 0.03) and E-values within 7% (mean 1.2%). The reference values live in [`tests/golden/corpus_reference.jsonl`](tests/golden/corpus_reference.jsonl). The Rust and Python golden tests check every row against it, and `npm --workspace packages/anarcism run parity` reports domain, residue, and score parity for the WebAssembly build.

## Development

The build needs Rust and Node.js 20.16+. `rust-toolchain.toml` pins the Rust release and installs the `wasm32-unknown-unknown` target. Install [Binaryen](https://github.com/WebAssembly/binaryen) for `wasm-opt`; without it the build ships the larger, unoptimized module and prints a warning. Python and ANARCI are only needed to regenerate the reference data; they are not runtime dependencies.

```sh
cargo test --workspace
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings

npm ci
npm --workspace packages/anarcism run build
npm --workspace packages/anarcism test
npm --workspace packages/anarcism run test:browser
npm --workspace packages/anarcism run parity
npm --workspace demo run dev
```

The package build compiles `anarcism-wasm` for `wasm32-unknown-unknown`, optimizes it with `wasm-opt -Oz`, copies the JavaScript entry points, declarations, and license files next to it in `packages/anarcism/dist/`, and prints the size report.

## Attribution and citation

anarcism reproduces the methods and data of the projects below. If you use it in published work, please cite them.

- **ANARCI** (Oxford Protein Informatics Group) defines the recognition, receptor classification, and numbering that anarcism reimplements. anarcism is validated against [oxpig/ANARCI](https://github.com/oxpig/ANARCI) at commit `79f6c575056dedef86cb8f405ebb039197923eec`, and its embedded profile HMMs and germline tables are generated with ANARCI's own build pipeline.

  Dunbar J, Deane CM. ANARCI: antigen receptor numbering and receptor classification. _Bioinformatics_ 32(2):298–300 (2016). [doi:10.1093/bioinformatics/btv552](https://doi.org/10.1093/bioinformatics/btv552)

- **HMMER** 3.4 (Sean R. Eddy and the HMMER developers) provides the profile-HMM search and scoring that ANARCI relies on. anarcism reimplements the parts ANARCI uses and is validated against HMMER 3.4.

  Eddy SR. Accelerated profile HMM searches. _PLoS Computational Biology_ 7(10):e1002195 (2011). [doi:10.1371/journal.pcbi.1002195](https://doi.org/10.1371/journal.pcbi.1002195)

- **IMGT®**, the international ImMunoGeneTics information system® (Marie-Paule Lefranc and colleagues), defines the IMGT numbering scheme and supplies the germline sequences, from IMGT/GENE-DB release `202638-7`.

  Lefranc M-P, Pommié C, Ruiz M, Giudicelli V, Foulquier E, Truong L, Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. _Developmental & Comparative Immunology_ 27(1):55–77 (2003). [doi:10.1016/S0145-305X(02)00039-3](<https://doi.org/10.1016/S0145-305X(02)00039-3>)

  Giudicelli V, Chaume D, Lefranc M-P. IMGT/GENE-DB: a comprehensive database for human and mouse immunoglobulin and T cell receptor genes. _Nucleic Acids Research_ 33:D256–D261 (2005). [doi:10.1093/nar/gki010](https://doi.org/10.1093/nar/gki010)

License texts and data provenance for all three are in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

## License

The project's own code is licensed under the Apache License 2.0; see [`LICENSE`](LICENSE). Third-party attributions and license texts are in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md). The embedded profile and germline data have additional provenance and licensing considerations. Public distribution remains blocked until the IMGT-derived data question is resolved in writing.
