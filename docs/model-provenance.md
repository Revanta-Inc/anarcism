# Model provenance and licensing audit

## Pinned source

PyPI's provenance attestation for ANARCI 2026.2.13.2 identifies tag `v2026.02.13.2` and commit `edcc29a08c40ac5acd49ce09f60a5ebfb7ccdd0c`. The locked wheel SHA-256 is `8f88d49857abc2c5427c6cc903ebdac7b35e6636635eaa9d60bd8f0d4e7d14b7`.

Every embedded input/output checksum and every model name is recorded in `assets/MANIFEST.toml`. The important chain/species inventory is:

| Species | Profiles |
|---|---|
| human | H, K, L, A, B, G, D |
| mouse | H, K, L, A, B, G, D |
| rat | K, L |
| rabbit | H, K, L |
| rhesus | H, K, L |
| pig | H, K, L |
| alpaca | H |
| cow | H, K, L |

All profiles have 128 HMM match states. ANARCI terminal numbering can end at 127 for some chains because of germline/J coverage; that behavior is separate from HMM profile length.

## Transformation chain

### Profiles

1. Read `anarci/dat/HMMs/ALL.hmm` from the pinned wheel.
2. Validate HMMER3/f, alphabet, node count, transitions, checksums, and Forward calibration.
3. Convert match negative-log probabilities to log odds using HMMER's protein background.
4. Precompute HMMER local-entry values from occupancy.
5. Omit protein insert emissions, whose configured log-odds values are zero.
6. Quantize retained scores at 1/32768 natural-log units in signed 24-bit fields and emit `ANRCPRF2`.

### Germlines

1. `tools/export_germlines.py` imports only the exact pinned `all_germlines` dictionary.
2. Validate segment, chain, 128-column length, symbols, and metadata.
3. Preserve source dictionary order for equal-identity tie behavior.
4. Encode each aligned sequence in 80 bytes and emit `ANRCGER1`.

The two output hashes are:

- `profiles.bin`: `385ffd02779d776b414250d72c53cff6b96c6c5001a0cfad3d0fda0c1f9dbc83`
- `germlines.bin`: `2498aac3bc732ec5abab31768c38854e9d6a1a8d127aa0821454fe19ccfd9d7e`

## License findings

### ANARCI

The pinned wheel includes the original BSD three-clause `LICENCE`, copyright 2019 Charlotte Deane, James Dunbar, Alexsandr Kovaltsuk, and Claire Marks. It permits modified source and binary redistribution when the copyright, conditions, and disclaimer are retained and names are not used for endorsement.

The transformed HMM and germline files are treated conservatively as binary redistributions/derivatives of ANARCI material. The full notice and a description of changes ship in `THIRD_PARTY_NOTICES.md` and are copied into the counted browser package.

### HMMER

HMMER 3.4 is BSD three-clause. This project does not compile its C or Easel code and does not copy its SIMD implementation. HMMER 3.4 was the behavioral source for independently expressed Plan7 recurrences and the build-time profile interpretation, so its notice is preserved as well.

### IMGT

ANARCI documents that its HMMs and germline tables are built from IMGT germline data, and the data use IMGT standardized numbering. As of 1 July 2026, IMGT states that its data and metadata are available to public and private users under CC BY 4.0; requests are required for use of IMGT tools. This project does not invoke or redistribute an IMGT tool. It redistributes transformed data obtained through ANARCI, provides attribution, links the CC BY 4.0 terms, identifies modifications, and avoids endorsement language.

Relevant upstream pages:

- <https://pypi.org/project/ANARCI/>
- <https://github.com/jnooree/ANARCI/tree/edcc29a08c40ac5acd49ce09f60a5ebfb7ccdd0c>
- <https://github.com/EddyRivasLab/hmmer/tree/hmmer-3.4>
- <https://www.imgt.org/about/termsofuse.php>
- <https://www.imgt.org/about/CitingIMGT.php>

## Redistribution conclusion

Under the stated upstream terms, modified/quantized profiles and packed germline records may be redistributed with the included conditions, attribution, and modification notices. This conclusion is limited to the audited artifacts and is not legal advice. It should be revisited when the ANARCI pin, source data, or IMGT terms change.

## Scientific citations

- Dunbar J, Deane CM. “ANARCI: antigen receptor numbering and receptor classification.” *Bioinformatics* 32(2), 2016. DOI: `10.1093/bioinformatics/btv552`.
- Lefranc M-P et al. “IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains.” *Developmental & Comparative Immunology* 27(1), 2003. DOI: `10.1016/S0145-305X(02)00039-3`.
- Eddy SR. “Accelerated Profile HMM Searches.” *PLoS Computational Biology* 7(10), 2011. DOI: `10.1371/journal.pcbi.1002195`.
