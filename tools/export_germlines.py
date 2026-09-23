#!/usr/bin/env python3
"""Export the versioned ANARCI/IMGT germline table for the Rust asset generator.

This is a build-time provenance tool. The generated browser package never
imports Python or ANARCI.
"""

from __future__ import annotations

import argparse
import importlib
import tomllib
from hashlib import sha256
from pathlib import Path

from anarci.germlines import all_germlines

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = tomllib.loads((ROOT / "assets/MANIFEST.toml").read_text())
ANARCI_REPOSITORY = MANIFEST["reference"]["anarci_repository"]
ANARCI_COMMIT = MANIFEST["reference"]["anarci_commit"]
IMGT_GENEDB_PROGRAM_VERSION = MANIFEST["reference"][
    "imgt_genedb_program_version"
]
IMGT_SNAPSHOT = MANIFEST["reference"]["imgt_snapshot"]
EXPECTED_GERMLINES_SHA256 = MANIFEST["germlines"]["source_sha256"]
VALID_CHAINS = frozenset("HKLABGD")
VALID_RESIDUES = frozenset("-ACDEFGHIKLMNPQRSTVWY")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    germline_module = importlib.import_module("anarci.germlines")
    source_path = Path(germline_module.__file__ or "")
    source_hash = sha256(source_path.read_bytes()).hexdigest()
    if source_hash != EXPECTED_GERMLINES_SHA256:
        raise SystemExit(
            f"expected {IMGT_SNAPSHOT} germlines.py SHA-256 "
            f"{EXPECTED_GERMLINES_SHA256}, found {source_hash} at {source_path}"
        )

    rows: list[tuple[str, str, str, str, str]] = []
    # Source order resolves equal-identity germline ties.
    for segment in ("V", "J"):
        for chain, by_species in all_germlines[segment].items():
            if chain not in VALID_CHAINS:
                raise ValueError(f"unsupported chain {chain!r}")
            for species, genes in by_species.items():
                for gene, sequence in genes.items():
                    if len(sequence) != 128:
                        raise ValueError(
                            f"{segment}/{chain}/{species}/{gene} has length {len(sequence)}"
                        )
                    if not set(sequence) <= VALID_RESIDUES:
                        raise ValueError(
                            f"{segment}/{chain}/{species}/{gene} has invalid residues"
                        )
                    if any("\t" in field or "\n" in field for field in (species, gene)):
                        raise ValueError(
                            "germline metadata cannot contain tabs or newlines"
                        )
                    rows.append((segment, chain, species, gene, sequence))

    with args.output.open("w", encoding="utf-8", newline="\n") as output:
        output.write(
            f"# ANARCI {ANARCI_REPOSITORY} commit {ANARCI_COMMIT}, IMGT/GENE-DB "
            f"{IMGT_GENEDB_PROGRAM_VERSION} snapshot {IMGT_SNAPSHOT} germlines.py\n"
        )
        output.write("# segment\tchain\tspecies\tgene\taligned_sequence\n")
        for row in rows:
            output.write("\t".join(row))
            output.write("\n")

    print(f"exported {len(rows)} germlines to {args.output}")


if __name__ == "__main__":
    main()
