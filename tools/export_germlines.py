#!/usr/bin/env python3
"""Export the pinned ANARCI germline table for the Rust asset generator.

This is a build-time provenance tool. The generated browser package never
imports Python or ANARCI.
"""

from __future__ import annotations

import argparse
from importlib.metadata import version
from pathlib import Path

from anarci.germlines import all_germlines


PINNED_ANARCI = "2026.2.13.2"
VALID_CHAINS = frozenset("HKLABGD")
VALID_RESIDUES = frozenset("-ACDEFGHIKLMNPQRSTVWY")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    observed_version = version("anarci")
    if observed_version != PINNED_ANARCI:
        raise SystemExit(
            f"expected ANARCI {PINNED_ANARCI}, found {observed_version}"
        )

    rows: list[tuple[str, str, str, str, str]] = []
    # Preserve the source dictionaries' order. Python's max() keeps the first
    # equal-identity germline, so this order is part of compatibility.
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
                        raise ValueError("germline metadata cannot contain tabs or newlines")
                    rows.append((segment, chain, species, gene, sequence))

    with args.output.open("w", encoding="utf-8", newline="\n") as output:
        output.write(f"# ANARCI {PINNED_ANARCI} germlines.py\n")
        output.write("# segment\tchain\tspecies\tgene\taligned_sequence\n")
        for row in rows:
            output.write("\t".join(row))
            output.write("\n")

    print(f"exported {len(rows)} germlines to {args.output}")


if __name__ == "__main__":
    main()
