#!/usr/bin/env python3
"""Generate the checked-in compatibility corpus with pinned ANARCI/HMMER.

Only this development tool invokes native ANARCI. Browser and Rust tests read
the resulting JSON and never need Python or HMMER.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from importlib.metadata import version
from pathlib import Path
from typing import Any

from anarci import anarci
from anarci.germlines import all_germlines


PINNED_ANARCI = "2026.2.13.2"
ANARCI_COMMIT = "edcc29a08c40ac5acd49ce09f60a5ebfb7ccdd0c"
VH = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA"
VL = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR"
SCFV = "DIQMTQSPSSLSASVGDRVTITCRTSGNIHNYLTWYQQKPGKAPQLLIYNAKTLADGVPSRFSGSGSGTQFTLTISSLQPEDFANYYCQHFWSLPFTFGQGTKVEIKRTGGGGSGGGGSGGGGSGGGGSEVQLVESGGGLVQPGGSLRLSCAASGFDFSRYDMSWVRQAPGKRLEWVAYISSGGGSTYFPDTVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCARQNKKLTWFDYWGQGTLVTVSSHHHHHH"
LYSOZYME = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3


def synthetic_domain(chain: str, species: str) -> str:
    v = next(iter(all_germlines["V"][chain][species].values()))
    j = next(iter(all_germlines["J"][chain][species].values()))
    aligned: list[str] = []
    for position, (v_residue, j_residue) in enumerate(zip(v, j), start=1):
        if v_residue != "-":
            aligned.append(v_residue)
        elif j_residue != "-":
            aligned.append(j_residue)
        elif 105 <= position <= 117:
            aligned.append("A")
        else:
            aligned.append("-")
    return "".join(aligned).replace("-", "")


def cases() -> list[dict[str, str]]:
    return [
        {"id": "mouse_vh", "category": "non-human VH", "sequence": VH},
        {"id": "mouse_kappa", "category": "kappa VL", "sequence": VL},
        {
            "id": "human_lambda",
            "category": "lambda VL",
            "sequence": synthetic_domain("L", "human"),
        },
        {
            "id": "cow_vh",
            "category": "non-human VH",
            "sequence": synthetic_domain("H", "cow"),
        },
        *[
            {
                "id": f"human_tcr_{chain.lower()}",
                "category": f"TCR {chain}",
                "sequence": synthetic_domain(chain, "human"),
            }
            for chain in "ABGD"
        ],
        {
            "id": "n_terminal_truncation",
            "category": "truncated domain",
            "sequence": VH[5:],
        },
        {
            "id": "c_terminal_truncation",
            "category": "truncated domain",
            "sequence": VH[:100],
        },
        {
            "id": "flanked_boundary",
            "category": "start/end boundary",
            "sequence": f"MPEPTIDE{VH}GG",
        },
        {
            "id": "cdr3_insertion",
            "category": "insertion",
            "sequence": VH.replace("CNAGHD", "CNAGAAAAHD"),
        },
        {
            "id": "cdr1_deletion",
            "category": "deletion",
            "sequence": VH[:29] + VH[33:],
        },
        {"id": "vl_vh_scfv", "category": "multiple domains", "sequence": SCFV},
        {"id": "lysozyme", "category": "non-antibody protein", "sequence": LYSOZYME},
    ]


def reference_domains(
    sequence: str,
    numbered: list[Any] | None,
    details: list[dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    if not numbered or not details:
        return []
    domains = []
    for domain_index, ((alignment, numbered_start, numbered_end), detail) in enumerate(
        zip(numbered, details)
    ):
        sequence_index = numbered_start
        residues = []
        for (position, insertion), amino_acid in alignment:
            if amino_acid == "-":
                continue
            residues.append(
                {
                    "sequenceIndex": sequence_index,
                    "aminoAcid": amino_acid,
                    "position": position,
                    "insertionCode": insertion.strip(),
                }
            )
            sequence_index += 1
        domains.append(
            {
                "domainIndex": domain_index,
                "profile": detail["id"],
                "chainType": detail["chain_type"],
                "species": detail["species"],
                "start": numbered_start,
                "end": numbered_end + 1,
                "bitScore": detail["bitscore"],
                "numbering": residues,
                "paddedImgtAlignment": "".join(amino for _, amino in alignment),
            }
        )
    return domains


def fnv1a(value: str) -> str:
    result = FNV_OFFSET
    for byte in value.encode("utf-8"):
        result ^= byte
        result = (result * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
    return f"{result:016x}"


def compact_domain(domain: dict[str, Any]) -> dict[str, Any]:
    numbering = "".join(
        f'{residue["sequenceIndex"]}|{residue["aminoAcid"]}|'
        f'{residue["position"]}|{residue["insertionCode"]}\n'
        for residue in domain["numbering"]
    )
    return {
        "chainType": domain["chainType"],
        "species": domain["species"],
        "start": domain["start"],
        "end": domain["end"],
        "bitScore": domain["bitScore"],
        "numberingLength": len(domain["numbering"]),
        "numberingFnv1a64": fnv1a(numbering),
        "paddedAlignmentFnv1a64": fnv1a(domain["paddedImgtAlignment"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--cases-file",
        type=Path,
        help="JSON array of {id, seq, category, note} cases instead of the built-in corpus",
    )
    parser.add_argument(
        "--compact",
        action="store_true",
        help="store hashes of alignments/numbering instead of the full reference vectors",
    )
    args = parser.parse_args()
    if version("anarci") != PINNED_ANARCI:
        raise SystemExit(f"this generator requires ANARCI {PINNED_ANARCI}")

    if args.cases_file:
        raw_cases = json.loads(args.cases_file.read_text(encoding="utf-8"))
        corpus_cases = [
            {
                "id": case["id"],
                "sequence": case["seq"],
                "category": case["category"],
                "note": case.get("note", ""),
            }
            for case in raw_cases
        ]
    else:
        corpus_cases = cases()
    numbered, details, _ = anarci(
        [(case["id"], case["sequence"]) for case in corpus_cases],
        scheme="imgt",
        allowed_species=None,
    )
    for case, case_numbered, case_details in zip(corpus_cases, numbered, details):
        case["referenceDomains"] = reference_domains(
            case["sequence"], case_numbered, case_details
        )

    hmmer_banner = subprocess.run(
        ["hmmscan", "-h"], check=True, capture_output=True, text=True
    ).stdout.splitlines()[1].strip()
    reference = {
        "anarciVersion": PINNED_ANARCI,
        "anarciCommit": ANARCI_COMMIT,
        "hmmer": hmmer_banner,
        "scheme": "imgt",
        "coordinates": "zero-based half-open numbered start/(inclusive end + 1)",
    }
    if args.compact:
        document = {
            "reference": reference,
            "hash": "FNV-1a 64 over canonical UTF-8 fields",
            "cases": [
                {
                    "id": case["id"],
                    "category": case["category"],
                    "domains": [
                        compact_domain(domain) for domain in case["referenceDomains"]
                    ],
                }
                for case in corpus_cases
            ],
        }
    else:
        document = {
            "reference": reference,
            "cases": corpus_cases,
            "pairs": [
                {"id": "valid", "vh": VH, "vl": VL, "ok": True},
                {"id": "swapped", "vh": VL, "vl": VH, "ok": False},
            ],
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {len(corpus_cases)} golden cases to {args.output}")


if __name__ == "__main__":
    main()
