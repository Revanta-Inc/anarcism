#!/usr/bin/env python3
"""Generate the checked-in ANARCI reference for the golden corpus.

Only this development tool invokes native ANARCI. JavaScript and Rust tests read
the resulting JSON and never need Python or HMMER.
"""

from __future__ import annotations

import argparse
import importlib
import json
import math
import subprocess
import tomllib
from hashlib import sha256
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "tests" / "golden"
MANIFEST = tomllib.loads((ROOT / "assets/MANIFEST.toml").read_text())
ANARCI_REPOSITORY = MANIFEST["reference"]["anarci_repository"]
ANARCI_COMMIT = MANIFEST["reference"]["anarci_commit"]
IMGT_GENEDB_PROGRAM_VERSION = MANIFEST["reference"][
    "imgt_genedb_program_version"
]
IMGT_SNAPSHOT = MANIFEST["reference"]["imgt_snapshot"]
EXPECTED_HMM_SHA256 = MANIFEST["profiles"]["source_sha256"]
EXPECTED_GERMLINES_SHA256 = MANIFEST["germlines"]["source_sha256"]
def anarci_with_score_components(
    sequences: list[tuple[str, str]],
) -> tuple[list[Any], list[Any], list[Any]]:
    """Run versioned ANARCI while retaining HMMER envelope coordinates.

    ANARCI's public details preserve the domain score and null2 bias,
    but discard the HSP envelope bounds needed to reconstruct the isolated
    Forward score. Its parser is wrapped only for this development tool.
    """
    anarci_module = importlib.import_module("anarci.anarci")
    original_parser = anarci_module._parse_hmmer_query

    def parse_with_envelopes(query: Any, *args: Any, **kwargs: Any) -> Any:
        result = original_parser(query, *args, **kwargs)
        details = result[2]
        used_hsps: set[int] = set()
        for detail in details:
            candidates = [
                (index, hsp)
                for index, hsp in enumerate(query.hsps)
                if index not in used_hsps
                and hsp.hit_id == detail["id"]
                and math.isclose(hsp.bitscore, detail["bitscore"], abs_tol=1.0e-6)
            ]
            if not candidates:
                raise RuntimeError(
                    f'could not recover HMMER envelope for {query.id}/{detail["id"]}'
                )
            index, hsp = min(
                candidates,
                key=lambda candidate: abs(
                    candidate[1].query_start - detail["query_start"]
                )
                + abs(candidate[1].query_end - detail["query_end"]),
            )
            used_hsps.add(index)
            detail["envelope_start"] = hsp.env_start
            detail["envelope_end"] = hsp.env_end
        return result

    anarci_module._parse_hmmer_query = parse_with_envelopes
    try:
        return anarci_module.anarci(
            sequences,
            scheme="imgt",
            allowed_species=None,
            assign_germline=True,
        )
    finally:
        anarci_module._parse_hmmer_query = original_parser


def verify_reference_inputs() -> None:
    anarci_module = importlib.import_module("anarci.anarci")
    germline_module = importlib.import_module("anarci.germlines")
    sources = [
        (
            Path(anarci_module.HMM_path) / "ALL.hmm",
            EXPECTED_HMM_SHA256,
            "ALL.hmm",
        ),
        (
            Path(germline_module.__file__ or ""),
            EXPECTED_GERMLINES_SHA256,
            "germlines.py",
        ),
    ]
    for source_path, expected_hash, label in sources:
        observed_hash = sha256(source_path.read_bytes()).hexdigest()
        if observed_hash != expected_hash:
            raise SystemExit(
                f"expected {IMGT_SNAPSHOT} {label} SHA-256 {expected_hash}, "
                f"found {observed_hash} at {source_path}"
            )


FORMAT = {
    "numbering": (
        "space-separated IMGT labels for residues start..end in sequence order; "
        "a-b is the run of plain positions a through b"
    ),
    "alternativeHits": ["profile", "bitScore", "eValue", "bias", "queryStart", "queryEnd"],
    "germline": ["species", "vGene", "vIdentity", "jGene", "jIdentity"],
    "eValueDisplaySignificantDigits": 2,
    "scoreDisplayPrecisionBits": 0.1,
}


def encode_numbering(labels: list[tuple[int, str]]) -> str:
    tokens = []
    index = 0
    while index < len(labels):
        position, insertion = labels[index]
        if insertion:
            tokens.append(f"{position}{insertion}")
            index += 1
            continue
        end = index
        while (
            end + 1 < len(labels)
            and not labels[end + 1][1]
            and labels[end + 1][0] == labels[end][0] + 1
        ):
            end += 1
        last = labels[end][0]
        tokens.append(str(position) if end == index else f"{position}-{last}")
        index = end + 1
    return " ".join(tokens)


def reference_domains(
    sequence: str,
    numbered: list[Any] | None,
    details: list[dict[str, Any]] | None,
    hit_table: list[list[Any]] | None,
) -> list[dict[str, Any]]:
    if not numbered or not details:
        return []
    domains = []
    for (alignment, numbered_start, numbered_end), detail in zip(numbered, details):
        labels = []
        sequence_index = numbered_start
        for (position, insertion), amino_acid in alignment:
            if amino_acid == "-":
                continue
            if sequence[sequence_index] != amino_acid:
                raise RuntimeError(
                    f"residue {sequence_index} is {amino_acid}, sequence has "
                    f"{sequence[sequence_index]}"
                )
            labels.append((position, insertion.strip()))
            sequence_index += 1
        # ANARCI reports an empty alignment when it cannot number the domain.
        if labels and sequence_index != numbered_end + 1:
            raise RuntimeError(
                f"numbered residues end at {sequence_index}, domain ends at {numbered_end + 1}"
            )
        domains.append(
            {
                "profile": detail["id"],
                "start": numbered_start,
                "end": numbered_end + 1,
                "bitScore": detail["bitscore"],
                "eValue": detail["evalue"],
                "bias": detail["bias"],
                "envelope": [detail["envelope_start"], detail["envelope_end"]],
                "alternativeHits": reference_alternative_hits(detail, hit_table),
                "germline": reference_germline(detail),
                "numbering": encode_numbering(labels),
                "alignment": "".join(amino for _, amino in alignment),
            }
        )
    return domains


def reference_alternative_hits(
    detail: dict[str, Any], hit_table: list[list[Any]] | None
) -> list[list[Any]]:
    """Return native ANARCI's first three overlapping non-winning hits."""
    if not hit_table:
        return []
    rows = [
        row
        for row in hit_table[1:]
        if row[5] < detail["query_end"] and detail["query_start"] < row[6]
    ]
    winner_index = next(
        (
            index
            for index, row in enumerate(rows)
            if row[0] == detail["id"]
            and math.isclose(row[3], detail["bitscore"], abs_tol=1.0e-6)
        ),
        None,
    )
    if winner_index is None:
        raise RuntimeError(f'could not find winning row for {detail["id"]}')
    del rows[winner_index]
    return [
        [profile, bit_score, e_value, bias, query_start, query_end]
        for profile, _description, e_value, bit_score, bias, query_start, query_end in rows[:3]
    ]


def reference_germline(detail: dict[str, Any]) -> list[Any] | None:
    germlines = detail.get("germlines")
    if not germlines:
        return None
    v_gene = germlines.get("v_gene")
    j_gene = germlines.get("j_gene")
    return [
        (v_gene or j_gene)[0][0],
        v_gene[0][1] if v_gene else None,
        round(v_gene[1], 6) if v_gene else None,
        j_gene[0][1] if j_gene else None,
        round(j_gene[1], 6) if j_gene else None,
    ]


def write_reference(path: Path, header: dict[str, Any], cases: list[dict[str, Any]]) -> None:
    """Write JSON Lines: a header record, then one record per case.

    Cases are sorted by id, so regenerating after corpus edits or reordering
    only changes the records whose values changed.
    """
    lines = [header, *sorted(cases, key=lambda case: case["id"])]
    path.write_text(
        "".join(
            json.dumps(line, separators=(",", ":"), ensure_ascii=False) + "\n"
            for line in lines
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cases-file",
        type=Path,
        default=GOLDEN / "corpus.json",
        help="JSON array of {id, seq, category, note} cases",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=GOLDEN / "corpus_reference.jsonl",
    )
    args = parser.parse_args()
    verify_reference_inputs()

    corpus = json.loads(args.cases_file.read_text(encoding="utf-8"))
    numbered, details, hit_tables = anarci_with_score_components(
        [(case["id"], case["seq"]) for case in corpus]
    )
    cases = [
        {
            "id": case["id"],
            "category": case["category"],
            "domains": reference_domains(
                case["seq"], case_numbered, case_details, case_hit_table
            ),
        }
        for case, case_numbered, case_details, case_hit_table in zip(
            corpus, numbered, details, hit_tables
        )
    ]

    hmmer_banner = subprocess.run(
        ["hmmscan", "-h"], check=True, capture_output=True, text=True
    ).stdout.splitlines()[1].strip()
    reference = {
        "anarciRepository": ANARCI_REPOSITORY,
        "anarciCommit": ANARCI_COMMIT,
        "hmmer": hmmer_banner,
        "imgtGenedbProgramVersion": IMGT_GENEDB_PROGRAM_VERSION,
        "imgtSnapshot": IMGT_SNAPSHOT,
        "hmmSourceSha256": EXPECTED_HMM_SHA256,
        "germlineSourceSha256": EXPECTED_GERMLINES_SHA256,
        "scheme": "imgt",
        "coordinates": "zero-based half-open [start, end) for domains, envelopes, and hits",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_reference(args.output, {"reference": reference, "format": FORMAT}, cases)
    print(f"wrote {len(cases)} reference cases to {args.output}")


if __name__ == "__main__":
    main()
