"""Golden-corpus parity for the Python bindings.

The Rust golden test already proves the engine matches pinned ANARCI. This test
proves the binding layer transfers every value faithfully, by recomputing the
same FNV-1a 64 fingerprints from Python-visible attributes and comparing them
against the committed reference. The canonical residue form must stay identical
to `numbering_fingerprint` in `crates/anarcism-core/tests/golden.rs`.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import anarcism

GOLDEN = Path(__file__).resolve().parents[2] / "tests" / "golden"

FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x00000100000001B3
MASK = 0xFFFFFFFFFFFFFFFF


def fnv1a64(data: bytes) -> str:
    digest = FNV_OFFSET
    for byte in data:
        digest ^= byte
        digest = (digest * FNV_PRIME) & MASK
    return f"{digest:016x}"


def numbering_fingerprint(domain: anarcism.DomainResult) -> str:
    canonical = "".join(
        f"{r.sequence_index}|{r.amino_acid}|{r.position}|{r.insertion_code}\n"
        for r in domain.numbering
    )
    return fnv1a64(canonical.encode())


def load(name: str):
    return json.loads((GOLDEN / name).read_text())


@pytest.fixture(scope="module")
def corpus() -> list[tuple[dict, dict]]:
    cases = {case["id"]: case for case in load("corpus_v2.json")}
    reference = load("corpus_v2_reference.json")["cases"]
    missing = [case["id"] for case in reference if case["id"] not in cases]
    assert not missing, f"reference ids absent from corpus_v2.json: {missing[:5]}"
    return [(cases[case["id"]], case) for case in reference]


def test_corpus_is_not_empty(corpus):
    assert len(corpus) > 1_000


def test_every_domain_matches_the_reference(corpus):
    inputs = [(case["id"], case["seq"]) for case, _ in corpus]
    # The engine caps a batch at 1,000 records, so the corpus is chunked.
    results = []
    for start in range(0, len(inputs), 1_000):
        results.extend(anarcism.number_sequences(inputs[start : start + 1_000]))
    assert len(results) == len(corpus)

    residues_checked = 0
    domains_checked = 0

    for (case, expected), result in zip(corpus, results, strict=True):
        assert result.id == case["id"]
        assert len(result.domains) == len(expected["domains"]), (
            f"{case['id']}: domain count"
        )

        for observed, want in zip(result.domains, expected["domains"], strict=True):
            label = f"{case['id']} domain {observed.domain_index}"
            assert observed.chain_type == want["chainType"], f"{label}: chain"
            assert observed.species == want["species"], f"{label}: species"
            assert observed.start == want["start"], f"{label}: start"
            assert observed.end == want["end"], f"{label}: end"
            assert len(observed.numbering) == want["numberingLength"], (
                f"{label}: numbering length"
            )
            assert numbering_fingerprint(observed) == want["numberingFnv1a64"], (
                f"{label}: numbering fingerprint"
            )
            assert (
                fnv1a64(observed.padded_imgt_alignment.encode())
                == want["paddedAlignmentFnv1a64"]
            ), f"{label}: padded alignment fingerprint"
            # Scores are quantized, so the reference records one decimal place.
            assert observed.bit_score == pytest.approx(want["bitScore"], abs=1.0), (
                f"{label}: bit score"
            )
            residues_checked += len(observed.numbering)
            domains_checked += 1

    assert domains_checked == sum(len(c["domains"]) for _, c in corpus)
    assert residues_checked > 100_000


def test_workers_match_the_serial_path(corpus):
    inputs = [(case["id"], case["seq"]) for case, _ in corpus[:200]]
    serial = anarcism.number_sequences(inputs)
    parallel = anarcism.number_sequences(inputs, workers=4)
    assert [r.to_dict() for r in parallel] == [r.to_dict() for r in serial]
