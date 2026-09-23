from __future__ import annotations

import json
from pathlib import Path

import pytest

import anarcism

GOLDEN = Path(__file__).resolve().parents[2] / "tests" / "golden"


def encode_numbering(case: dict, domain: anarcism.DomainResult) -> str:
    """Encode numbering like the reference: labels in residue order, `a-b` runs."""
    tokens: list[str] = []
    run: list[int] | None = None
    for offset, residue in enumerate(domain.numbering):
        assert residue.sequence_index == domain.start + offset, case["id"]
        assert residue.amino_acid == case["seq"][residue.sequence_index], case["id"]
        if not residue.insertion_code and run and residue.position == run[1] + 1:
            run[1] = residue.position
            continue
        if run:
            tokens.append(str(run[0]) if run[0] == run[1] else f"{run[0]}-{run[1]}")
            run = None
        if residue.insertion_code:
            tokens.append(f"{residue.position}{residue.insertion_code}")
        else:
            run = [residue.position, residue.position]
    if run:
        tokens.append(str(run[0]) if run[0] == run[1] else f"{run[0]}-{run[1]}")
    return " ".join(tokens)


def load_reference() -> list[dict]:
    lines = (GOLDEN / "corpus_reference.jsonl").read_text().splitlines()
    return [json.loads(line) for line in lines[1:]]


@pytest.fixture(scope="module")
def corpus() -> list[tuple[dict, dict]]:
    cases = {case["id"]: case for case in json.loads((GOLDEN / "corpus.json").read_text())}
    reference = load_reference()
    missing = [case["id"] for case in reference if case["id"] not in cases]
    assert not missing, f"reference ids absent from corpus.json: {missing[:5]}"
    return [(cases[case["id"]], case) for case in reference]


def test_corpus_is_not_empty(corpus):
    assert len(corpus) > 1_000


def test_every_domain_matches_the_reference(corpus):
    inputs = [(case["id"], case["seq"]) for case, _ in corpus]
    results = anarcism.number_sequences(inputs, assign_germline=True)
    assert len(results) == len(corpus)

    residues_checked = 0
    domains_checked = 0

    for (case, expected), result in zip(corpus, results, strict=True):
        assert result.id == case["id"]
        assert len(result.domains) == len(
            expected["domains"]
        ), f"{case['id']}: domain count"

        for observed, want in zip(result.domains, expected["domains"], strict=True):
            label = f"{case['id']} domain {observed.domain_index}"
            profile = f"{observed.species}_{observed.chain_type}"
            assert profile == want["profile"], f"{label}: profile"
            assert observed.start == want["start"], f"{label}: start"
            assert observed.end == want["end"], f"{label}: end"
            assert (
                encode_numbering(case, observed) == want["numbering"]
            ), f"{label}: numbering"
            assert (
                observed.padded_imgt_alignment == want["alignment"]
            ), f"{label}: padded alignment"
            assert observed.bit_score == pytest.approx(
                want["bitScore"], abs=0.2
            ), f"{label}: bit score"
            assert (
                abs(observed.e_value - want["eValue"]) / want["eValue"] <= 0.08
            ), f"{label}: E-value"
            assert observed.bias == pytest.approx(want["bias"], abs=0.2), f"{label}: bias"
            germline = observed.germline
            expected_germline = want["germline"]
            assert (germline is None) == (
                expected_germline is None
            ), f"{label}: germline presence"
            if germline is not None:
                species, v_gene, v_identity, j_gene, j_identity = expected_germline
                assert germline.species == species
                assert germline.v_gene == v_gene
                assert germline.j_gene == j_gene
                assert germline.v_identity == pytest.approx(v_identity, abs=0.08)
                assert germline.j_identity == pytest.approx(j_identity, abs=0.08)
            residues_checked += len(observed.numbering)
            domains_checked += 1

    assert domains_checked == sum(len(c["domains"]) for _, c in corpus)
    assert residues_checked > 100_000


def test_alternative_hits_match_versioned_anarci_across_chain_families(corpus):
    strict_ids = {
        "trastuzumab_vh",
        "trastuzumab_vl",
        "human_trav12_2_traj33",
        "human_trbv19_trbj2_7",
        "human_trgv9_trgjp",
        "human_trdv2_trdj1",
    }
    selected = [
        (case, reference) for case, reference in corpus if case["id"] in strict_ids
    ]
    assert len(selected) == len(strict_ids)

    for case, reference in selected:
        result = anarcism.number_sequence(
            case["seq"],
            id=case["id"],
            alternative_hit_count=28,
        )
        for observed, expected in zip(
            result.domains, reference["domains"], strict=True
        ):
            for actual_hit, expected_hit in zip(
                observed.alternative_hits, expected["alternativeHits"], strict=False
            ):
                profile, bit_score, e_value, bias, query_start, query_end = expected_hit
                label = f"{case['id']} alternative {profile}"
                assert actual_hit.profile == profile, label
                assert f"{actual_hit.species}_{actual_hit.chain_type}" == profile, label
                assert actual_hit.bit_score == pytest.approx(bit_score, abs=0.2), label
                assert abs(actual_hit.e_value - e_value) / e_value <= 0.08, label
                assert actual_hit.bias == pytest.approx(bias, abs=0.2), label
                assert actual_hit.query_start == query_start, label
                assert actual_hit.query_end == query_end, label


def test_workers_match_the_serial_path(corpus):
    inputs = [(case["id"], case["seq"]) for case, _ in corpus[:200]]
    serial = anarcism.number_sequences(inputs)
    parallel = anarcism.number_sequences(inputs, workers=4)
    assert [r.to_dict() for r in parallel] == [r.to_dict() for r in serial]
