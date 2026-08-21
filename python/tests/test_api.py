"""Behaviour of the Python surface: options, errors, types, and GIL release."""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor

import pytest

import anarcism

VH = (
    "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTAD"
    "TSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA"
)
VL = (
    "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTL"
    "TISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR"
)


def test_module_metadata():
    assert anarcism.__version__
    assert anarcism.chains() == ["H", "K", "L", "A", "B", "G", "D"]
    assert "human" in anarcism.species()
    assert "mouse" in anarcism.species()


def test_number_sequence_defaults():
    result = anarcism.number_sequence(VH)
    assert result.id == "sequence"
    assert len(result.domains) == 1
    domain = result.domains[0]
    assert domain.chain_type == "H"
    assert domain.receptor_type == "IG"
    assert domain.start == 0
    # The trailing residue of this construct falls outside FR4, so the detected
    # domain deliberately stops short of the full input.
    assert domain.end == 120 < len(VH)
    assert len(domain.padded_imgt_alignment) >= 128
    # E-values are deliberately absent in this release.
    assert domain.e_value is None
    assert domain.germline is None


def test_custom_id_is_returned_and_used_in_errors():
    assert anarcism.number_sequence(VH, id="heavy").id == "heavy"
    with pytest.raises(anarcism.AnarcismError) as excinfo:
        anarcism.number_sequence("XXXX", id="bad-record")
    assert excinfo.value.input_id == "bad-record"


def test_regions_cover_the_imgt_boundaries():
    domain = anarcism.number_sequence(VH).domains[0]
    regions = {r.position: r.region for r in domain.numbering}
    assert regions[1] == "FR1"
    assert regions[27] == "CDR1"
    assert regions[39] == "FR2"
    assert regions[105] == "CDR3"
    assert regions[128] == "FR4"


def test_alternative_hit_count_bounds_the_hit_list():
    assert len(anarcism.number_sequence(VH, alternative_hit_count=1).domains[0].alternative_hits) == 1
    assert len(anarcism.number_sequence(VH, alternative_hit_count=0).domains[0].alternative_hits) == 0


def test_allowed_chains_is_strict_and_case_insensitive():
    assert anarcism.number_sequence(VH, allowed_chains=["K", "L"]).domains == []
    assert anarcism.number_sequence(VH, allowed_chains=["h"]).domains[0].chain_type == "H"


def test_allowed_species_is_strict():
    assert anarcism.number_sequence(VH, allowed_species=["human"]).domains[0].species == "human"


def test_min_bit_score_filters_domains():
    assert anarcism.number_sequence(VH, min_bit_score=1e6).domains == []


def test_assign_germline_populates_the_assignment():
    domain = anarcism.number_sequence(VH, assign_germline=True).domains[0]
    assert domain.germline is not None
    assert domain.germline.species == domain.species
    assert domain.germline.v_gene
    assert 0.0 <= domain.germline.v_identity <= 1.0


def test_number_sequences_preserves_order():
    results = anarcism.number_sequences([("a", VH), ("b", VL), ("c", VH)])
    assert [r.id for r in results] == ["a", "b", "c"]
    assert [r.domains[0].chain_type for r in results] == ["H", "K", "H"]


def test_number_fasta_reads_records():
    results = anarcism.number_fasta(f">heavy\n{VH}\n>light\n{VL}\n")
    assert [r.id for r in results] == ["heavy", "light"]


def test_validate_antibody_pair():
    pair = anarcism.validate_antibody_pair(VH, VL)
    assert pair.ok is True
    assert bool(pair) is True
    assert pair.errors == []
    assert pair.vh.chain_type == "H"
    assert pair.vl.chain_type in {"K", "L"}


def test_validate_antibody_pair_rejects_a_swapped_pair():
    pair = anarcism.validate_antibody_pair(VL, VH)
    assert pair.ok is False
    assert bool(pair) is False
    assert pair.errors


def test_results_are_immutable():
    domain = anarcism.number_sequence(VH).domains[0]
    with pytest.raises(AttributeError):
        domain.chain_type = "K"


def test_to_dict_is_json_serializable_and_snake_case():
    payload = anarcism.number_sequence(VH, assign_germline=True).to_dict()
    encoded = json.dumps(payload)
    assert "normalized_sequence" in payload
    assert "padded_imgt_alignment" in payload["domains"][0]
    assert "paddedImgtAlignment" not in encoded
    residue = payload["domains"][0]["numbering"][0]
    assert set(residue) == {
        "sequence_index",
        "amino_acid",
        "position",
        "insertion_code",
        "region",
    }


@pytest.mark.parametrize(
    ("call", "code"),
    [
        (lambda: anarcism.number_sequence("QQQ!QQQ"), "INVALID_SEQUENCE"),
        (lambda: anarcism.number_sequence("A" * 10_001), "SEQUENCE_TOO_LONG"),
        (lambda: anarcism.number_sequences([("x", VH)] * 1_001), "BATCH_TOO_LARGE"),
        (lambda: anarcism.number_fasta("not a fasta document"), "INVALID_FASTA"),
    ],
)
def test_errors_carry_a_stable_code(call, code):
    with pytest.raises(anarcism.AnarcismError) as excinfo:
        call()
    assert excinfo.value.code == code


def test_unknown_chain_is_a_value_error_not_an_engine_error():
    with pytest.raises(ValueError, match="unknown chain type"):
        anarcism.number_sequence(VH, allowed_chains=["Q"])


def test_workers_must_be_positive():
    with pytest.raises(ValueError, match="positive integer"):
        anarcism.number_sequences([("a", VH)], workers=0)


def test_gil_is_released_so_threads_overlap():
    """A pool must beat the serial path; if the GIL were held it would not."""
    inputs = [(f"s{i}", VH) for i in range(8)]

    def one(item):
        return anarcism.number_sequence(item[1], id=item[0])

    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(one, inputs))

    assert len(results) == 8
    assert all(r.domains[0].chain_type == "H" for r in results)
