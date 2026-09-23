from __future__ import annotations

import re
from io import StringIO
from pathlib import Path

import pytest

from anarcism.cli import main

FIXTURES = Path(__file__).with_name("fixtures")

VH = (
    "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTAD"
    "TSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA"
)
VL = (
    "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTL"
    "TISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR"
)
HIT_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTIS"
    "ADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
LONG_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTIS"
    "ADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGYYYDSSGYYYRGSGWGGDGFYAMDYWGQGTLVTVSS"
)


def test_vertical_output_matches_anarci_layout(capsys):
    assert main(["-i", VH]) == 0
    captured = capsys.readouterr()
    assert captured.err == ""
    lines = captured.out.splitlines()
    assert lines[:5] == [
        "# Input sequence",
        "# ANARCI numbered",
        "# Domain 1 of 1",
        "# Most significant HMM hit",
        "#|species|chain_type|e-value|score|seqstart_index|seqend_index|",
    ]
    assert re.fullmatch(
        r"#\|(human|mouse)\|H\|[0-9.]+e[+-]?[0-9]+\|[0-9]+\.[0-9]\|0\|119\|",
        lines[5],
    )
    assert lines[6] == "# Scheme = imgt"
    assert "H 1       E" in lines
    assert "H 10      -" in lines
    assert lines[-2] == "H 128     A"
    assert lines[-1] == "//"


def test_no_domain_output_is_exact(capsys):
    assert main(["-i", "ACDEFGHIKLMNPQRSTVWY"]) == 0
    captured = capsys.readouterr()
    assert captured.out == "# Input sequence\n//\n"
    assert captured.err == ""


def test_fasta_ids_and_outfile(tmp_path: Path, capsys):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(f">heavy record\n{VH}\n>negative\nACDEFGHIKLMNPQRSTVWY\n")
    output = tmp_path / "numbered.anarci"

    assert main(["-i", str(fasta), "-o", str(output), "-p", "2"]) == 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
    rendered = output.read_text()
    assert rendered.startswith("# heavy record\n# ANARCI numbered\n")
    assert rendered.endswith("# negative\n//\n")


def test_restrict_aliases_filter_chain_types(capsys):
    assert main(["-i", VH, "-r", "light"]) == 0
    assert capsys.readouterr().out == "# Input sequence\n//\n"

    assert main(["-i", VH, "-r", "ig", "-s", "i"]) == 0
    assert "#|mouse|H|" in capsys.readouterr().out


def test_bit_score_threshold_filters_domains(capsys):
    assert main(["-i", VH, "--bit_score_threshold", "1000"]) == 0
    assert capsys.readouterr().out == "# Input sequence\n//\n"


def test_fasta_stdin(monkeypatch, capsys):
    monkeypatch.setattr("sys.stdin", StringIO(f">stdin-heavy\n{VH}\n"))
    assert main(["-i", "-"]) == 0
    assert capsys.readouterr().out.startswith("# stdin-heavy\n# ANARCI numbered\n")


def test_batch_boundaries_do_not_change_vertical_or_hit_output(
    tmp_path: Path, capsys, monkeypatch
):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(
        f">heavy\n{HIT_VH}\n>long-heavy\n{LONG_VH}\n>negative\nACDEFGHIKLMNPQRSTVWY\n"
    )
    streamed = tmp_path / "streamed.anarci"
    streamed_hits = tmp_path / "streamed.hits"
    combined = tmp_path / "combined.anarci"
    combined_hits = tmp_path / "combined.hits"

    monkeypatch.setattr("anarcism.cli.DEFAULT_BATCH_SIZE", 1)
    monkeypatch.setattr("anarcism.cli.DEFAULT_BATCH_RESIDUES", 1)
    assert (
        main(
            [
                "-i",
                str(fasta),
                "-o",
                str(streamed),
                "-ht",
                str(streamed_hits),
            ]
        )
        == 0
    )
    monkeypatch.setattr("anarcism.cli.DEFAULT_BATCH_SIZE", 100)
    monkeypatch.setattr("anarcism.cli.DEFAULT_BATCH_RESIDUES", 100_000)
    assert (
        main(
            [
                "-i",
                str(fasta),
                "-o",
                str(combined),
                "-ht",
                str(combined_hits),
            ]
        )
        == 0
    )
    assert capsys.readouterr().err == ""
    assert streamed.read_bytes() == combined.read_bytes()
    assert streamed_hits.read_bytes() == combined_hits.read_bytes()


def test_germline_section_matches_anarci_layout(capsys):
    assert main(["-i", VH, "--assign_germline", "--use_species", "mouse"]) == 0
    lines = capsys.readouterr().out.splitlines()
    index = lines.index("# Most sequence-identical germlines")
    assert lines[index + 1] == "#|species|v_gene|v_identity|j_gene|j_identity|"
    assert re.fullmatch(
        r"#\|mouse\|[^|]+\|0\.\d{2}\|[^|]+\|0\.\d{2}\|", lines[index + 2]
    )


def test_csv_output_matches_anarci_files_and_columns(tmp_path: Path, capsys):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(f">heavy,one\n{VH}\n>light\n{VL}\n")
    output_root = tmp_path / "numbered"

    assert main(["-i", str(fasta), "--csv", "-o", str(output_root)]) == 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""

    heavy = Path(f"{output_root}_H.csv")
    light = Path(f"{output_root}_KL.csv")
    assert heavy.is_file()
    assert light.is_file()
    header, row = heavy.read_text().splitlines()
    assert header.startswith(
        "Id,domain_no,hmm_species,chain_type,e-value,score,seqstart_index,"
        "seqend_index,identity_species,v_gene,v_identity,j_gene,j_identity,1,2,3"
    )
    fields = row.split(",")
    assert fields[0] == "heavy one"
    assert fields[1:4] == ["0", "mouse", "H"]
    assert re.fullmatch(r"[0-9.]+e[+-]?[0-9]+", fields[4])
    assert re.fullmatch(r"[0-9]+\.[0-9]", fields[5])
    assert fields[6:8] == ["0", "119"]
    assert fields[8:13] == ["", "", "0.00", "", "0.00"]
    assert light.read_text().splitlines()[1].split(",")[3] in {"K", "L"}


def test_csv_schema_is_accumulated_across_streamed_batches(
    tmp_path: Path, capsys, monkeypatch
):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(f">ordinary\n{HIT_VH}\n>long-cdr3\n{LONG_VH}\n>light\n{VL}\n")
    streamed_root = tmp_path / "streamed"
    combined_root = tmp_path / "combined"

    monkeypatch.setattr("anarcism.cli.DEFAULT_BATCH_SIZE", 1)
    assert (
        main(
            [
                "-i",
                str(fasta),
                "--csv",
                "-o",
                str(streamed_root),
            ]
        )
        == 0
    )
    monkeypatch.setattr("anarcism.cli.DEFAULT_BATCH_SIZE", 100)
    assert (
        main(
            [
                "-i",
                str(fasta),
                "--csv",
                "-o",
                str(combined_root),
            ]
        )
        == 0
    )
    assert capsys.readouterr().err == ""
    for chain_class in ("H", "KL"):
        streamed = Path(f"{streamed_root}_{chain_class}.csv")
        combined = Path(f"{combined_root}_{chain_class}.csv")
        assert streamed.read_bytes() == combined.read_bytes()
    assert any(
        re.fullmatch(r"\d+[A-Z]+", field)
        for field in Path(f"{streamed_root}_H.csv")
        .read_text()
        .splitlines()[0]
        .split(",")
    )


def test_csv_requires_an_output_root(capsys):
    assert main(["-i", VH, "--csv"]) == 1
    assert capsys.readouterr().err == (
        "Error: When --csv option is used an ouput file name must be given.\n"
    )


def test_hit_output_matches_versioned_anarci_layout(tmp_path: Path, capsys):
    numbered = tmp_path / "numbered.anarci"
    hits = tmp_path / "hits.txt"

    assert main(["-i", HIT_VH, "-o", str(numbered), "-ht", str(hits)]) == 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""

    lines = hits.read_text().splitlines()
    assert lines[:5] == [
        "# Hit file for ANARCI",
        "NAME     Input sequence",
        f"SEQUENCE {HIT_VH[:71]}",
        f"SEQUENCE {HIT_VH[71:]}",
        "         id description      evalue    bitscore        bias query_start   query_end",
    ]
    rows = [line.split() for line in lines[5:-1]]
    assert rows == [
        ["human_H", "2.1e-60", "193.6", "0.5", "0", "120"],
        ["mouse_H", "2.7e-56", "180.2", "0.1", "0", "120"],
    ]
    assert lines[-1] == "//"


def test_versioned_anarci_output_files_match_byte_for_byte(tmp_path: Path, capsys):
    vertical = tmp_path / "numbered.anarci"
    csv_root = tmp_path / "numbered"
    hits = tmp_path / "hits.txt"
    hit_numbering = tmp_path / "hit-numbered.anarci"

    assert (
        main(
            [
                "-i",
                VH,
                "-o",
                str(vertical),
                "--assign_germline",
                "--use_species",
                "mouse",
            ]
        )
        == 0
    )
    assert (
        main(
            [
                "-i",
                VH,
                "--csv",
                "-o",
                str(csv_root),
                "--assign_germline",
                "--use_species",
                "mouse",
            ]
        )
        == 0
    )
    assert main(["-i", HIT_VH, "-o", str(hit_numbering), "-ht", str(hits)]) == 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""

    assert vertical.read_bytes() == (FIXTURES / "anarci_vertical.txt").read_bytes()
    assert (
        Path(f"{csv_root}_H.csv").read_bytes()
        == (FIXTURES / "anarci_H.csv").read_bytes()
    )
    assert hits.read_bytes() == (FIXTURES / "anarci_hits.txt").read_bytes()


def test_hit_output_precedes_chain_restriction(tmp_path: Path, capsys):
    numbered = tmp_path / "numbered.anarci"
    hits = tmp_path / "hits.txt"

    assert main(["-i", VH, "-r", "light", "-o", str(numbered), "-ht", str(hits)]) == 0
    assert capsys.readouterr().err == ""
    assert numbered.read_text() == "# Input sequence\n//\n"
    assert "mouse_H" in hits.read_text()


def test_missing_hit_output_directory_matches_anarci_error(capsys, tmp_path: Path):
    hits = tmp_path / "missing" / "hits.txt"
    assert main(["-i", VH, "-ht", str(hits)]) == 1
    assert capsys.readouterr().err == "Error: Hit output file path does not exist\n"


def test_missing_output_directory_matches_anarci_error(capsys, tmp_path: Path):
    output = tmp_path / "missing" / "numbered.anarci"
    assert main(["-i", VH, "-o", str(output)]) == 1
    assert capsys.readouterr().err == "Error: Output file path does not exist\n"


def test_late_fasta_error_does_not_replace_output_files(
    tmp_path: Path, capsys, monkeypatch
):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(f">valid\n{VH}\n>invalid\n!!!!\n")
    numbered = tmp_path / "numbered.anarci"
    hits = tmp_path / "hits.txt"
    numbered.write_text("existing numbered output\n")
    hits.write_text("existing hit output\n")

    monkeypatch.setattr("anarcism.cli.DEFAULT_BATCH_SIZE", 1)
    assert (
        main(
            [
                "-i",
                str(fasta),
                "-o",
                str(numbered),
                "-ht",
                str(hits),
            ]
        )
        == 1
    )
    assert "unsupported residue symbol '!'" in capsys.readouterr().err
    assert numbered.read_text() == "existing numbered output\n"
    assert hits.read_text() == "existing hit output\n"


def test_late_fasta_error_does_not_replace_csv(tmp_path: Path, capsys, monkeypatch):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(f">valid\n{VH}\n>invalid\n!!!!\n")
    output_root = tmp_path / "numbered"
    heavy = Path(f"{output_root}_H.csv")
    heavy.write_text("existing CSV output\n")

    monkeypatch.setattr("anarcism.cli.DEFAULT_BATCH_SIZE", 1)
    assert (
        main(
            [
                "-i",
                str(fasta),
                "--csv",
                "-o",
                str(output_root),
            ]
        )
        == 1
    )
    assert "unsupported residue symbol '!'" in capsys.readouterr().err
    assert heavy.read_text() == "existing CSV output\n"


def test_non_imgt_schemes_are_rejected():
    with pytest.raises(SystemExit, match="2"):
        main(["-i", VH, "-s", "kabat"])


def test_zero_workers_are_rejected():
    with pytest.raises(SystemExit, match="2"):
        main(["-i", VH, "-p", "0"])


@pytest.mark.parametrize("option", ["--batch_size", "--batch_residues"])
def test_batch_options_are_not_exposed(option):
    with pytest.raises(SystemExit, match="2"):
        main(["-i", VH, option, "1"])
