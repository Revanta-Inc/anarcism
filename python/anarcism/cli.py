"""ANARCI-compatible command-line interface for the embedded IMGT backend."""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import textwrap
from collections.abc import Iterable, Iterator, Sequence
from contextlib import ExitStack, contextmanager, nullcontext
from pathlib import Path
from typing import TextIO

from . import AnarcismError, DomainResult, SequenceResult, number_sequences, species
from .streaming import (
    DEFAULT_BATCH_RESIDUES,
    DEFAULT_BATCH_SIZE,
    iter_fasta,
    iter_sequence_batches,
)

_CHAIN_ORDER = ("H", "K", "L", "A", "B", "G", "D")
_CHAIN_GROUPS = {
    "ig": ("H", "K", "L"),
    "tr": ("A", "B", "G", "D"),
    "heavy": ("H",),
    "light": ("K", "L"),
}
_CSV_CHAIN_ORDER = ("H", "KL", "A", "B", "G", "D")
_LIGHT_CHAINS = {"K", "L"}

_DESCRIPTION = r"""
ANARCI                                                 \\    //
Antibody Numbering and Antigen Receptor ClassIfication  \\  //
                                                          ||
                                                          WA

ANARCI-compatible, embedded IMGT numbering. No HMMER executable or external
model files are required.
"""

_EPILOGUE = """Only the IMGT numbering scheme is supported by this backend."""


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ANARCI",
        description=_DESCRIPTION,
        epilog=_EPILOGUE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--sequence",
        "-i",
        dest="inputsequence",
        help="A sequence, an input FASTA file, or - for FASTA on stdin",
    )
    parser.add_argument(
        "--outfile",
        "-o",
        default=False,
        help="The output file to use. Default is stdout",
    )
    parser.add_argument(
        "--scheme",
        "-s",
        choices=("imgt", "i"),
        default="imgt",
        help=(
            "Numbering scheme. This backend supports IMGT (or shorthand i). "
            "Default IMGT"
        ),
    )
    parser.add_argument(
        "--restrict",
        "-r",
        nargs="+",
        choices=("ig", "tr", "heavy", "light", "H", "K", "L", "A", "B"),
        default=False,
        help="Restrict recognition to specified receptor chain types.",
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        default=False,
        help=(
            "Write CSV output. Outfile must be specified; one file is written "
            "per chain type and kappa/lambda are combined."
        ),
    )
    parser.add_argument(
        "--outfile_hits",
        "-ht",
        dest="hitfile",
        default=False,
        help="Output file for domain hit tables for each sequence.",
    )
    parser.add_argument(
        "--ncpu",
        "-p",
        type=int,
        default=1,
        help="Number of native worker threads to use. Default is 1.",
    )
    parser.add_argument(
        "--assign_germline",
        action="store_true",
        default=False,
        help="Assign the most sequence-identical V and J germlines.",
    )
    parser.add_argument(
        "--use_species",
        choices=tuple(species()),
        help=(
            "Restrict profile recognition and germline assignment to one "
            "species. By default human and mouse are considered."
        ),
    )
    parser.add_argument(
        "--bit_score_threshold",
        type=int,
        default=80,
        help="Change the bit-score threshold used to accept a domain.",
    )
    return parser


def _allowed_chains(restrictions: Sequence[str] | bool) -> list[str] | None:
    if not restrictions:
        return None
    selected: set[str] = set()
    for restriction in restrictions:
        selected.update(_CHAIN_GROUPS.get(restriction, (restriction,)))
    return [chain for chain in _CHAIN_ORDER if chain in selected]


@contextmanager
def _input_records(value: str) -> Iterator[Iterable[tuple[str, str]]]:
    if value == "-":
        yield iter_fasta(sys.stdin)
        return

    path = Path(value)
    if path.is_file():
        with path.open("r", encoding="utf-8") as fasta:
            yield iter_fasta(fasta)
        return

    yield iter((("Input sequence", value),))


def _number_records(
    records: Sequence[tuple[str, str]],
    *,
    workers: int,
    allowed_chains: list[str] | None,
    allowed_species: list[str],
    min_bit_score: float,
    alternative_hit_count: int = 3,
    assign_germline: bool,
) -> list[SequenceResult]:
    return number_sequences(
        list(records),
        workers=workers,
        allowed_chains=allowed_chains,
        allowed_species=allowed_species,
        min_bit_score=min_bit_score,
        alternative_hit_count=alternative_hit_count,
        assign_germline=assign_germline,
    )


def _format_e_value(value: float) -> str:
    # Match HMMER's displayed precision.
    return format(value, ".2g")


def _alignment(domain: DomainResult) -> list[tuple[tuple[int, str], str]]:
    if not domain.numbering:
        return []

    alignment: list[tuple[tuple[int, str], str]] = []
    next_position = 1
    for residue in domain.numbering:
        while next_position < residue.position:
            alignment.append(((next_position, ""), "-"))
            next_position += 1
        alignment.append(
            ((residue.position, residue.insertion_code), residue.amino_acid)
        )
        if not residue.insertion_code:
            next_position = max(next_position, residue.position + 1)

    while next_position <= 117:
        alignment.append(((next_position, ""), "-"))
        next_position += 1

    assert "".join(amino_acid for _, amino_acid in alignment) == (
        domain.padded_imgt_alignment
    ), "backend numbering and padded alignment disagree"
    return alignment


def _print_germline(domain: DomainResult, outfile: TextIO) -> None:
    germline = domain.germline
    species_name = germline.species if germline is not None else ""
    v_gene = germline.v_gene if germline is not None else None
    v_identity = germline.v_identity if germline is not None else None
    j_gene = germline.j_gene if germline is not None else None
    j_identity = germline.j_identity if germline is not None else None
    print("# Most sequence-identical germlines", file=outfile)
    print("#|species|v_gene|v_identity|j_gene|j_identity|", file=outfile)
    print(
        f"#|{species_name}|{v_gene or 'unknown'}|{v_identity or 0.0:.2f}|"
        f"{j_gene or 'unknown'}|{j_identity or 0.0:.2f}|",
        file=outfile,
    )


def write_anarci_output(
    results: Iterable[SequenceResult],
    outfile: TextIO,
    *,
    assign_germline: bool = False,
) -> None:
    """Write ANARCI's vertical numbering format."""
    for result in results:
        print(f"# {result.id}", file=outfile)
        if result.domains:
            print("# ANARCI numbered", file=outfile)
            for domain_number, domain in enumerate(result.domains, start=1):
                print(
                    f"# Domain {domain_number} of {len(result.domains)}", file=outfile
                )
                print("# Most significant HMM hit", file=outfile)
                print(
                    "#|species|chain_type|e-value|score|seqstart_index|seqend_index|",
                    file=outfile,
                )
                print(
                    f"#|{domain.species}|{domain.chain_type}|"
                    f"{_format_e_value(domain.e_value)}|{domain.bit_score:.1f}|"
                    f"{domain.start}|{domain.end - 1}|",
                    file=outfile,
                )
                if assign_germline:
                    _print_germline(domain, outfile)
                print("# Scheme = imgt", file=outfile)
                if not domain.numbering:
                    print(
                        "# Warning: imgt scheme could not be applied to this sequence.",
                        file=outfile,
                    )
                chain_class = (
                    "L" if domain.chain_type in _LIGHT_CHAINS else domain.chain_type
                )
                for (position, insertion), amino_acid in _alignment(domain):
                    print(
                        chain_class,
                        str(position).ljust(5),
                        insertion or " ",
                        amino_acid,
                        file=outfile,
                    )
        print("//", file=outfile)


def _csv_metadata(result: SequenceResult, domain: DomainResult) -> list[str]:
    germline = domain.germline
    identity_species = germline.species if germline is not None else ""
    v_gene = germline.v_gene if germline is not None else None
    v_identity = germline.v_identity if germline is not None else None
    j_gene = germline.j_gene if germline is not None else None
    j_identity = germline.j_identity if germline is not None else None
    return [
        result.id.replace(",", " "),
        str(domain.domain_index),
        domain.species,
        domain.chain_type,
        _format_e_value(domain.e_value),
        f"{domain.bit_score:.1f}",
        str(domain.start),
        str(domain.end - 1),
        identity_species,
        v_gene or "",
        f"{v_identity or 0.0:.2f}",
        j_gene or "",
        f"{j_identity or 0.0:.2f}",
    ]


_CSV_METADATA_FIELDS = [
    "Id",
    "domain_no",
    "hmm_species",
    "chain_type",
    "e-value",
    "score",
    "seqstart_index",
    "seqend_index",
    "identity_species",
    "v_gene",
    "v_identity",
    "j_gene",
    "j_identity",
]


@contextmanager
def _atomic_text_output(path: Path) -> Iterator[TextIO]:
    temporary = tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    )
    temporary_path = Path(temporary.name)
    try:
        with temporary:
            yield temporary
        os.replace(temporary_path, path)
    except BaseException:
        temporary.close()
        temporary_path.unlink(missing_ok=True)
        raise


class _CsvOutputSpool:
    """Collect a dynamic CSV schema while keeping rows in temporary files."""

    def __init__(self, output_root: str) -> None:
        self.output_root = output_root
        self.directory = Path(output_root).parent
        self.spools: dict[str, TextIO] = {}
        self.ranks: dict[str, dict[tuple[int, str], int]] = {}
        self.positions: dict[str, set[tuple[int, str]]] = {}

    def __enter__(self) -> _CsvOutputSpool:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def close(self) -> None:
        for spool in self.spools.values():
            spool.close()

    def write(self, results: Iterable[SequenceResult]) -> None:
        for result in results:
            for domain in result.domains:
                chain_class = (
                    "KL" if domain.chain_type in _LIGHT_CHAINS else domain.chain_type
                )
                alignment = _alignment(domain)
                self._record_positions(chain_class, alignment)
                spool = self.spools.get(chain_class)
                if spool is None:
                    spool = tempfile.TemporaryFile(
                        "w+",
                        encoding="utf-8",
                        newline="",
                        dir=self.directory,
                    )
                    self.spools[chain_class] = spool
                payload = [
                    _csv_metadata(result, domain),
                    [
                        [position, insertion, amino_acid]
                        for (position, insertion), amino_acid in alignment
                    ],
                ]
                json.dump(payload, spool, separators=(",", ":"))
                spool.write("\n")

    def _record_positions(
        self,
        chain_class: str,
        alignment: list[tuple[tuple[int, str], str]],
    ) -> None:
        chain_ranks = self.ranks.setdefault(chain_class, {})
        chain_positions = self.positions.setdefault(chain_class, set())
        last_position = -1
        rank = 0
        for position, _ in alignment:
            if position[0] != last_position:
                last_position = position[0]
                rank = 0
            else:
                rank += 1
            chain_ranks[position] = max(rank, chain_ranks.get(position, rank))
            chain_positions.add(position)

    def finish(self) -> None:
        for chain_class in _CSV_CHAIN_ORDER:
            spool = self.spools.get(chain_class)
            if spool is None:
                continue
            chain_positions = sorted(
                self.positions[chain_class],
                key=lambda position: (
                    position[0],
                    self.ranks[chain_class][position],
                ),
            )
            path = Path(f"{self.output_root}_{chain_class}.csv")
            with _atomic_text_output(path) as outfile:
                fields = _CSV_METADATA_FIELDS + [
                    f"{position}{insertion}" for position, insertion in chain_positions
                ]
                print(",".join(fields), file=outfile)
                spool.seek(0)
                for serialized in spool:
                    metadata, serialized_alignment = json.loads(serialized)
                    numbered = {
                        (position, insertion): amino_acid
                        for position, insertion, amino_acid in serialized_alignment
                    }
                    metadata.extend(
                        numbered.get(position, "-") for position in chain_positions
                    )
                    print(",".join(metadata), file=outfile)


def write_csv_output(results: Iterable[SequenceResult], output_root: str) -> None:
    """Write ANARCI's one-file-per-chain CSV representation."""
    with _CsvOutputSpool(output_root) as output:
        output.write(results)
        output.finish()


_HIT_TABLE_HEADER = (
    "id",
    "description",
    "evalue",
    "bitscore",
    "bias",
    "query_start",
    "query_end",
)


def _hit_table_rows(result: SequenceResult) -> list[tuple[float, list[str]]]:
    rows: list[tuple[float, list[str]]] = []
    for domain in result.domains:
        rows.append(
            (
                domain.e_value,
                [
                    f"{domain.species}_{domain.chain_type}",
                    "",
                    _format_e_value(domain.e_value),
                    f"{domain.bit_score:.1f}",
                    f"{domain.bias:.1f}",
                    str(domain.query_start),
                    str(domain.query_end),
                ],
            )
        )
        for hit in domain.alternative_hits:
            rows.append(
                (
                    hit.e_value,
                    [
                        hit.profile,
                        "",
                        _format_e_value(hit.e_value),
                        f"{hit.bit_score:.1f}",
                        f"{hit.bias:.1f}",
                        str(hit.query_start),
                        str(hit.query_end),
                    ],
                )
            )
    rows.sort(key=lambda row: (row[0], row[1][0], int(row[1][-2])))
    return rows


def _write_hit_header(outfile: TextIO) -> None:
    print("# Hit file for ANARCI", file=outfile)


def _write_hit_results(
    results: Iterable[SequenceResult], outfile: TextIO
) -> None:
    padding = max(len(field) for field in _HIT_TABLE_HEADER)
    for result in results:
        print("NAME    ", result.id, file=outfile)
        for block in textwrap.wrap(result.normalized_sequence, width=71):
            print("SEQUENCE", block, file=outfile)
        print(
            " ".join(field.rjust(padding) for field in _HIT_TABLE_HEADER),
            file=outfile,
        )
        for _, row in _hit_table_rows(result):
            print(" ".join(field.rjust(padding) for field in row), file=outfile)
        print("//", file=outfile)


def write_hit_output(results: Iterable[SequenceResult], outfile: TextIO) -> None:
    """Write ANARCI's seven-column per-domain HMM hit file."""
    _write_hit_header(outfile)
    _write_hit_results(results, outfile)


def _validate_output_path(value: str, label: str) -> bool:
    parent = Path(value).parent
    if str(parent) in ("", ".") or parent.exists():
        return True
    print(f"Error: {label} path does not exist", file=sys.stderr)
    return False


def _stream_records(
    records: Iterable[tuple[str, str]],
    *,
    outfile: str | bool,
    hitfile: str | bool,
    csv: bool,
    workers: int,
    batch_size: int,
    batch_residues: int,
    allowed_chains: list[str] | None,
    allowed_species: list[str],
    min_bit_score: float,
    assign_germline: bool,
) -> None:
    with ExitStack() as outputs:
        csv_output: _CsvOutputSpool | None = None
        number_output: TextIO | None = None
        if csv:
            csv_output = outputs.enter_context(_CsvOutputSpool(str(outfile)))
        else:
            output_context = (
                _atomic_text_output(Path(str(outfile)))
                if outfile
                else nullcontext(sys.stdout)
            )
            number_output = outputs.enter_context(output_context)

        hit_output: TextIO | None = None
        if hitfile:
            hit_output = outputs.enter_context(
                _atomic_text_output(Path(str(hitfile)))
            )
            _write_hit_header(hit_output)

        batches = iter_sequence_batches(
            records,
            batch_size=batch_size,
            batch_residues=batch_residues,
        )
        for batch in batches:
            results = _number_records(
                batch,
                workers=workers,
                allowed_chains=allowed_chains,
                allowed_species=allowed_species,
                min_bit_score=min_bit_score,
                alternative_hit_count=(
                    28 if hit_output is not None and allowed_chains is None else 3
                ),
                assign_germline=assign_germline,
            )
            if csv_output is not None:
                csv_output.write(results)
            else:
                assert number_output is not None
                write_anarci_output(
                    results,
                    number_output,
                    assign_germline=assign_germline,
                )

            if hit_output is None:
                del results
                del batch
                continue
            if allowed_chains is None:
                _write_hit_results(results, hit_output)
                del results
                del batch
                continue

            # ANARCI applies --restrict after constructing its hit table.
            del results
            hit_results = _number_records(
                batch,
                workers=workers,
                allowed_chains=None,
                allowed_species=allowed_species,
                min_bit_score=min_bit_score,
                alternative_hit_count=28,
                assign_germline=False,
            )
            _write_hit_results(hit_results, hit_output)
            del hit_results
            del batch

        if csv_output is not None:
            csv_output.finish()


def main(argv: Iterable[str] | None = None) -> int:
    arguments = list(argv) if argv is not None else sys.argv[1:]
    parser = _parser()
    if not arguments:
        parser.print_help()
        return 0
    args = parser.parse_args(arguments)

    if args.inputsequence is None:
        parser.error("the following arguments are required: --sequence/-i")
    if args.ncpu < 1:
        parser.error("--ncpu must be a positive integer")
    if args.csv and not args.outfile:
        print(
            "Error: When --csv option is used an ouput file name must be given.",
            file=sys.stderr,
        )
        return 1
    if args.outfile and not _validate_output_path(args.outfile, "Output file"):
        return 1
    if args.hitfile and not _validate_output_path(args.hitfile, "Hit output file"):
        return 1

    try:
        allowed_chains = _allowed_chains(args.restrict)
        allowed_species = [args.use_species] if args.use_species else ["human", "mouse"]
        with _input_records(args.inputsequence) as records:
            _stream_records(
                records,
                outfile=args.outfile,
                hitfile=args.hitfile,
                csv=args.csv,
                workers=args.ncpu,
                batch_size=DEFAULT_BATCH_SIZE,
                batch_residues=DEFAULT_BATCH_RESIDUES,
                allowed_chains=allowed_chains,
                allowed_species=allowed_species,
                min_bit_score=args.bit_score_threshold,
                assign_germline=args.assign_germline,
            )
    except (AnarcismError, OSError, ValueError) as error:
        print("Error: ", error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover - exercised through __main__
    raise SystemExit(main())
