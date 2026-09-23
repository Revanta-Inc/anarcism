"""Bounded-memory FASTA parsing and numbering helpers."""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from ._anarcism import SequenceResult, number_sequences

DEFAULT_BATCH_SIZE = 1_024
DEFAULT_BATCH_RESIDUES = 1_000_000


def iter_fasta(lines: Iterable[str]) -> Iterator[tuple[str, str]]:
    """Yield ``(identifier, sequence)`` records from FASTA text lines."""
    identifier: str | None = None
    sequence: list[str] = []
    record_count = 0

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\r\n")
        if line.startswith(">"):
            if identifier is not None:
                if not sequence:
                    raise ValueError(
                        f"FASTA record {identifier!r} has an empty sequence"
                    )
                yield identifier, "".join(sequence)
                record_count += 1
            identifier = line[1:].strip()
            if not identifier:
                raise ValueError(f"FASTA header on line {line_number} is empty")
            sequence = []
        elif line.strip():
            if identifier is None:
                raise ValueError(
                    "FASTA sequence data appears before the first header "
                    f"on line {line_number}"
                )
            sequence.append(line)

    if identifier is not None:
        if not sequence:
            raise ValueError(f"FASTA record {identifier!r} has an empty sequence")
        yield identifier, "".join(sequence)
        record_count += 1

    if record_count == 0:
        raise ValueError("no FASTA records found")


def iter_sequence_batches(
    inputs: Iterable[tuple[str, str]],
    *,
    batch_size: int = DEFAULT_BATCH_SIZE,
    batch_residues: int = DEFAULT_BATCH_RESIDUES,
) -> Iterator[list[tuple[str, str]]]:
    """Group records by count and target residues without splitting records."""
    if batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    if batch_residues < 1:
        raise ValueError("batch_residues must be a positive integer")

    batch: list[tuple[str, str]] = []
    residue_count = 0
    for record in inputs:
        record_residues = len(record[1])
        if batch and (
            len(batch) >= batch_size
            or residue_count + record_residues > batch_residues
        ):
            yield batch
            batch = []
            residue_count = 0
        batch.append(record)
        residue_count += record_residues

    if batch:
        yield batch


def iter_number_sequences(
    inputs: Iterable[tuple[str, str]],
    *,
    batch_size: int = DEFAULT_BATCH_SIZE,
    batch_residues: int = DEFAULT_BATCH_RESIDUES,
    workers: int | None = None,
    allowed_chains: list[str] | None = None,
    allowed_species: list[str] | None = None,
    min_bit_score: float | None = None,
    alternative_hit_count: int | None = None,
    assign_germline: bool = False,
) -> Iterator[SequenceResult]:
    """Number an arbitrary record iterable using bounded native batches."""
    batches = iter_sequence_batches(
        inputs,
        batch_size=batch_size,
        batch_residues=batch_residues,
    )
    for batch in batches:
        yield from number_sequences(
            batch,
            workers=workers,
            allowed_chains=allowed_chains,
            allowed_species=allowed_species,
            min_bit_score=min_bit_score,
            alternative_hit_count=alternative_hit_count,
            assign_germline=assign_germline,
        )


def iter_number_fasta(
    lines: Iterable[str],
    *,
    batch_size: int = DEFAULT_BATCH_SIZE,
    batch_residues: int = DEFAULT_BATCH_RESIDUES,
    workers: int | None = None,
    allowed_chains: list[str] | None = None,
    allowed_species: list[str] | None = None,
    min_bit_score: float | None = None,
    alternative_hit_count: int | None = None,
    assign_germline: bool = False,
) -> Iterator[SequenceResult]:
    """Parse and number FASTA lines without retaining the whole document."""
    return iter_number_sequences(
        iter_fasta(lines),
        batch_size=batch_size,
        batch_residues=batch_residues,
        workers=workers,
        allowed_chains=allowed_chains,
        allowed_species=allowed_species,
        min_bit_score=min_bit_score,
        alternative_hit_count=alternative_hit_count,
        assign_germline=assign_germline,
    )
