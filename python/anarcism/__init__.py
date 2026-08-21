"""ANARCI-compatible IMGT numbering for antibody and TCR variable domains.

The analysis surface is implemented in Rust and exposed by the compiled
``anarcism._anarcism`` extension module, which embeds every profile and germline
model. Nothing here reads a data file, starts a subprocess, or sends a sequence
anywhere.
"""

from ._anarcism import (
    AnarcismError,
    __version__,
    chains,
    number_fasta,
    number_sequence,
    number_sequences,
    number_sequences_parallel,
    species,
    validate_antibody_pair,
)

__all__ = [
    "AnarcismError",
    "__version__",
    "chains",
    "number_fasta",
    "number_sequence",
    "number_sequences",
    "number_sequences_parallel",
    "species",
    "validate_antibody_pair",
]
