"""ANARCI-compatible IMGT numbering for antibody and TCR variable domains.

The analysis surface is implemented in Rust and exposed by the compiled
``anarcism._anarcism`` extension module, which embeds every profile and germline
model. Nothing here reads a data file, starts a subprocess, or sends a sequence
anywhere.

Every analysis call releases the GIL, so a :class:`concurrent.futures.ThreadPoolExecutor`
scales across cores. :func:`number_sequences` also accepts ``workers`` to spread
one batch over native threads.
"""

from ._anarcism import (
    AnarcismError,
    DomainResult,
    GermlineAssignment,
    NumberedResidue,
    PairValidationResult,
    ProfileHit,
    SequenceResult,
    __version__,
    chains,
    number_fasta,
    number_sequence,
    number_sequences,
    species,
    validate_antibody_pair,
)
from .streaming import iter_fasta, iter_number_fasta, iter_number_sequences

__all__ = [
    "AnarcismError",
    "DomainResult",
    "GermlineAssignment",
    "NumberedResidue",
    "PairValidationResult",
    "ProfileHit",
    "SequenceResult",
    "__version__",
    "chains",
    "iter_fasta",
    "iter_number_fasta",
    "iter_number_sequences",
    "number_fasta",
    "number_sequence",
    "number_sequences",
    "species",
    "validate_antibody_pair",
]
