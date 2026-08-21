"""Type stubs for the compiled extension module.

Keep this in step with `crates/anarcism-python/src/lib.rs` and with the browser
declarations in `browser/src/index.d.ts`; the three describe one contract.
"""

from typing import Any, Final, Literal

ChainType = Literal["H", "K", "L", "A", "B", "G", "D"]
ReceptorType = Literal["IG", "TR"]
Region = Literal["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
ErrorCode = Literal[
    "INVALID_SEQUENCE",
    "SEQUENCE_TOO_LONG",
    "BATCH_TOO_LARGE",
    "FASTA_TOO_LARGE",
    "INVALID_FASTA",
    "INVALID_OPTIONS",
    "CORRUPT_MODEL_DATA",
    "INTERNAL",
]

__version__: Final[str]

class AnarcismError(Exception):
    """Raised when an input is rejected or the embedded models are unreadable."""

    code: ErrorCode
    input_id: str | None

class NumberedResidue:
    @property
    def sequence_index(self) -> int: ...
    @property
    def amino_acid(self) -> str: ...
    @property
    def position(self) -> int: ...
    @property
    def insertion_code(self) -> str: ...
    @property
    def region(self) -> Region: ...
    def to_dict(self) -> dict[str, Any]: ...

class ProfileHit:
    @property
    def profile(self) -> str: ...
    @property
    def chain_type(self) -> ChainType: ...
    @property
    def species(self) -> str: ...
    @property
    def bit_score(self) -> float: ...
    @property
    def e_value(self) -> float | None: ...
    def to_dict(self) -> dict[str, Any]: ...

class GermlineAssignment:
    @property
    def species(self) -> str: ...
    @property
    def v_gene(self) -> str | None: ...
    @property
    def v_identity(self) -> float | None: ...
    @property
    def j_gene(self) -> str | None: ...
    @property
    def j_identity(self) -> float | None: ...
    def to_dict(self) -> dict[str, Any]: ...

class DomainResult:
    @property
    def domain_index(self) -> int: ...
    @property
    def receptor_type(self) -> ReceptorType: ...
    @property
    def chain_type(self) -> ChainType: ...
    @property
    def species(self) -> str: ...
    @property
    def start(self) -> int: ...
    @property
    def end(self) -> int: ...
    @property
    def bit_score(self) -> float: ...
    @property
    def e_value(self) -> float | None: ...
    @property
    def numbering(self) -> list[NumberedResidue]: ...
    @property
    def padded_imgt_alignment(self) -> str: ...
    @property
    def alternative_hits(self) -> list[ProfileHit]: ...
    @property
    def germline(self) -> GermlineAssignment | None: ...
    def to_dict(self) -> dict[str, Any]: ...

class SequenceResult:
    @property
    def id(self) -> str: ...
    @property
    def normalized_sequence(self) -> str: ...
    @property
    def domains(self) -> list[DomainResult]: ...
    @property
    def warnings(self) -> list[str]: ...
    def to_dict(self) -> dict[str, Any]: ...

class PairValidationResult:
    @property
    def ok(self) -> bool: ...
    @property
    def vh(self) -> DomainResult | None: ...
    @property
    def vl(self) -> DomainResult | None: ...
    @property
    def errors(self) -> list[str]: ...
    def to_dict(self) -> dict[str, Any]: ...
    def __bool__(self) -> bool: ...

def number_sequence(
    sequence: str,
    *,
    id: str = ...,
    allowed_chains: list[ChainType] | None = ...,
    allowed_species: list[str] | None = ...,
    min_bit_score: float | None = ...,
    alternative_hit_count: int | None = ...,
    assign_germline: bool = ...,
) -> SequenceResult: ...
def number_sequences(
    inputs: list[tuple[str, str]],
    *,
    workers: int | None = ...,
    allowed_chains: list[ChainType] | None = ...,
    allowed_species: list[str] | None = ...,
    min_bit_score: float | None = ...,
    alternative_hit_count: int | None = ...,
    assign_germline: bool = ...,
) -> list[SequenceResult]: ...
def number_fasta(
    fasta: str,
    *,
    allowed_chains: list[ChainType] | None = ...,
    allowed_species: list[str] | None = ...,
    min_bit_score: float | None = ...,
    alternative_hit_count: int | None = ...,
    assign_germline: bool = ...,
) -> list[SequenceResult]: ...
def validate_antibody_pair(
    vh: str,
    vl: str,
    *,
    start_max: int | None = ...,
    end_min: int | None = ...,
    allowed_chains: list[ChainType] | None = ...,
    allowed_species: list[str] | None = ...,
    min_bit_score: float | None = ...,
    alternative_hit_count: int | None = ...,
    assign_germline: bool = ...,
) -> PairValidationResult: ...
def chains() -> list[ChainType]: ...
def species() -> list[str]: ...
