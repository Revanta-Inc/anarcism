#![forbid(unsafe_code)]

mod engine;
mod error;
mod fasta;
mod germlines;
mod hmm;
mod models;
mod numbering;
mod sequence;
mod types;

#[cfg(not(target_family = "wasm"))]
pub use engine::number_sequences_parallel;
pub use engine::{
    number_fasta, number_sequence, number_sequence_with_id, number_sequences,
    validate_antibody_pair,
};
pub use error::{Error, ErrorCode, Result};
pub use germlines::{GermlineDatabase, embedded_germlines};
pub use models::{Profile, ProfileDatabase, embedded_profiles};
pub use types::{
    ChainType, DomainResult, GermlineAssignment, NumberedResidue, NumberingOptions,
    PairValidationOptions, PairValidationResult, ProfileHit, ReceptorType, Region, SequenceInput,
    SequenceResult,
};
