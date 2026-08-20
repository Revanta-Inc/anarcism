use std::fmt::{self, Display};

use serde::{Deserialize, Serialize};

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum ErrorCode {
    InvalidSequence,
    SequenceTooLong,
    BatchTooLarge,
    FastaTooLarge,
    InvalidFasta,
    InvalidOptions,
    CorruptModelData,
    Internal,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Error {
    pub code: ErrorCode,
    pub message: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_id: Option<String>,
}

impl Error {
    pub(crate) fn model(message: impl Into<String>) -> Self {
        Self {
            code: ErrorCode::CorruptModelData,
            message: message.into(),
            input_id: None,
        }
    }

    pub(crate) fn sequence(
        code: ErrorCode,
        message: impl Into<String>,
        input_id: Option<&str>,
    ) -> Self {
        Self {
            code,
            message: message.into(),
            input_id: input_id.map(str::to_owned),
        }
    }
}

impl Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(input_id) = &self.input_id {
            write!(f, "{} ({input_id})", self.message)
        } else {
            f.write_str(&self.message)
        }
    }
}

impl std::error::Error for Error {}
