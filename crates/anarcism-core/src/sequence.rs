use crate::{Error, ErrorCode, Result, ValidationLimits};

pub(crate) const AMINO_ALPHABET: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

#[derive(Debug)]
pub(crate) struct NormalizedSequence {
    pub text: String,
    pub encoded: Vec<u8>,
    pub warnings: Vec<String>,
}

pub(crate) fn normalize_sequence(
    id: &str,
    input: &str,
    limits: ValidationLimits,
) -> Result<NormalizedSequence> {
    let mut text = String::with_capacity(input.len());
    let mut encoded = Vec::with_capacity(input.len());
    let mut removed_whitespace = false;
    let mut normalized_case = false;

    for byte in input.bytes() {
        if byte.is_ascii_whitespace() {
            removed_whitespace = true;
            continue;
        }
        let upper = byte.to_ascii_uppercase();
        normalized_case |= upper != byte;
        let residue = residue_index(upper).ok_or_else(|| {
            Error::sequence(
                ErrorCode::InvalidSequence,
                format!(
                    "sequence contains unsupported residue symbol {:?}",
                    char::from(byte)
                ),
                Some(id),
            )
        })?;
        if encoded.len() == limits.max_sequence_length {
            return Err(Error::sequence(
                ErrorCode::SequenceTooLong,
                format!(
                    "sequence exceeds the configured {} residue limit",
                    limits.max_sequence_length
                ),
                Some(id),
            ));
        }
        text.push(char::from(upper));
        encoded.push(residue);
    }

    if encoded.is_empty() {
        return Err(Error::sequence(
            ErrorCode::InvalidSequence,
            "sequence is empty",
            Some(id),
        ));
    }

    let mut warnings = Vec::new();
    if removed_whitespace {
        warnings.push("ASCII whitespace was removed during normalization".to_owned());
    }
    if normalized_case {
        warnings.push("lowercase residues were normalized to uppercase".to_owned());
    }
    Ok(NormalizedSequence {
        text,
        encoded,
        warnings,
    })
}

pub(crate) fn residue_index(residue: u8) -> Option<u8> {
    AMINO_ALPHABET
        .iter()
        .position(|candidate| *candidate == residue)
        .and_then(|index| index.try_into().ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_case_and_whitespace_without_accepting_ambiguous_residues() {
        let normalized = normalize_sequence("x", " acd\nEF ", ValidationLimits::default()).unwrap();
        assert_eq!(normalized.text, "ACDEF");
        assert_eq!(normalized.warnings.len(), 2);
        assert!(normalize_sequence("x", "ACX", ValidationLimits::default()).is_err());
    }
}
