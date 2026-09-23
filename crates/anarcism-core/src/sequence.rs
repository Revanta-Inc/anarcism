use crate::{Error, ErrorCode, Result};

pub(crate) const AMINO_ALPHABET: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
pub(crate) const CANONICAL_RESIDUE_COUNT: usize = AMINO_ALPHABET.len();
pub(crate) const UNKNOWN_RESIDUE_INDEX: u8 = CANONICAL_RESIDUE_COUNT as u8;
pub(crate) const RESIDUE_CODE_COUNT: usize = CANONICAL_RESIDUE_COUNT + 1;
pub(crate) const MAX_SEQUENCE_LENGTH: usize = 10_000;
const INVALID_RESIDUE: u8 = u8::MAX;
const RESIDUE_INDICES: [u8; 256] = {
    let mut indices = [INVALID_RESIDUE; 256];
    let mut index = 0;
    while index < AMINO_ALPHABET.len() {
        indices[AMINO_ALPHABET[index] as usize] = index as u8;
        index += 1;
    }
    indices[b'X' as usize] = UNKNOWN_RESIDUE_INDEX;
    indices
};

#[derive(Debug)]
pub(crate) struct NormalizedSequence {
    pub text: String,
    pub encoded: Vec<u8>,
    pub warnings: Vec<String>,
}

pub(crate) fn normalize_sequence(id: &str, input: &str) -> Result<NormalizedSequence> {
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
        if encoded.len() == MAX_SEQUENCE_LENGTH {
            return Err(Error::sequence(
                ErrorCode::SequenceTooLong,
                format!(
                    "sequence exceeds the configured {} residue limit",
                    MAX_SEQUENCE_LENGTH
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
    let index = RESIDUE_INDICES[usize::from(residue)];
    (index != INVALID_RESIDUE).then_some(index)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_case_whitespace_and_unknown_residues() {
        let normalized = normalize_sequence("x", " acd\nxF ").unwrap();
        assert_eq!(normalized.text, "ACDXF");
        assert_eq!(normalized.encoded[3], UNKNOWN_RESIDUE_INDEX);
        assert_eq!(normalized.warnings.len(), 2);
        assert!(normalize_sequence("x", "ACB").is_err());
    }
}
