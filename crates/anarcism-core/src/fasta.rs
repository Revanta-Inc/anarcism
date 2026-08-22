use crate::{Error, ErrorCode, Result, SequenceInput};

const MAX_FASTA_BYTES: usize = 10 * 1024 * 1024;

pub(crate) fn parse_fasta(input: &str) -> Result<Vec<SequenceInput>> {
    if input.len() > MAX_FASTA_BYTES {
        return Err(Error::sequence(
            ErrorCode::FastaTooLarge,
            format!(
                "FASTA input contains {} bytes; configured limit is {}",
                input.len(),
                MAX_FASTA_BYTES
            ),
            None,
        ));
    }

    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_sequence = String::new();
    for (line_index, raw_line) in input.lines().enumerate() {
        let line = raw_line.trim_end_matches('\r');
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                push_record(&mut records, id, &mut current_sequence)?;
            }
            let id = header.trim();
            if id.is_empty() {
                return Err(Error::sequence(
                    ErrorCode::InvalidFasta,
                    format!("FASTA header on line {} is empty", line_index + 1),
                    None,
                ));
            }
            current_id = Some(id.to_owned());
        } else if line.trim().is_empty() {
            continue;
        } else if current_id.is_none() {
            return Err(Error::sequence(
                ErrorCode::InvalidFasta,
                format!(
                    "FASTA sequence data appears before the first header on line {}",
                    line_index + 1
                ),
                None,
            ));
        } else {
            current_sequence.push_str(line);
        }
    }
    if let Some(id) = current_id {
        push_record(&mut records, id, &mut current_sequence)?;
    }
    if records.is_empty() {
        return Err(Error::sequence(
            ErrorCode::InvalidFasta,
            "FASTA input contains no records",
            None,
        ));
    }
    Ok(records)
}

fn push_record(records: &mut Vec<SequenceInput>, id: String, sequence: &mut String) -> Result<()> {
    if sequence.is_empty() {
        return Err(Error::sequence(
            ErrorCode::InvalidFasta,
            "FASTA record has an empty sequence",
            Some(&id),
        ));
    }
    records.push(SequenceInput {
        id,
        sequence: std::mem::take(sequence),
    });
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multiline_and_crlf_fasta() {
        let parsed = parse_fasta(">one description\r\nACD\r\nEF\r\n>two\nGH\n").unwrap();
        assert_eq!(parsed[0].id, "one description");
        assert_eq!(parsed[0].sequence, "ACDEF");
        assert_eq!(parsed[1].sequence, "GH");
    }

    #[test]
    fn rejects_unheaded_and_empty_records() {
        assert!(parse_fasta("ACD").is_err());
        assert!(parse_fasta(">empty\n").is_err());
    }

    #[test]
    fn record_count_is_not_limited_by_the_core_parser() {
        let fasta: String = (0..1_001)
            .map(|index| format!(">short-{index}\nA\n"))
            .collect();
        let parsed = parse_fasta(&fasta).unwrap();
        assert_eq!(parsed.len(), 1_001);
        assert_eq!(parsed.last().unwrap().id, "short-1000");
    }
}
