use crate::hmm::{RawDomain, TraceState, TraceStep};
use crate::{ChainType, NumberedResidue, Region};

#[derive(Clone, Debug)]
struct LabelledResidue {
    position: u16,
    insertion_code: String,
    sequence_index: usize,
}

pub(crate) fn number_imgt(
    sequence: &[u8],
    domain: &RawDomain,
    chain_type: ChainType,
    is_first_domain: bool,
    is_only_domain: bool,
) -> (Vec<NumberedResidue>, String, usize, usize) {
    let steps = extend_terminal_steps(
        sequence,
        domain,
        chain_type,
        is_first_domain,
        is_only_domain,
    );
    let mut segments: [Vec<LabelledResidue>; 7] = std::array::from_fn(|_| Vec::new());
    let mut cdr_indices: [Vec<usize>; 3] = std::array::from_fn(|_| Vec::new());
    let mut insertion_counts = [0_usize; 129];

    for step in &steps {
        let model_position = usize::from(step.model_position);
        match (step.state, step.sequence_index) {
            (TraceState::Delete, _) | (_, None) => {}
            (TraceState::Match, Some(sequence_index)) => {
                if let Some(cdr) = cdr_index(model_position) {
                    cdr_indices[cdr].push(sequence_index);
                } else {
                    let segment = framework_segment(model_position);
                    segments[segment].push(LabelledResidue {
                        position: step.model_position,
                        insertion_code: String::new(),
                        sequence_index,
                    });
                }
            }
            (TraceState::Insert, Some(sequence_index)) => {
                if let Some(cdr) = insertion_cdr_index(model_position) {
                    cdr_indices[cdr].push(sequence_index);
                } else {
                    let position = model_position.clamp(1, 128);
                    let insertion_index = insertion_counts[position];
                    insertion_counts[position] += 1;
                    segments[framework_segment(position)].push(LabelledResidue {
                        position: position as u16,
                        insertion_code: insertion_code(insertion_index),
                        sequence_index,
                    });
                }
            }
        }
    }

    segments[1] = label_cdr(&cdr_indices[0], 12, 27, 39);
    segments[3] = label_cdr(&cdr_indices[1], 10, 56, 66);
    segments[5] = label_cdr(&cdr_indices[2], 13, 105, 118);

    let labelled: Vec<LabelledResidue> = segments.into_iter().flatten().collect();
    let start = labelled
        .iter()
        .map(|residue| residue.sequence_index)
        .min()
        .unwrap_or(domain.start);
    let end = labelled
        .iter()
        .map(|residue| residue.sequence_index)
        .max()
        .map_or(start, |index| index + 1);

    let numbering = labelled
        .iter()
        .map(|residue| NumberedResidue {
            sequence_index: residue.sequence_index,
            amino_acid: char::from(sequence[residue.sequence_index]),
            position: residue.position,
            insertion_code: residue.insertion_code.clone(),
            region: Region::for_imgt_position(residue.position),
        })
        .collect();
    let padded_alignment = padded_alignment(sequence, &labelled);
    (numbering, padded_alignment, start, end)
}

fn extend_terminal_steps(
    sequence: &[u8],
    domain: &RawDomain,
    chain_type: ChainType,
    is_first_domain: bool,
    is_only_domain: bool,
) -> Vec<TraceStep> {
    let mut steps = domain.steps.clone();
    let Some(first_emitting) = steps.iter().position(|step| step.sequence_index.is_some()) else {
        return steps;
    };
    let Some(last_emitting) = steps.iter().rposition(|step| step.sequence_index.is_some()) else {
        return steps;
    };

    if is_first_domain {
        let first = steps[first_emitting];
        let missing_model_positions = usize::from(first.model_position.saturating_sub(1));
        if missing_model_positions > 0 && missing_model_positions < 5 {
            let first_sequence_index = first.sequence_index.unwrap_or(0);
            let extension = missing_model_positions.min(first_sequence_index);
            let mut prefix = Vec::with_capacity(extension);
            for offset in (1..=extension).rev() {
                prefix.push(TraceStep {
                    state: TraceState::Match,
                    model_position: first.model_position - offset as u16,
                    sequence_index: Some(first_sequence_index - offset),
                });
            }
            steps.splice(first_emitting..first_emitting, prefix);
        }
    }

    if is_only_domain {
        let last_emitting = steps
            .iter()
            .rposition(|step| step.sequence_index.is_some())
            .unwrap_or(last_emitting);
        let last = steps[last_emitting];
        let last_sequence_index = last.sequence_index.unwrap_or(sequence.len());
        // HMMER's optimal-accuracy display retains a short, J-less CDR3 tail
        // after the conserved position-104 cysteine while leaving its last two
        // unsupported query residues outside the numbered domain. Reproduce
        // that behavior for isolated C-terminal truncations.
        if last.model_position == 104 && sequence.len() > last_sequence_index + 3 {
            let extension = (sequence.len() - last_sequence_index - 3).min(13);
            let suffix = (1..=extension).map(|offset| TraceStep {
                state: TraceState::Match,
                model_position: 104 + offset as u16,
                sequence_index: Some(last_sequence_index + offset),
            });
            steps.splice(last_emitting + 1..last_emitting + 1, suffix);
        }

        let last_emitting = steps
            .iter()
            .rposition(|step| step.sequence_index.is_some())
            .unwrap_or(last_emitting);
        let last = steps[last_emitting];
        let last_sequence_index = last.sequence_index.unwrap_or(sequence.len());
        let effective_model_end = match chain_type {
            ChainType::K | ChainType::L | ChainType::B => 127,
            _ => 128,
        };
        if last.model_position > 123
            && last.model_position < effective_model_end
            && last_sequence_index + 1 < sequence.len()
        {
            let extension = usize::from(effective_model_end - last.model_position)
                .min(sequence.len() - last_sequence_index - 1);
            let suffix = (1..=extension).map(|offset| TraceStep {
                state: TraceState::Match,
                model_position: last.model_position + offset as u16,
                sequence_index: Some(last_sequence_index + offset),
            });
            steps.splice(last_emitting + 1..last_emitting + 1, suffix);
        }
    }
    steps
}

fn cdr_index(position: usize) -> Option<usize> {
    match position {
        27..=38 => Some(0),
        56..=65 => Some(1),
        105..=117 => Some(2),
        _ => None,
    }
}

fn insertion_cdr_index(position: usize) -> Option<usize> {
    match position {
        26..=38 => Some(0),
        55..=65 => Some(1),
        104..=117 => Some(2),
        _ => None,
    }
}

fn framework_segment(position: usize) -> usize {
    match position {
        1..=26 => 0,
        27..=38 => 1,
        39..=55 => 2,
        56..=65 => 3,
        66..=104 => 4,
        105..=117 => 5,
        _ => 6,
    }
}

fn label_cdr(
    indices: &[usize],
    maximum_length: usize,
    start: u16,
    end: u16,
) -> Vec<LabelledResidue> {
    let annotations = imgt_cdr_annotations(indices.len(), maximum_length, start, end);
    annotations
        .into_iter()
        .flatten()
        .zip(indices.iter().copied())
        .map(
            |((position, insertion_code), sequence_index)| LabelledResidue {
                position,
                insertion_code,
                sequence_index,
            },
        )
        .collect()
}

fn imgt_cdr_annotations(
    length: usize,
    maximum_length: usize,
    start: u16,
    end: u16,
) -> Vec<Option<(u16, String)>> {
    let mut annotations = vec![None; length.max(maximum_length)];
    if length == 0 {
        return annotations;
    }
    if length == 1 {
        annotations[0] = Some((start, String::new()));
        return annotations;
    }

    let mut front = 0;
    let mut back = annotations.len() - 1;
    for index in 0..length.min(maximum_length) {
        if index % 2 == 0 {
            annotations[front] = Some((start + front as u16, String::new()));
            front += 1;
        } else {
            annotations[back] = Some((
                end - 1 - (annotations.len() - 1 - back) as u16,
                String::new(),
            ));
            back = back.saturating_sub(1);
        }
    }

    if length <= maximum_length {
        return annotations;
    }

    let centre_left = annotations[..front]
        .iter()
        .rev()
        .flatten()
        .next()
        .map(|annotation| annotation.0)
        .unwrap_or(start);
    let centre_right = annotations[back + 1..]
        .iter()
        .flatten()
        .next()
        .map(|annotation| annotation.0)
        .unwrap_or(end - 1);
    let mut left_insertion = 0;
    let mut right_insertion = 0;

    for insertion in 0..(length - maximum_length) {
        if insertion % 2 == 0 {
            annotations[back] = Some((centre_right, insertion_code(right_insertion)));
            right_insertion += 1;
            back = back.saturating_sub(1);
        } else {
            annotations[front] = Some((centre_left, insertion_code(left_insertion)));
            left_insertion += 1;
            front += 1;
        }
    }
    annotations
}

fn padded_alignment(sequence: &[u8], labelled: &[LabelledResidue]) -> String {
    let mut alignment = String::with_capacity(labelled.len().max(128));
    let mut previous_position = 0_u16;
    for residue in labelled {
        if residue.position > previous_position + 1 {
            alignment.extend(std::iter::repeat_n(
                '-',
                usize::from(residue.position - previous_position - 1),
            ));
        }
        alignment.push(char::from(sequence[residue.sequence_index]));
        previous_position = previous_position.max(residue.position);
    }
    alignment
}

fn insertion_code(mut index: usize) -> String {
    let mut result = String::new();
    loop {
        let digit = (index % 26) as u8;
        result.insert(0, char::from(b'A' + digit));
        if index < 26 {
            return result;
        }
        index = index / 26 - 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn imgt_cdr3_insertions_are_symmetric() {
        let annotations: Vec<String> = imgt_cdr_annotations(16, 13, 105, 118)
            .into_iter()
            .flatten()
            .map(|(position, insertion)| format!("{position}{insertion}"))
            .collect();
        assert_eq!(
            annotations,
            [
                "105", "106", "107", "108", "109", "110", "111", "111A", "112B", "112A", "112",
                "113", "114", "115", "116", "117"
            ]
        );
    }

    #[test]
    fn insertion_codes_continue_after_z() {
        assert_eq!(insertion_code(0), "A");
        assert_eq!(insertion_code(25), "Z");
        assert_eq!(insertion_code(26), "AA");
    }
}
