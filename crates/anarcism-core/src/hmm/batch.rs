use std::collections::BTreeMap;

use wide::f32x8;

use super::{
    DD, DM, II, IM, LN_2, MD, MI, MM, SEQUENCE_BATCH_LANES, length_scores, null_one_score,
};
use crate::Profile;

/// Score one profile against many targets in sequence-major SIMD lanes.
///
/// Targets are grouped by exact length so every lane follows the same DP
/// control flow. A short final group is padded by repeating its last target;
/// only scores for real lanes are returned. The model is shared across lanes,
/// keeping its transitions and emissions hot while LLVM vectorizes the
/// independent lane arithmetic. Only two DP rows are retained because this
/// stage ranks profiles; selected profiles are traced with `ViterbiMatrix`.
pub(crate) fn viterbi_filter_bit_scores_batch(
    profile: &Profile<'_>,
    sequences: &[&[u8]],
) -> Vec<f32> {
    let mut scores = vec![f32::NEG_INFINITY; sequences.len()];
    let mut by_length: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (index, sequence) in sequences.iter().enumerate() {
        if !sequence.is_empty() {
            by_length.entry(sequence.len()).or_default().push(index);
        }
    }

    let mut workspace = SequenceBatchViterbiWorkspace::new(profile.consensus().len());
    for indices in by_length.values() {
        for chunk in indices.chunks(SEQUENCE_BATCH_LANES) {
            let mut lanes = [sequences[chunk[0]]; SEQUENCE_BATCH_LANES];
            for (lane, sequence_index) in chunk.iter().copied().enumerate() {
                lanes[lane] = sequences[sequence_index];
            }
            let lane_scores = workspace.score_equal_length(profile, &lanes);
            for (lane, sequence_index) in chunk.iter().copied().enumerate() {
                scores[sequence_index] =
                    (lane_scores[lane] - null_one_score(sequences[sequence_index].len())) / LN_2;
            }
        }
    }
    scores
}

struct SequenceBatchViterbiWorkspace {
    previous_match: Vec<f32x8>,
    previous_insert: Vec<f32x8>,
    previous_delete: Vec<f32x8>,
    current_match: Vec<f32x8>,
    current_insert: Vec<f32x8>,
    current_delete: Vec<f32x8>,
}

impl SequenceBatchViterbiWorkspace {
    fn new(model_length: usize) -> Self {
        let row = vec![f32x8::NEG_INFINITY; model_length + 1];
        Self {
            previous_match: row.clone(),
            previous_insert: row.clone(),
            previous_delete: row.clone(),
            current_match: row.clone(),
            current_insert: row.clone(),
            current_delete: row,
        }
    }

    fn score_equal_length(
        &mut self,
        profile: &Profile<'_>,
        sequences: &[&[u8]; SEQUENCE_BATCH_LANES],
    ) -> [f32; SEQUENCE_BATCH_LANES] {
        debug_assert!(
            sequences
                .iter()
                .all(|sequence| sequence.len() == sequences[0].len())
        );
        let model_length = profile.consensus().len();
        let sequence_length = sequences[0].len();
        let (loop_score, move_score) = length_scores(sequence_length);
        let local_entry_scores = profile.local_entry_scores();
        let match_score_rows = profile.match_score_rows();
        let transition_score_rows = profile.transition_score_rows();

        self.previous_match.fill(f32x8::NEG_INFINITY);
        self.previous_insert.fill(f32x8::NEG_INFINITY);
        self.previous_delete.fill(f32x8::NEG_INFINITY);
        let loop_score = f32x8::splat(loop_score);
        let move_score = f32x8::splat(move_score);
        let mut prefix = f32x8::ZERO;
        let mut join = f32x8::NEG_INFINITY;
        let mut suffix = f32x8::NEG_INFINITY;
        let mut begin = move_score;

        for (sequence_position, _) in sequences[0].iter().enumerate() {
            self.current_match.fill(f32x8::NEG_INFINITY);
            self.current_insert.fill(f32x8::NEG_INFINITY);
            self.current_delete.fill(f32x8::NEG_INFINITY);
            let mut end = f32x8::NEG_INFINITY;

            for model_position in 1..=model_length {
                let previous_transitions =
                    (model_position > 1).then(|| &transition_score_rows[model_position - 2]);
                let transitions = (model_position < model_length)
                    .then(|| &transition_score_rows[model_position - 1]);
                let (from_match, from_insert, from_delete) = previous_transitions.map_or(
                    (
                        f32x8::NEG_INFINITY,
                        f32x8::NEG_INFINITY,
                        f32x8::NEG_INFINITY,
                    ),
                    |scores| {
                        (
                            self.previous_match[model_position - 1] + f32x8::splat(scores[MM]),
                            self.previous_insert[model_position - 1] + f32x8::splat(scores[IM]),
                            self.previous_delete[model_position - 1] + f32x8::splat(scores[DM]),
                        )
                    },
                );
                let from_begin = begin + f32x8::splat(local_entry_scores[model_position - 1]);
                let emission = f32x8::new(std::array::from_fn(|lane| {
                    let residue = usize::from(sequences[lane][sequence_position]);
                    match_score_rows[model_position - 1][residue]
                }));
                let match_score =
                    max4_wide(from_match, from_insert, from_delete, from_begin) + emission;
                self.current_match[model_position] = match_score;
                end = end.fast_max(match_score);

                if let Some(scores) = transitions {
                    self.current_insert[model_position] = (self.previous_match[model_position]
                        + f32x8::splat(scores[MI]))
                    .fast_max(self.previous_insert[model_position] + f32x8::splat(scores[II]));
                }

                if let Some(scores) = previous_transitions {
                    self.current_delete[model_position] = (self.current_match[model_position - 1]
                        + f32x8::splat(scores[MD]))
                    .fast_max(self.current_delete[model_position - 1] + f32x8::splat(scores[DD]));
                    if model_position == model_length {
                        end = end.fast_max(self.current_delete[model_position]);
                    }
                }
            }

            let end_transition = end - f32x8::splat(LN_2);
            join = (join + loop_score).fast_max(end_transition);
            suffix = (suffix + loop_score).fast_max(end_transition);
            prefix += loop_score;
            begin = (prefix + move_score).fast_max(join + move_score);
            std::mem::swap(&mut self.previous_match, &mut self.current_match);
            std::mem::swap(&mut self.previous_insert, &mut self.current_insert);
            std::mem::swap(&mut self.previous_delete, &mut self.current_delete);
        }

        *(suffix + move_score).as_array()
    }
}

#[inline]
fn max4_wide(a: f32x8, b: f32x8, c: f32x8, d: f32x8) -> f32x8 {
    a.fast_max(b).fast_max(c).fast_max(d)
}
