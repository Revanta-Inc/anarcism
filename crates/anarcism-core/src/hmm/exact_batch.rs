use wide::f32x8;

use super::{
    B, C, DELETE, DomainScoreComponents, E, INSERT, J, LN_2, Logsum, MATCH, N, PosteriorRegion,
    Profile, RawDomain, SPECIAL_COUNT, STATE_COUNT, length_scores, null_one_score,
    regions_from_posterior_totals, unihit_length_scores,
};

pub(crate) const DOMAIN_SCORE_LANES: usize = 8;

#[derive(Clone, Copy)]
pub(crate) struct DomainScoreRequest<'a> {
    pub(crate) seed: &'a RawDomain,
    pub(crate) sequence: &'a [u8],
    pub(crate) trace_null2_bias: Option<f32>,
}

pub(crate) struct DomainScoreResult {
    pub(crate) components: DomainScoreComponents,
    pub(crate) aligned_domain: Option<RawDomain>,
}

pub(crate) fn score_domains(
    profile: &Profile<'_>,
    requests: &[DomainScoreRequest<'_>],
) -> Vec<DomainScoreResult> {
    if requests.is_empty() {
        return Vec::new();
    }
    assert!(requests.len() <= DOMAIN_SCORE_LANES);
    let lanes: [DomainScoreRequest<'_>; DOMAIN_SCORE_LANES] =
        std::array::from_fn(|lane| requests.get(lane).copied().unwrap_or(requests[0]));
    let envelopes: [&[u8]; DOMAIN_SCORE_LANES] = std::array::from_fn(|lane| {
        &lanes[lane].sequence[lanes[lane].seed.start..lanes[lane].seed.end]
    });
    debug_assert!(!envelopes[0].is_empty());
    debug_assert!(
        envelopes
            .iter()
            .all(|envelope| envelope.len() == envelopes[0].len())
    );
    let needs_posterior = lanes[0].trace_null2_bias.is_none();
    debug_assert!(
        lanes
            .iter()
            .all(|request| request.trace_null2_bias.is_none() == needs_posterior)
    );
    let target_lengths = std::array::from_fn(|lane| lanes[lane].sequence.len());
    let forward = BatchForwardMatrix::fill(profile, &envelopes, target_lengths);
    let forward_scores = *forward.score().as_array();
    let posterior_lanes = needs_posterior.then(|| {
        let backward = BatchBackwardMatrix::fill(profile, &envelopes, target_lengths);
        BatchPosteriorMatrix::decode(&forward, backward).lanes(requests.len())
    });

    requests
        .iter()
        .enumerate()
        .map(|(lane, request)| {
            let envelope = envelopes[lane];
            let posterior = posterior_lanes.as_ref().map(|matrices| &matrices[lane]);
            let null2_bias_bits = request.trace_null2_bias.unwrap_or_else(|| {
                posterior
                    .expect("posterior scoring requires a posterior matrix")
                    .null2_bias_bits(profile, envelope)
            });
            let outside_envelope_nats = (request.sequence.len() - envelope.len()) as f32
                * (request.sequence.len() as f32 / (request.sequence.len() as f32 + 3.0)).ln();
            let null_one_nats = null_one_score(request.sequence.len());
            let bit_score = (forward_scores[lane] + outside_envelope_nats - null_one_nats) / LN_2
                - null2_bias_bits;
            let components = DomainScoreComponents {
                envelope_forward_nats: forward_scores[lane],
                outside_envelope_nats,
                null_one_nats,
                null2_bias_bits,
                bit_score,
            };
            let aligned_domain = posterior.and_then(|posterior| {
                posterior.optimal_accuracy_trace(profile, request.seed.start)
            });
            DomainScoreResult {
                components,
                aligned_domain,
            }
        })
        .collect()
}

pub(crate) fn posterior_regions(
    profile: &Profile<'_>,
    sequences: &[&[u8]],
) -> Vec<Vec<PosteriorRegion>> {
    if sequences.is_empty() {
        return Vec::new();
    }
    assert!(sequences.len() <= DOMAIN_SCORE_LANES);
    debug_assert!(!sequences[0].is_empty());
    debug_assert!(
        sequences
            .iter()
            .all(|sequence| sequence.len() == sequences[0].len())
    );
    let lanes: [&[u8]; DOMAIN_SCORE_LANES] =
        std::array::from_fn(|lane| sequences.get(lane).copied().unwrap_or(sequences[0]));
    let target_lengths = std::array::from_fn(|lane| lanes[lane].len());
    let forward = BatchForwardMatrix::fill_multihit(profile, &lanes, target_lengths);
    let backward = BatchBackwardMatrix::fill_multihit(profile, &lanes, target_lengths);
    let overall_score = forward.score();

    (0..sequences.len())
        .map(|lane| {
            let sequence_length = sequences[lane].len();
            let mut begin_total = vec![0.0_f32; sequence_length + 1];
            let mut end_total = vec![0.0_f32; sequence_length + 1];
            let mut model_occupancy = vec![0.0_f32; sequence_length + 1];
            for sequence_position in 1..=sequence_length {
                let begin_log_probability = special(&forward.specials, sequence_position - 1, B)
                    .as_array()[lane]
                    + special(&backward.specials, sequence_position - 1, B).as_array()[lane]
                    - overall_score.as_array()[lane];
                let end_log_probability = special(&forward.specials, sequence_position, E)
                    .as_array()[lane]
                    + special(&backward.specials, sequence_position, E).as_array()[lane]
                    - overall_score.as_array()[lane];
                begin_total[sequence_position] = begin_total[sequence_position - 1]
                    + f64::from(begin_log_probability).exp() as f32;
                end_total[sequence_position] =
                    end_total[sequence_position - 1] + f64::from(end_log_probability).exp() as f32;

                let flank_probability: f32 = [N, J, C]
                    .into_iter()
                    .map(|state| {
                        (special(&forward.specials, sequence_position - 1, state).as_array()[lane]
                            + special(&backward.specials, sequence_position, state).as_array()
                                [lane]
                            + forward.loop_score.as_array()[lane]
                            - overall_score.as_array()[lane])
                            .exp()
                    })
                    .sum();
                model_occupancy[sequence_position] = (1.0 - flank_probability).clamp(0.0, 1.0);
            }
            regions_from_posterior_totals(&begin_total, &end_total, &model_occupancy)
        })
        .collect()
}

struct BatchForwardMatrix {
    model_length: usize,
    sequence_length: usize,
    cells: Vec<f32x8>,
    specials: Vec<f32x8>,
    loop_score: f32x8,
    move_score: f32x8,
}

impl BatchForwardMatrix {
    fn fill(
        profile: &Profile<'_>,
        sequences: &[&[u8]; DOMAIN_SCORE_LANES],
        target_lengths: [usize; DOMAIN_SCORE_LANES],
    ) -> Self {
        Self::fill_with_mode(profile, sequences, target_lengths, false)
    }

    fn fill_multihit(
        profile: &Profile<'_>,
        sequences: &[&[u8]; DOMAIN_SCORE_LANES],
        target_lengths: [usize; DOMAIN_SCORE_LANES],
    ) -> Self {
        Self::fill_with_mode(profile, sequences, target_lengths, true)
    }

    fn fill_with_mode(
        profile: &Profile<'_>,
        sequences: &[&[u8]; DOMAIN_SCORE_LANES],
        target_lengths: [usize; DOMAIN_SCORE_LANES],
        multihit: bool,
    ) -> Self {
        let model_length = profile.consensus().len();
        let sequence_length = sequences[0].len();
        let logsum = WideLogsum::shared();
        let row_width = (model_length + 1) * STATE_COUNT;
        let mut cells = vec![f32x8::NEG_INFINITY; (sequence_length + 1) * row_width];
        let mut specials = vec![f32x8::NEG_INFINITY; (sequence_length + 1) * SPECIAL_COUNT];
        let length_score = if multihit {
            length_scores
        } else {
            unihit_length_scores
        };
        let loop_score = f32x8::new(std::array::from_fn(|lane| {
            length_score(target_lengths[lane]).0
        }));
        let move_score = f32x8::new(std::array::from_fn(|lane| {
            length_score(target_lengths[lane]).1
        }));
        let end_loop_score = f32x8::splat(if multihit { -LN_2 } else { f32::NEG_INFINITY });
        let end_move_score = f32x8::splat(if multihit { -LN_2 } else { 0.0 });
        set_special(&mut specials, 0, N, f32x8::ZERO);
        set_special(&mut specials, 0, B, move_score);

        let local_entry_scores = profile.local_entry_scores();
        let match_score_rows = profile.match_score_rows();
        let transition_score_rows = profile.transition_score_rows();
        let state_row_width = model_length + 1;
        let (state_cells, remainder) = cells.as_chunks_mut::<STATE_COUNT>();
        debug_assert!(remainder.is_empty());

        for sequence_position in 1..=sequence_length {
            let current_start = sequence_position * state_row_width;
            let (completed_rows, current_and_later) = state_cells.split_at_mut(current_start);
            let previous_row = &completed_rows[completed_rows.len() - state_row_width..];
            let current_row = &mut current_and_later[..state_row_width];
            let residues: [usize; DOMAIN_SCORE_LANES] =
                std::array::from_fn(|lane| usize::from(sequences[lane][sequence_position - 1]));
            let begin_state_score = special(&specials, sequence_position - 1, B);
            let mut end_score = f32x8::NEG_INFINITY;
            for model_position in 1..=model_length {
                let (from_match, from_insert, from_delete) = if model_position > 1 {
                    let transitions = &transition_score_rows[model_position - 2];
                    (
                        previous_row[model_position - 1][MATCH]
                            + f32x8::splat(transitions[super::MM]),
                        previous_row[model_position - 1][INSERT]
                            + f32x8::splat(transitions[super::IM]),
                        previous_row[model_position - 1][DELETE]
                            + f32x8::splat(transitions[super::DM]),
                    )
                } else {
                    (
                        f32x8::NEG_INFINITY,
                        f32x8::NEG_INFINITY,
                        f32x8::NEG_INFINITY,
                    )
                };
                let from_begin =
                    begin_state_score + f32x8::splat(local_entry_scores[model_position - 1]);
                let emission = f32x8::new(std::array::from_fn(|lane| {
                    match_score_rows[model_position - 1][residues[lane]]
                }));
                let match_score =
                    logsum.four(from_match, from_insert, from_begin, from_delete) + emission;
                current_row[model_position][MATCH] = match_score;

                let insert_score = if model_position < model_length {
                    let transitions = &transition_score_rows[model_position - 1];
                    logsum.pair(
                        previous_row[model_position][MATCH] + f32x8::splat(transitions[super::MI]),
                        previous_row[model_position][INSERT] + f32x8::splat(transitions[super::II]),
                    )
                } else {
                    f32x8::NEG_INFINITY
                };
                current_row[model_position][INSERT] = insert_score;

                let delete_score = if model_position > 1 {
                    let transitions = &transition_score_rows[model_position - 2];
                    logsum.pair(
                        current_row[model_position - 1][MATCH]
                            + f32x8::splat(transitions[super::MD]),
                        current_row[model_position - 1][DELETE]
                            + f32x8::splat(transitions[super::DD]),
                    )
                } else {
                    f32x8::NEG_INFINITY
                };
                current_row[model_position][DELETE] = delete_score;
                end_score = logsum.three(end_score, match_score, delete_score);
            }

            set_special(&mut specials, sequence_position, E, end_score);
            let join = logsum.pair(
                special(&specials, sequence_position - 1, J) + loop_score,
                end_score + end_loop_score,
            );
            let suffix = logsum.pair(
                special(&specials, sequence_position - 1, C) + loop_score,
                end_score + end_move_score,
            );
            set_special(&mut specials, sequence_position, C, suffix);
            set_special(&mut specials, sequence_position, J, join);
            let prefix = special(&specials, sequence_position - 1, N) + loop_score;
            set_special(&mut specials, sequence_position, N, prefix);
            set_special(
                &mut specials,
                sequence_position,
                B,
                logsum.pair(prefix + move_score, join + move_score),
            );
        }

        Self {
            model_length,
            sequence_length,
            cells,
            specials,
            loop_score,
            move_score,
        }
    }

    fn score(&self) -> f32x8 {
        special(&self.specials, self.sequence_length, C) + self.move_score
    }
}

struct BatchBackwardMatrix {
    cells: Vec<f32x8>,
    specials: Vec<f32x8>,
}

impl BatchBackwardMatrix {
    fn fill(
        profile: &Profile<'_>,
        sequences: &[&[u8]; DOMAIN_SCORE_LANES],
        target_lengths: [usize; DOMAIN_SCORE_LANES],
    ) -> Self {
        Self::fill_with_mode(profile, sequences, target_lengths, false)
    }

    fn fill_multihit(
        profile: &Profile<'_>,
        sequences: &[&[u8]; DOMAIN_SCORE_LANES],
        target_lengths: [usize; DOMAIN_SCORE_LANES],
    ) -> Self {
        Self::fill_with_mode(profile, sequences, target_lengths, true)
    }

    fn fill_with_mode(
        profile: &Profile<'_>,
        sequences: &[&[u8]; DOMAIN_SCORE_LANES],
        target_lengths: [usize; DOMAIN_SCORE_LANES],
        multihit: bool,
    ) -> Self {
        let model_length = profile.consensus().len();
        let sequence_length = sequences[0].len();
        let logsum = WideLogsum::shared();
        let row_width = (model_length + 1) * STATE_COUNT;
        let mut cells = vec![f32x8::NEG_INFINITY; (sequence_length + 1) * row_width];
        let mut specials = vec![f32x8::NEG_INFINITY; (sequence_length + 1) * SPECIAL_COUNT];
        let length_score = if multihit {
            length_scores
        } else {
            unihit_length_scores
        };
        let loop_score = f32x8::new(std::array::from_fn(|lane| {
            length_score(target_lengths[lane]).0
        }));
        let move_score = f32x8::new(std::array::from_fn(|lane| {
            length_score(target_lengths[lane]).1
        }));
        let end_loop_score = f32x8::splat(if multihit { -LN_2 } else { f32::NEG_INFINITY });
        let end_move_score = f32x8::splat(if multihit { -LN_2 } else { 0.0 });

        set_special(&mut specials, sequence_length, C, move_score);
        let terminal_end = move_score + end_move_score;
        set_special(&mut specials, sequence_length, E, terminal_end);
        let local_entry_scores = profile.local_entry_scores();
        let match_score_rows = profile.match_score_rows();
        let transition_score_rows = profile.transition_score_rows();
        let state_row_width = model_length + 1;
        let (state_cells, remainder) = cells.as_chunks_mut::<STATE_COUNT>();
        debug_assert!(remainder.is_empty());

        let terminal_start = sequence_length * state_row_width;
        let terminal_row = &mut state_cells[terminal_start..terminal_start + state_row_width];
        terminal_row[model_length][MATCH] = terminal_end;
        terminal_row[model_length][DELETE] = terminal_end;
        for model_position in (1..model_length).rev() {
            let transitions = &transition_score_rows[model_position - 1];
            let next_delete = terminal_row[model_position + 1][DELETE];
            terminal_row[model_position][MATCH] = logsum.pair(
                terminal_end,
                next_delete + f32x8::splat(transitions[super::MD]),
            );
            terminal_row[model_position][DELETE] = logsum.pair(
                terminal_end,
                next_delete + f32x8::splat(transitions[super::DD]),
            );
        }

        for sequence_position in (1..sequence_length).rev() {
            let next_start = (sequence_position + 1) * state_row_width;
            let (through_current, next_and_later) = state_cells.split_at_mut(next_start);
            let current_start = sequence_position * state_row_width;
            let current_row = &mut through_current[current_start..current_start + state_row_width];
            let next_row = &next_and_later[..state_row_width];
            let residues: [usize; DOMAIN_SCORE_LANES] =
                std::array::from_fn(|lane| usize::from(sequences[lane][sequence_position]));
            let mut begin_score = f32x8::NEG_INFINITY;
            for model_position in 1..=model_length {
                let emission = f32x8::new(std::array::from_fn(|lane| {
                    match_score_rows[model_position - 1][residues[lane]]
                }));
                begin_score = logsum.pair(
                    begin_score,
                    next_row[model_position][MATCH]
                        + f32x8::splat(local_entry_scores[model_position - 1])
                        + emission,
                );
            }
            set_special(&mut specials, sequence_position, B, begin_score);
            let join = logsum.pair(
                special(&specials, sequence_position + 1, J) + loop_score,
                begin_score + move_score,
            );
            set_special(&mut specials, sequence_position, J, join);
            let suffix = special(&specials, sequence_position + 1, C) + loop_score;
            set_special(&mut specials, sequence_position, C, suffix);
            let end_score = logsum.pair(join + end_loop_score, suffix + end_move_score);
            set_special(&mut specials, sequence_position, E, end_score);
            let prefix = logsum.pair(
                special(&specials, sequence_position + 1, N) + loop_score,
                begin_score + move_score,
            );
            set_special(&mut specials, sequence_position, N, prefix);

            current_row[model_length][MATCH] = end_score;
            current_row[model_length][DELETE] = end_score;
            for model_position in (1..model_length).rev() {
                let transitions = &transition_score_rows[model_position - 1];
                let emission = f32x8::new(std::array::from_fn(|lane| {
                    match_score_rows[model_position][residues[lane]]
                }));
                let next_match = next_row[model_position + 1][MATCH] + emission;
                current_row[model_position][MATCH] = logsum.four(
                    next_match + f32x8::splat(transitions[super::MM]),
                    next_row[model_position][INSERT] + f32x8::splat(transitions[super::MI]),
                    end_score,
                    current_row[model_position + 1][DELETE] + f32x8::splat(transitions[super::MD]),
                );
                current_row[model_position][INSERT] = logsum.pair(
                    next_match + f32x8::splat(transitions[super::IM]),
                    next_row[model_position][INSERT] + f32x8::splat(transitions[super::II]),
                );
                current_row[model_position][DELETE] = logsum.three(
                    next_match + f32x8::splat(transitions[super::DM]),
                    current_row[model_position + 1][DELETE] + f32x8::splat(transitions[super::DD]),
                    end_score,
                );
            }
        }

        if sequence_length > 0 {
            let first_row = &state_cells[state_row_width..state_row_width * 2];
            let residues: [usize; DOMAIN_SCORE_LANES] =
                std::array::from_fn(|lane| usize::from(sequences[lane][0]));
            let mut begin_score = f32x8::NEG_INFINITY;
            for model_position in 1..=model_length {
                let emission = f32x8::new(std::array::from_fn(|lane| {
                    match_score_rows[model_position - 1][residues[lane]]
                }));
                begin_score = logsum.pair(
                    begin_score,
                    first_row[model_position][MATCH]
                        + f32x8::splat(local_entry_scores[model_position - 1])
                        + emission,
                );
            }
            set_special(&mut specials, 0, B, begin_score);
            let prefix = logsum.pair(
                special(&specials, 1, N) + loop_score,
                begin_score + move_score,
            );
            set_special(&mut specials, 0, N, prefix);
        }

        Self { cells, specials }
    }
}

struct BatchPosteriorMatrix {
    model_length: usize,
    sequence_length: usize,
    cells: Vec<f32x8>,
    specials: Vec<f32x8>,
}

impl BatchPosteriorMatrix {
    fn decode(forward: &BatchForwardMatrix, mut backward: BatchBackwardMatrix) -> Self {
        let overall_score = forward.score();
        for sequence_position in 1..=forward.sequence_length {
            for model_position in 1..=forward.model_length {
                let match_probability = exp_lanes(
                    matrix_cell(
                        &forward.cells,
                        forward.model_length,
                        sequence_position,
                        model_position,
                        MATCH,
                    ) + matrix_cell(
                        &backward.cells,
                        forward.model_length,
                        sequence_position,
                        model_position,
                        MATCH,
                    ) - overall_score,
                );
                let insert_probability = if model_position < forward.model_length {
                    exp_lanes(
                        matrix_cell(
                            &forward.cells,
                            forward.model_length,
                            sequence_position,
                            model_position,
                            INSERT,
                        ) + matrix_cell(
                            &backward.cells,
                            forward.model_length,
                            sequence_position,
                            model_position,
                            INSERT,
                        ) - overall_score,
                    )
                } else {
                    f32x8::ZERO
                };
                set_matrix_cell(
                    &mut backward.cells,
                    forward.model_length,
                    sequence_position,
                    model_position,
                    MATCH,
                    match_probability,
                );
                set_matrix_cell(
                    &mut backward.cells,
                    forward.model_length,
                    sequence_position,
                    model_position,
                    INSERT,
                    insert_probability,
                );
                set_matrix_cell(
                    &mut backward.cells,
                    forward.model_length,
                    sequence_position,
                    model_position,
                    DELETE,
                    f32x8::ZERO,
                );
            }

            for state in [N, J, C] {
                let probability = exp_lanes(
                    special(&forward.specials, sequence_position - 1, state)
                        + special(&backward.specials, sequence_position, state)
                        + forward.loop_score
                        - overall_score,
                );
                set_special(
                    &mut backward.specials,
                    sequence_position,
                    state,
                    probability,
                );
            }
            set_special(&mut backward.specials, sequence_position, E, f32x8::ZERO);
            set_special(&mut backward.specials, sequence_position, B, f32x8::ZERO);
        }
        backward.cells[..(forward.model_length + 1) * STATE_COUNT].fill(f32x8::ZERO);
        backward.specials[..SPECIAL_COUNT].fill(f32x8::ZERO);
        Self {
            model_length: forward.model_length,
            sequence_length: forward.sequence_length,
            cells: backward.cells,
            specials: backward.specials,
        }
    }

    fn lanes(&self, lane_count: usize) -> Vec<super::PosteriorMatrix> {
        let mut cells: Vec<Vec<f32>> = (0..lane_count)
            .map(|_| Vec::with_capacity(self.cells.len()))
            .collect();
        for value in &self.cells {
            for (lane, cells) in cells.iter_mut().enumerate() {
                cells.push(value.as_array()[lane]);
            }
        }
        let mut specials: Vec<Vec<f32>> = (0..lane_count)
            .map(|_| Vec::with_capacity(self.specials.len()))
            .collect();
        for value in &self.specials {
            for (lane, specials) in specials.iter_mut().enumerate() {
                specials.push(value.as_array()[lane]);
            }
        }
        cells
            .into_iter()
            .zip(specials)
            .map(|(cells, specials)| super::PosteriorMatrix {
                model_length: self.model_length,
                sequence_length: self.sequence_length,
                cells,
                specials,
            })
            .collect()
    }
}

#[derive(Clone, Copy)]
struct WideLogsum(Logsum);

impl WideLogsum {
    fn shared() -> Self {
        Self(Logsum::shared())
    }

    #[inline(always)]
    fn pair(self, a: f32x8, b: f32x8) -> f32x8 {
        let difference = (a - b).abs();
        let correction = f32x8::new(std::array::from_fn(|lane| {
            let difference = difference.as_array()[lane];
            if difference >= 15.7 {
                0.0
            } else {
                self.0.0[(difference * super::LOGSUM_SCALE) as usize]
            }
        }));
        a.fast_max(b) + correction
    }

    #[inline(always)]
    fn three(self, a: f32x8, b: f32x8, c: f32x8) -> f32x8 {
        self.pair(self.pair(a, b), c)
    }

    #[inline(always)]
    fn four(self, a: f32x8, b: f32x8, c: f32x8, d: f32x8) -> f32x8 {
        self.pair(self.pair(a, b), self.pair(c, d))
    }
}

#[inline(always)]
fn exp_lanes(value: f32x8) -> f32x8 {
    let value = value.as_array();
    f32x8::new(std::array::from_fn(|lane| value[lane].exp()))
}

fn matrix_cell(
    cells: &[f32x8],
    model_length: usize,
    sequence_position: usize,
    model_position: usize,
    state: usize,
) -> f32x8 {
    cells[((sequence_position * (model_length + 1) + model_position) * STATE_COUNT) + state]
}

fn set_matrix_cell(
    cells: &mut [f32x8],
    model_length: usize,
    sequence_position: usize,
    model_position: usize,
    state: usize,
    value: f32x8,
) {
    let index = ((sequence_position * (model_length + 1) + model_position) * STATE_COUNT) + state;
    cells[index] = value;
}

fn special(specials: &[f32x8], sequence_position: usize, state: usize) -> f32x8 {
    specials[sequence_position * SPECIAL_COUNT + state]
}

fn set_special(specials: &mut [f32x8], sequence_position: usize, state: usize, value: f32x8) {
    specials[sequence_position * SPECIAL_COUNT + state] = value;
}

#[cfg(test)]
mod tests {
    use crate::{
        embedded_profiles,
        sequence::{CANONICAL_RESIDUE_COUNT, UNKNOWN_RESIDUE_INDEX, residue_index},
    };

    use super::*;
    use crate::hmm::{
        define_domain, domain_score_components_with_posterior,
        posterior_regions as scalar_posterior_regions, viterbi_domains,
    };

    const VH: &str = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";

    fn encode(sequence: &str) -> Vec<u8> {
        sequence
            .bytes()
            .map(|residue| residue_index(residue).unwrap())
            .collect()
    }

    fn assert_components_equal(
        observed: DomainScoreComponents,
        expected: DomainScoreComponents,
        lane: usize,
    ) {
        assert_eq!(
            observed.envelope_forward_nats, expected.envelope_forward_nats,
            "Forward score in lane {lane}"
        );
        assert_eq!(
            observed.outside_envelope_nats, expected.outside_envelope_nats,
            "outside-envelope score in lane {lane}"
        );
        assert_eq!(
            observed.null_one_nats, expected.null_one_nats,
            "null1 score in lane {lane}"
        );
        assert_eq!(
            observed.null2_bias_bits, expected.null2_bias_bits,
            "null2 score in lane {lane}"
        );
        assert_eq!(
            observed.bit_score, expected.bit_score,
            "bit score in lane {lane}"
        );
    }

    fn assert_domain_equal(
        observed: Option<&RawDomain>,
        expected: Option<&RawDomain>,
        lane: usize,
    ) {
        match (observed, expected) {
            (Some(observed), Some(expected)) => {
                assert_eq!(observed.start, expected.start, "start in lane {lane}");
                assert_eq!(observed.end, expected.end, "end in lane {lane}");
                assert_eq!(observed.steps, expected.steps, "trace in lane {lane}");
            }
            (None, None) => {}
            _ => panic!("alignment presence differs in lane {lane}"),
        }
    }

    #[test]
    fn domain_batches_match_scalar_scoring_and_alignment() {
        let database = embedded_profiles().unwrap();
        let profile = &database.profiles()[0];
        let base = encode(VH);
        let mut seed = viterbi_domains(profile, &base).remove(0);
        let definition = define_domain(profile, &seed, &base);
        seed.start = definition.start;
        seed.end = definition.end;

        let mut sequences = Vec::with_capacity(DOMAIN_SCORE_LANES);
        for lane in 0..DOMAIN_SCORE_LANES {
            let mut sequence = base.clone();
            sequence[seed.start + 4 + lane] = if lane == 0 {
                UNKNOWN_RESIDUE_INDEX
            } else {
                ((usize::from(sequence[seed.start + 4 + lane]) + lane + 1)
                    % CANONICAL_RESIDUE_COUNT) as u8
            };
            sequence.extend(std::iter::repeat_n(0, lane));
            sequences.push(sequence);
        }
        let seeds = vec![seed; DOMAIN_SCORE_LANES];

        for count in [1, 3, DOMAIN_SCORE_LANES] {
            let requests: Vec<_> = (0..count)
                .map(|lane| DomainScoreRequest {
                    seed: &seeds[lane],
                    sequence: &sequences[lane],
                    trace_null2_bias: None,
                })
                .collect();
            let observed = score_domains(profile, &requests);

            for lane in 0..count {
                let (expected_components, expected_alignment) =
                    domain_score_components_with_posterior(
                        profile,
                        &seeds[lane],
                        &sequences[lane],
                        None,
                    );
                let expected_domain =
                    expected_alignment.and_then(|alignment| alignment.align(profile));
                assert_components_equal(observed[lane].components, expected_components, lane);
                assert_domain_equal(
                    observed[lane].aligned_domain.as_ref(),
                    expected_domain.as_ref(),
                    lane,
                );
            }
        }
    }

    #[test]
    fn domain_batches_match_scalar_trace_bias_scoring() {
        let database = embedded_profiles().unwrap();
        let profile = &database.profiles()[0];
        let sequence = encode(VH);
        let mut seed = viterbi_domains(profile, &sequence).remove(0);
        let definition = define_domain(profile, &seed, &sequence);
        seed.start = definition.start;
        seed.end = definition.end;
        let bias = 4.25;
        let requests = vec![
            DomainScoreRequest {
                seed: &seed,
                sequence: &sequence,
                trace_null2_bias: Some(bias),
            };
            DOMAIN_SCORE_LANES
        ];

        let observed = score_domains(profile, &requests);
        let (expected, alignment) =
            domain_score_components_with_posterior(profile, &seed, &sequence, Some(bias));
        assert!(alignment.is_none());
        for (lane, result) in observed.iter().enumerate() {
            assert_components_equal(result.components, expected, lane);
            assert!(result.aligned_domain.is_none());
        }
    }

    #[test]
    fn posterior_region_batches_match_scalar_domain_definition() {
        let database = embedded_profiles().unwrap();
        let profile = &database.profiles()[0];
        let base = encode(VH);
        let sequences: Vec<_> = (0..DOMAIN_SCORE_LANES)
            .map(|lane| {
                let mut sequence = base.clone();
                sequence[10 + lane] =
                    ((usize::from(sequence[10 + lane]) + lane + 1) % CANONICAL_RESIDUE_COUNT) as u8;
                sequence
            })
            .collect();
        let references: Vec<_> = sequences.iter().map(Vec::as_slice).collect();

        for count in [1, 3, DOMAIN_SCORE_LANES] {
            let observed = posterior_regions(profile, &references[..count]);
            for lane in 0..count {
                let expected = scalar_posterior_regions(profile, &sequences[lane]);
                assert_eq!(observed[lane].len(), expected.len(), "lane {lane}");
                for (observed, expected) in observed[lane].iter().zip(expected) {
                    assert_eq!(observed.start, expected.start, "start in lane {lane}");
                    assert_eq!(observed.end, expected.end, "end in lane {lane}");
                    assert_eq!(
                        observed.is_multidomain, expected.is_multidomain,
                        "multidomain flag in lane {lane}"
                    );
                }
            }
        }
    }
}
