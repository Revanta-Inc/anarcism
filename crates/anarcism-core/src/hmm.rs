use crate::Profile;

const MATCH: usize = 0;
const INSERT: usize = 1;
const DELETE: usize = 2;
const STATE_COUNT: usize = 3;

const MM: usize = 0;
const MI: usize = 1;
const MD: usize = 2;
const IM: usize = 3;
const II: usize = 4;
const DM: usize = 5;
const DD: usize = 6;

const E: usize = 0;
const N: usize = 1;
const J: usize = 2;
const B: usize = 3;
const C: usize = 4;
const SPECIAL_COUNT: usize = 5;

const LN_2: f32 = std::f32::consts::LN_2;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TraceState {
    Match,
    Insert,
    Delete,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct TraceStep {
    pub state: TraceState,
    pub model_position: u16,
    pub sequence_index: Option<usize>,
}

#[derive(Clone, Debug)]
pub(crate) struct RawDomain {
    pub steps: Vec<TraceStep>,
    pub start: usize,
    pub end: usize,
}

pub(crate) fn viterbi_domains(profile: &Profile<'_>, sequence: &[u8]) -> Vec<RawDomain> {
    if sequence.is_empty() {
        return Vec::new();
    }
    let matrix = ViterbiMatrix::fill(profile, sequence);
    matrix.trace(profile, sequence)
}

/// Cheap profile ordering score from a Viterbi trace. The exact public bit
/// score is still computed with Forward/null1/null2 for the retained profiles.
pub(crate) fn trace_emission_score(
    profile: &Profile<'_>,
    domain: &RawDomain,
    sequence: &[u8],
) -> f32 {
    domain
        .steps
        .iter()
        .filter(|step| step.state == TraceState::Match)
        .filter_map(|step| {
            step.sequence_index
                .map(|index| (step.model_position, index))
        })
        .filter(|(_, index)| *index < sequence.len())
        .map(|(position, index)| {
            profile.match_score(usize::from(position), usize::from(sequence[index]))
        })
        .sum()
}

/// Ungapped local match score used only to decide which profiles merit the
/// substantially more expensive gapped Viterbi/Forward pass.
pub(crate) fn ungapped_filter_score(profile: &Profile<'_>, sequence: &[u8]) -> f32 {
    let model_length = profile.consensus().len();
    let mut previous = vec![0.0_f32; model_length + 1];
    let mut current = vec![0.0_f32; model_length + 1];
    let mut best = 0.0_f32;
    for residue in sequence {
        for model_position in 1..=model_length {
            let score = previous[model_position - 1]
                + profile.match_score(model_position, usize::from(*residue));
            current[model_position] = score.max(0.0);
            best = best.max(current[model_position]);
        }
        std::mem::swap(&mut previous, &mut current);
    }
    best
}

/// HMMER3 generic Forward score without the optional null2 composition-bias
/// correction. The profile and null length models match `hmmscan` local,
/// multihit configuration.
pub(crate) fn forward_bit_score(profile: &Profile<'_>, sequence: &[u8]) -> f32 {
    if sequence.is_empty() {
        return f32::NEG_INFINITY;
    }
    let model_length = profile.consensus().len();
    let sequence_length = sequence.len();
    let (loop_score, move_score) = length_scores(sequence_length);
    let end_loop = -LN_2;
    let end_move = -LN_2;

    let mut previous = vec![f32::NEG_INFINITY; (model_length + 1) * STATE_COUNT];
    let mut current = vec![f32::NEG_INFINITY; (model_length + 1) * STATE_COUNT];
    let mut previous_special = [f32::NEG_INFINITY; SPECIAL_COUNT];
    previous_special[N] = 0.0;
    previous_special[B] = move_score;

    for residue in sequence {
        current.fill(f32::NEG_INFINITY);
        let mut special = [f32::NEG_INFINITY; SPECIAL_COUNT];

        for model_position in 1..model_length {
            let emission = profile.match_score(model_position, usize::from(*residue));
            let from_match = if model_position > 1 {
                cell(&previous, model_length, model_position - 1, MATCH)
                    + profile.transition_score(model_position - 1, MM)
            } else {
                f32::NEG_INFINITY
            };
            let from_insert = if model_position > 1 {
                cell(&previous, model_length, model_position - 1, INSERT)
                    + profile.transition_score(model_position - 1, IM)
            } else {
                f32::NEG_INFINITY
            };
            let from_delete = if model_position > 1 {
                cell(&previous, model_length, model_position - 1, DELETE)
                    + profile.transition_score(model_position - 1, DM)
            } else {
                f32::NEG_INFINITY
            };
            let from_begin = previous_special[B] + profile.local_entry_score(model_position);
            let match_score = logsum4(from_match, from_insert, from_begin, from_delete) + emission;
            set_cell(
                &mut current,
                model_length,
                model_position,
                MATCH,
                match_score,
            );

            let insert_score = logsum(
                cell(&previous, model_length, model_position, MATCH)
                    + profile.transition_score(model_position, MI),
                cell(&previous, model_length, model_position, INSERT)
                    + profile.transition_score(model_position, II),
            );
            set_cell(
                &mut current,
                model_length,
                model_position,
                INSERT,
                insert_score,
            );

            let delete_score = if model_position > 1 {
                logsum(
                    cell(&current, model_length, model_position - 1, MATCH)
                        + profile.transition_score(model_position - 1, MD),
                    cell(&current, model_length, model_position - 1, DELETE)
                        + profile.transition_score(model_position - 1, DD),
                )
            } else {
                f32::NEG_INFINITY
            };
            set_cell(
                &mut current,
                model_length,
                model_position,
                DELETE,
                delete_score,
            );
            special[E] = logsum3(special[E], match_score, delete_score);
        }

        let model_position = model_length;
        let emission = profile.match_score(model_position, usize::from(*residue));
        let match_score = logsum4(
            cell(&previous, model_length, model_position - 1, MATCH)
                + profile.transition_score(model_position - 1, MM),
            cell(&previous, model_length, model_position - 1, INSERT)
                + profile.transition_score(model_position - 1, IM),
            previous_special[B] + profile.local_entry_score(model_position),
            cell(&previous, model_length, model_position - 1, DELETE)
                + profile.transition_score(model_position - 1, DM),
        ) + emission;
        set_cell(
            &mut current,
            model_length,
            model_position,
            MATCH,
            match_score,
        );
        let delete_score = logsum(
            cell(&current, model_length, model_position - 1, MATCH)
                + profile.transition_score(model_position - 1, MD),
            cell(&current, model_length, model_position - 1, DELETE)
                + profile.transition_score(model_position - 1, DD),
        );
        set_cell(
            &mut current,
            model_length,
            model_position,
            DELETE,
            delete_score,
        );
        special[E] = logsum3(special[E], match_score, delete_score);

        special[J] = logsum(previous_special[J] + loop_score, special[E] + end_loop);
        special[C] = logsum(previous_special[C] + loop_score, special[E] + end_move);
        special[N] = previous_special[N] + loop_score;
        special[B] = logsum(special[N] + move_score, special[J] + move_score);

        std::mem::swap(&mut previous, &mut current);
        previous_special = special;
    }

    let raw_score = previous_special[C] + move_score;
    let null_score = null_one_score(sequence_length);
    (raw_score - null_score) / LN_2
}

/// Deterministic null2 composition correction estimated from the Viterbi
/// domain trace. HMMER uses posterior expected state occupancy for simple
/// envelopes; this trace form is retained as a compact fallback and is close
/// enough to make the remaining parity error measurable during the gate.
pub(crate) fn trace_null2_bias_bits(
    profile: &Profile<'_>,
    domain: &RawDomain,
    sequence: &[u8],
) -> f32 {
    let model_length = profile.consensus().len();
    let mut match_usage = vec![0.0_f32; model_length + 1];
    let mut insert_usage = 0.0_f32;
    let mut emitted = 0.0_f32;
    for step in &domain.steps {
        if step.sequence_index.is_none() {
            continue;
        }
        emitted += 1.0;
        match step.state {
            TraceState::Match => match_usage[usize::from(step.model_position)] += 1.0,
            TraceState::Insert => insert_usage += 1.0,
            TraceState::Delete => {}
        }
    }
    if emitted == 0.0 {
        return 0.0;
    }

    let mut odds = [0.0_f32; 20];
    for (residue, value) in odds.iter_mut().enumerate() {
        let mut weighted = insert_usage;
        for (model_position, usage) in match_usage.iter().enumerate().skip(1) {
            if *usage != 0.0 {
                weighted += usage * profile.match_score(model_position, residue).exp();
            }
        }
        *value = weighted / emitted;
    }
    let correction: f32 = sequence[domain.start..domain.end]
        .iter()
        .map(|residue| odds[usize::from(*residue)].ln())
        .sum();
    (1.0 + (1.0 / 256.0) * correction.exp()).ln() / LN_2
}

fn null_one_score(length: usize) -> f32 {
    let p1 = length as f32 / (length as f32 + 1.0);
    length as f32 * p1.ln() + (1.0 - p1).ln()
}

fn length_scores(length: usize) -> (f32, f32) {
    let move_probability = 3.0 / (length as f32 + 3.0);
    ((1.0 - move_probability).ln(), move_probability.ln())
}

struct ViterbiMatrix {
    model_length: usize,
    sequence_length: usize,
    cells: Vec<f32>,
    specials: Vec<f32>,
    loop_score: f32,
    move_score: f32,
}

impl ViterbiMatrix {
    fn fill(profile: &Profile<'_>, sequence: &[u8]) -> Self {
        let model_length = profile.consensus().len();
        let sequence_length = sequence.len();
        let row_width = (model_length + 1) * STATE_COUNT;
        let mut cells = vec![f32::NEG_INFINITY; (sequence_length + 1) * row_width];
        let mut specials = vec![f32::NEG_INFINITY; (sequence_length + 1) * SPECIAL_COUNT];
        let (loop_score, move_score) = length_scores(sequence_length);
        set_special(&mut specials, 0, N, 0.0);
        set_special(&mut specials, 0, B, move_score);

        for sequence_position in 1..=sequence_length {
            let residue = usize::from(sequence[sequence_position - 1]);
            let mut end_score = f32::NEG_INFINITY;

            for model_position in 1..model_length {
                let from_match = if model_position > 1 {
                    matrix_cell(
                        &cells,
                        model_length,
                        sequence_position - 1,
                        model_position - 1,
                        MATCH,
                    ) + profile.transition_score(model_position - 1, MM)
                } else {
                    f32::NEG_INFINITY
                };
                let from_insert = if model_position > 1 {
                    matrix_cell(
                        &cells,
                        model_length,
                        sequence_position - 1,
                        model_position - 1,
                        INSERT,
                    ) + profile.transition_score(model_position - 1, IM)
                } else {
                    f32::NEG_INFINITY
                };
                let from_delete = if model_position > 1 {
                    matrix_cell(
                        &cells,
                        model_length,
                        sequence_position - 1,
                        model_position - 1,
                        DELETE,
                    ) + profile.transition_score(model_position - 1, DM)
                } else {
                    f32::NEG_INFINITY
                };
                let from_begin = special(&specials, sequence_position - 1, B)
                    + profile.local_entry_score(model_position);
                let match_score = max4(from_match, from_insert, from_delete, from_begin)
                    + profile.match_score(model_position, residue);
                set_matrix_cell(
                    &mut cells,
                    model_length,
                    sequence_position,
                    model_position,
                    MATCH,
                    match_score,
                );
                end_score = end_score.max(match_score);

                let insert_score = (matrix_cell(
                    &cells,
                    model_length,
                    sequence_position - 1,
                    model_position,
                    MATCH,
                ) + profile.transition_score(model_position, MI))
                .max(
                    matrix_cell(
                        &cells,
                        model_length,
                        sequence_position - 1,
                        model_position,
                        INSERT,
                    ) + profile.transition_score(model_position, II),
                );
                set_matrix_cell(
                    &mut cells,
                    model_length,
                    sequence_position,
                    model_position,
                    INSERT,
                    insert_score,
                );

                let delete_score = if model_position > 1 {
                    (matrix_cell(
                        &cells,
                        model_length,
                        sequence_position,
                        model_position - 1,
                        MATCH,
                    ) + profile.transition_score(model_position - 1, MD))
                    .max(
                        matrix_cell(
                            &cells,
                            model_length,
                            sequence_position,
                            model_position - 1,
                            DELETE,
                        ) + profile.transition_score(model_position - 1, DD),
                    )
                } else {
                    f32::NEG_INFINITY
                };
                set_matrix_cell(
                    &mut cells,
                    model_length,
                    sequence_position,
                    model_position,
                    DELETE,
                    delete_score,
                );
            }

            let model_position = model_length;
            let match_score = max4(
                matrix_cell(
                    &cells,
                    model_length,
                    sequence_position - 1,
                    model_position - 1,
                    MATCH,
                ) + profile.transition_score(model_position - 1, MM),
                matrix_cell(
                    &cells,
                    model_length,
                    sequence_position - 1,
                    model_position - 1,
                    INSERT,
                ) + profile.transition_score(model_position - 1, IM),
                matrix_cell(
                    &cells,
                    model_length,
                    sequence_position - 1,
                    model_position - 1,
                    DELETE,
                ) + profile.transition_score(model_position - 1, DM),
                special(&specials, sequence_position - 1, B)
                    + profile.local_entry_score(model_position),
            ) + profile.match_score(model_position, residue);
            set_matrix_cell(
                &mut cells,
                model_length,
                sequence_position,
                model_position,
                MATCH,
                match_score,
            );
            let delete_score = (matrix_cell(
                &cells,
                model_length,
                sequence_position,
                model_position - 1,
                MATCH,
            ) + profile.transition_score(model_position - 1, MD))
            .max(
                matrix_cell(
                    &cells,
                    model_length,
                    sequence_position,
                    model_position - 1,
                    DELETE,
                ) + profile.transition_score(model_position - 1, DD),
            );
            set_matrix_cell(
                &mut cells,
                model_length,
                sequence_position,
                model_position,
                DELETE,
                delete_score,
            );
            end_score = end_score.max(match_score).max(delete_score);
            set_special(&mut specials, sequence_position, E, end_score);

            let join =
                (special(&specials, sequence_position - 1, J) + loop_score).max(end_score - LN_2);
            let suffix =
                (special(&specials, sequence_position - 1, C) + loop_score).max(end_score - LN_2);
            let prefix = special(&specials, sequence_position - 1, N) + loop_score;
            let begin = (prefix + move_score).max(join + move_score);
            set_special(&mut specials, sequence_position, J, join);
            set_special(&mut specials, sequence_position, C, suffix);
            set_special(&mut specials, sequence_position, N, prefix);
            set_special(&mut specials, sequence_position, B, begin);
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

    fn trace(&self, profile: &Profile<'_>, sequence: &[u8]) -> Vec<RawDomain> {
        #[derive(Clone, Copy)]
        enum State {
            C,
            E,
            Match,
            Insert,
            Delete,
            Begin,
            Join,
            Prefix,
        }

        let mut sequence_position = self.sequence_length;
        let mut model_position = 0;
        let mut state = State::C;
        let mut reversed_steps = Vec::new();
        let mut domains = Vec::new();

        loop {
            match state {
                State::C => {
                    let current = special(&self.specials, sequence_position, C);
                    if sequence_position > 0
                        && near(
                            current,
                            special(&self.specials, sequence_position - 1, C) + self.loop_score,
                        )
                    {
                        sequence_position -= 1;
                    } else {
                        state = State::E;
                    }
                }
                State::E => {
                    let current = special(&self.specials, sequence_position, E);
                    model_position = (1..=self.model_length)
                        .rev()
                        .find(|position| {
                            near(
                                current,
                                matrix_cell(
                                    &self.cells,
                                    self.model_length,
                                    sequence_position,
                                    *position,
                                    MATCH,
                                ),
                            )
                        })
                        .unwrap_or_else(|| {
                            (1..=self.model_length)
                                .max_by(|left, right| {
                                    matrix_cell(
                                        &self.cells,
                                        self.model_length,
                                        sequence_position,
                                        *left,
                                        MATCH,
                                    )
                                    .total_cmp(&matrix_cell(
                                        &self.cells,
                                        self.model_length,
                                        sequence_position,
                                        *right,
                                        MATCH,
                                    ))
                                })
                                .unwrap_or(1)
                        });
                    state = State::Match;
                }
                State::Match => {
                    let sequence_index = sequence_position - 1;
                    reversed_steps.push(TraceStep {
                        state: TraceState::Match,
                        model_position: model_position as u16,
                        sequence_index: Some(sequence_index),
                    });
                    let current = matrix_cell(
                        &self.cells,
                        self.model_length,
                        sequence_position,
                        model_position,
                        MATCH,
                    );
                    let emission =
                        profile.match_score(model_position, usize::from(sequence[sequence_index]));
                    let from_begin = special(&self.specials, sequence_position - 1, B)
                        + profile.local_entry_score(model_position)
                        + emission;
                    let next = if near(current, from_begin) {
                        State::Begin
                    } else if model_position > 1
                        && near(
                            current,
                            matrix_cell(
                                &self.cells,
                                self.model_length,
                                sequence_position - 1,
                                model_position - 1,
                                MATCH,
                            ) + profile.transition_score(model_position - 1, MM)
                                + emission,
                        )
                    {
                        State::Match
                    } else if model_position > 1
                        && near(
                            current,
                            matrix_cell(
                                &self.cells,
                                self.model_length,
                                sequence_position - 1,
                                model_position - 1,
                                INSERT,
                            ) + profile.transition_score(model_position - 1, IM)
                                + emission,
                        )
                    {
                        State::Insert
                    } else {
                        State::Delete
                    };
                    model_position -= 1;
                    sequence_position -= 1;
                    state = next;
                }
                State::Insert => {
                    let sequence_index = sequence_position - 1;
                    reversed_steps.push(TraceStep {
                        state: TraceState::Insert,
                        model_position: model_position as u16,
                        sequence_index: Some(sequence_index),
                    });
                    let current = matrix_cell(
                        &self.cells,
                        self.model_length,
                        sequence_position,
                        model_position,
                        INSERT,
                    );
                    let from_match = matrix_cell(
                        &self.cells,
                        self.model_length,
                        sequence_position - 1,
                        model_position,
                        MATCH,
                    ) + profile.transition_score(model_position, MI);
                    state = if near(current, from_match) {
                        State::Match
                    } else {
                        State::Insert
                    };
                    sequence_position -= 1;
                }
                State::Delete => {
                    reversed_steps.push(TraceStep {
                        state: TraceState::Delete,
                        model_position: model_position as u16,
                        sequence_index: None,
                    });
                    let current = matrix_cell(
                        &self.cells,
                        self.model_length,
                        sequence_position,
                        model_position,
                        DELETE,
                    );
                    let from_match = matrix_cell(
                        &self.cells,
                        self.model_length,
                        sequence_position,
                        model_position - 1,
                        MATCH,
                    ) + profile.transition_score(model_position - 1, MD);
                    state = if near(current, from_match) {
                        State::Match
                    } else {
                        State::Delete
                    };
                    model_position -= 1;
                }
                State::Begin => {
                    if !reversed_steps.is_empty() {
                        reversed_steps.reverse();
                        let start = reversed_steps
                            .iter()
                            .filter_map(|step| step.sequence_index)
                            .min()
                            .unwrap_or(sequence_position);
                        let end = reversed_steps
                            .iter()
                            .filter_map(|step| step.sequence_index)
                            .max()
                            .map_or(start, |index| index + 1);
                        domains.push(RawDomain {
                            steps: std::mem::take(&mut reversed_steps),
                            start,
                            end,
                        });
                    }
                    let current = special(&self.specials, sequence_position, B);
                    state = if near(
                        current,
                        special(&self.specials, sequence_position, N) + self.move_score,
                    ) {
                        State::Prefix
                    } else {
                        State::Join
                    };
                }
                State::Join => {
                    let current = special(&self.specials, sequence_position, J);
                    if sequence_position > 0
                        && near(
                            current,
                            special(&self.specials, sequence_position - 1, J) + self.loop_score,
                        )
                    {
                        sequence_position -= 1;
                    } else {
                        state = State::E;
                    }
                }
                State::Prefix => {
                    break;
                }
            }
        }

        domains.reverse();
        domains
    }
}

fn cell(row: &[f32], _model_length: usize, model_position: usize, state: usize) -> f32 {
    row[model_position * STATE_COUNT + state]
}

fn set_cell(
    row: &mut [f32],
    _model_length: usize,
    model_position: usize,
    state: usize,
    value: f32,
) {
    row[model_position * STATE_COUNT + state] = value;
}

fn matrix_cell(
    cells: &[f32],
    model_length: usize,
    sequence_position: usize,
    model_position: usize,
    state: usize,
) -> f32 {
    cells[((sequence_position * (model_length + 1) + model_position) * STATE_COUNT) + state]
}

fn set_matrix_cell(
    cells: &mut [f32],
    model_length: usize,
    sequence_position: usize,
    model_position: usize,
    state: usize,
    value: f32,
) {
    let index = ((sequence_position * (model_length + 1) + model_position) * STATE_COUNT) + state;
    cells[index] = value;
}

fn special(specials: &[f32], sequence_position: usize, state: usize) -> f32 {
    specials[sequence_position * SPECIAL_COUNT + state]
}

fn set_special(specials: &mut [f32], sequence_position: usize, state: usize, value: f32) {
    specials[sequence_position * SPECIAL_COUNT + state] = value;
}

fn max4(a: f32, b: f32, c: f32, d: f32) -> f32 {
    a.max(b).max(c).max(d)
}

fn logsum(a: f32, b: f32) -> f32 {
    if a.is_infinite() && a.is_sign_negative() {
        return b;
    }
    if b.is_infinite() && b.is_sign_negative() {
        return a;
    }
    let maximum = a.max(b);
    maximum + (-(a - b).abs()).exp().ln_1p()
}

fn logsum3(a: f32, b: f32, c: f32) -> f32 {
    logsum(logsum(a, b), c)
}

fn logsum4(a: f32, b: f32, c: f32, d: f32) -> f32 {
    logsum(logsum(a, b), logsum(c, d))
}

fn near(left: f32, right: f32) -> bool {
    if left == right {
        return true;
    }
    if !left.is_finite() || !right.is_finite() {
        return false;
    }
    (left - right).abs() <= 1.0e-4 * (1.0 + left.abs().max(right.abs()))
}

#[cfg(test)]
mod tests {
    use crate::{embedded_profiles, sequence::residue_index};

    use super::*;

    const VH: &str = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";

    fn encode(sequence: &str) -> Vec<u8> {
        sequence
            .bytes()
            .map(|residue| residue_index(residue).unwrap())
            .collect()
    }

    #[test]
    fn human_heavy_profile_recovers_a_full_domain_trace() {
        let database = embedded_profiles().unwrap();
        let profile = &database.profiles()[0];
        let sequence = encode(VH);
        let domains = viterbi_domains(profile, &sequence);
        assert_eq!(domains.len(), 1);
        assert!(domains[0].start <= 2);
        assert!(domains[0].end >= 115);
        assert!(forward_bit_score(profile, &sequence) > 100.0);
    }
}
