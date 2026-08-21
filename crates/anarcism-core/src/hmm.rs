use std::sync::OnceLock;

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
// HMMER3's p7_FLogsum discretizes log(1 + exp(-difference)) at
// 0.001-natural-log intervals over [0, 16). Values at a difference of
// 15.7 or greater return the larger operand directly.
const LOGSUM_SCALE: f32 = 1_000.0;
const LOGSUM_TABLE_SIZE: usize = 16_000;
static LOGSUM_LOOKUP: OnceLock<Box<[f32]>> = OnceLock::new();

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

#[derive(Clone, Copy, Debug)]
pub(crate) struct PosteriorRegion {
    pub start: usize,
    pub end: usize,
    pub is_multidomain: bool,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct DomainDefinition {
    pub start: usize,
    pub end: usize,
    pub trace_null2_bias: Option<f32>,
}

#[derive(Clone, Copy, Debug)]
struct SampleSegment {
    sample_index: usize,
    sequence_start: usize,
    sequence_end: usize,
    model_start: u16,
    model_end: u16,
}

#[derive(Clone, Copy, Debug)]
struct SampleCluster {
    start: usize,
    end: usize,
    probability: f32,
}

struct TraceEnsemble {
    clusters: Vec<SampleCluster>,
    log_null2_odds: Vec<f32>,
}

pub(crate) fn viterbi_domains(profile: &Profile<'_>, sequence: &[u8]) -> Vec<RawDomain> {
    if sequence.is_empty() {
        return Vec::new();
    }
    let matrix = ViterbiMatrix::fill(profile, sequence);
    matrix.trace(profile, sequence)
}

/// Reproduce HMMER's inexpensive posterior region trigger. The initial
/// search model is multihit-local; its B/E/core occupancies decide whether a
/// contiguous hit should be isolated as one envelope or sent to domain
/// clustering before unihit rescoring and OA alignment.
pub(crate) fn posterior_regions(profile: &Profile<'_>, sequence: &[u8]) -> Vec<PosteriorRegion> {
    if sequence.is_empty() {
        return Vec::new();
    }
    let forward = ForwardMatrix::fill_multihit(profile, sequence, sequence.len());
    let backward = BackwardMatrix::fill_multihit(profile, sequence, sequence.len());
    let overall_score = forward.score();
    let (loop_score, _) = length_scores(sequence.len());
    let mut begin_total = vec![0.0_f32; sequence.len() + 1];
    let mut end_total = vec![0.0_f32; sequence.len() + 1];
    let mut model_occupancy = vec![0.0_f32; sequence.len() + 1];

    for sequence_position in 1..=sequence.len() {
        let begin_log_probability = special(&forward.specials, sequence_position - 1, B)
            + special(&backward.specials, sequence_position - 1, B)
            - overall_score;
        let end_log_probability = special(&forward.specials, sequence_position, E)
            + special(&backward.specials, sequence_position, E)
            - overall_score;
        begin_total[sequence_position] =
            begin_total[sequence_position - 1] + f64::from(begin_log_probability).exp() as f32;
        end_total[sequence_position] =
            end_total[sequence_position - 1] + f64::from(end_log_probability).exp() as f32;

        let flank_probability: f32 = [N, J, C]
            .into_iter()
            .map(|state| {
                (special(&forward.specials, sequence_position - 1, state)
                    + special(&backward.specials, sequence_position, state)
                    + loop_score
                    - overall_score)
                    .exp()
            })
            .sum();
        model_occupancy[sequence_position] = (1.0 - flank_probability).clamp(0.0, 1.0);
    }

    let mut regions = Vec::new();
    let mut region_start = None;
    let mut triggered = false;
    for sequence_position in 1..=sequence.len() {
        if !triggered {
            let begin_probability =
                begin_total[sequence_position] - begin_total[sequence_position - 1];
            if model_occupancy[sequence_position] - begin_probability < 0.10
                || region_start.is_none()
            {
                region_start = Some(sequence_position);
            }
            if model_occupancy[sequence_position] >= 0.25 {
                triggered = true;
            }
        } else {
            let end_probability = end_total[sequence_position] - end_total[sequence_position - 1];
            if model_occupancy[sequence_position] - end_probability < 0.10 {
                let start = region_start.unwrap_or(sequence_position);
                let maximum_split_expectation = (start..=sequence_position)
                    .map(|split| {
                        let ends_before = end_total[split] - end_total[start - 1];
                        let begins_after = begin_total[sequence_position] - begin_total[split - 1];
                        ends_before.min(begins_after)
                    })
                    .fold(0.0_f32, f32::max);
                regions.push(PosteriorRegion {
                    start: start - 1,
                    end: sequence_position,
                    is_multidomain: maximum_split_expectation >= 0.20,
                });
                region_start = None;
                triggered = false;
            }
        }
    }
    regions
}

pub(crate) fn define_domain(
    profile: &Profile<'_>,
    seed: &RawDomain,
    sequence: &[u8],
) -> DomainDefinition {
    let regions = posterior_regions(profile, sequence);
    let Some(region) = regions.iter().max_by_key(|region| {
        region
            .end
            .min(seed.end)
            .saturating_sub(region.start.max(seed.start))
    }) else {
        return DomainDefinition {
            start: seed.start,
            end: seed.end,
            trace_null2_bias: None,
        };
    };

    if !region.is_multidomain {
        return DomainDefinition {
            start: region.start,
            end: region.end,
            trace_null2_bias: None,
        };
    }

    let ensemble = trace_ensemble(profile, sequence, *region);
    let selected = ensemble.clusters.iter().max_by(|left, right| {
        let left_overlap = left
            .end
            .min(seed.end)
            .saturating_sub(left.start.max(seed.start));
        let right_overlap = right
            .end
            .min(seed.end)
            .saturating_sub(right.start.max(seed.start));
        left_overlap.cmp(&right_overlap).then_with(|| {
            left.probability
                .total_cmp(&right.probability)
                .then_with(|| right.start.cmp(&left.start))
        })
    });
    let Some(selected) = selected else {
        return DomainDefinition {
            start: region.start,
            end: region.end,
            trace_null2_bias: None,
        };
    };
    let correction: f32 = ensemble.log_null2_odds
        [selected.start - region.start..selected.end - region.start]
        .iter()
        .sum();
    DomainDefinition {
        start: selected.start,
        end: selected.end,
        trace_null2_bias: Some(logsum(0.0, -(256.0_f32).ln() + correction) / LN_2),
    }
}

fn trace_ensemble(
    profile: &Profile<'_>,
    sequence: &[u8],
    region: PosteriorRegion,
) -> TraceEnsemble {
    let optimized = trace_ensemble_with_order(profile, sequence, region, true);
    let optimized_span = optimized
        .clusters
        .iter()
        .map(|cluster| cluster.end - cluster.start)
        .max()
        .unwrap_or(0);
    if optimized_span == region.end - region.start {
        return optimized;
    }

    // The optimized HMMER implementation enumerates striped SIMD lanes in
    // stochastic traceback. Tiny arithmetic differences between its
    // probability-space matrix and this scalar log-space matrix can move a
    // seeded roll across a boundary. Probe the mathematically equivalent
    // generic ordering as well and retain the better-supported envelope.
    let generic = trace_ensemble_with_order(profile, sequence, region, false);
    let generic_span = generic
        .clusters
        .iter()
        .map(|cluster| cluster.end - cluster.start)
        .max()
        .unwrap_or(0);
    if generic_span > optimized_span {
        generic
    } else {
        optimized
    }
}

fn trace_ensemble_with_order(
    profile: &Profile<'_>,
    sequence: &[u8],
    region: PosteriorRegion,
    optimized_end_order: bool,
) -> TraceEnsemble {
    const SAMPLE_COUNT: usize = 200;

    let region_sequence = &sequence[region.start..region.end];
    let forward = ForwardMatrix::fill_multihit(profile, region_sequence, sequence.len());
    let mut rng = EaselFastRng::new(42);
    let mut choices = vec![f32::NEG_INFINITY; 2 * forward.model_length + 1];
    let mut segments = Vec::with_capacity(SAMPLE_COUNT * 2);
    let mut null2_odds_sum = vec![0.0_f32; region_sequence.len()];

    for sample_index in 0..SAMPLE_COUNT {
        let domains = forward.stochastic_trace(
            profile,
            region.start,
            &mut rng,
            &mut choices,
            optimized_end_order,
        );
        let mut position = region.start;
        for domain in &domains {
            let Some((sequence_start, sequence_end, model_start, model_end)) = match_bounds(domain)
            else {
                continue;
            };
            segments.push(SampleSegment {
                sample_index,
                sequence_start,
                sequence_end,
                model_start,
                model_end,
            });

            for sequence_index in position..=sequence_start.min(region.end - 1) {
                null2_odds_sum[sequence_index - region.start] += 1.0;
            }
            let odds = trace_null2_odds(profile, domain);
            for sequence_index in sequence_start.saturating_add(1)..=sequence_end {
                if sequence_index >= region.end {
                    break;
                }
                null2_odds_sum[sequence_index - region.start] +=
                    odds[usize::from(sequence[sequence_index])];
            }
            position = sequence_end.saturating_add(1).min(region.end);
        }
        for sequence_index in position..region.end {
            null2_odds_sum[sequence_index - region.start] += 1.0;
        }
    }

    for odds in &mut null2_odds_sum {
        *odds = (*odds / SAMPLE_COUNT as f32).ln();
    }
    TraceEnsemble {
        clusters: cluster_sample_segments(&segments, SAMPLE_COUNT),
        log_null2_odds: null2_odds_sum,
    }
}

fn match_bounds(domain: &RawDomain) -> Option<(usize, usize, u16, u16)> {
    let first = domain
        .steps
        .iter()
        .find(|step| step.state == TraceState::Match && step.sequence_index.is_some())?;
    let last = domain
        .steps
        .iter()
        .rev()
        .find(|step| step.state == TraceState::Match && step.sequence_index.is_some())?;
    Some((
        first.sequence_index?,
        last.sequence_index?,
        first.model_position,
        last.model_position,
    ))
}

fn trace_null2_odds(profile: &Profile<'_>, domain: &RawDomain) -> [f32; 20] {
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
    let mut odds = [1.0_f32; 20];
    if emitted == 0.0 {
        return odds;
    }
    for (residue, value) in odds.iter_mut().enumerate() {
        let mut weighted = insert_usage;
        for (model_position, usage) in match_usage.iter().enumerate().skip(1) {
            weighted += usage * profile.match_score(model_position, residue).exp();
        }
        *value = weighted / emitted;
    }
    odds
}

fn cluster_sample_segments(segments: &[SampleSegment], sample_count: usize) -> Vec<SampleCluster> {
    if segments.is_empty() {
        return Vec::new();
    }
    let mut parents: Vec<usize> = (0..segments.len()).collect();
    for left in 0..segments.len() {
        for right in left + 1..segments.len() {
            if sample_segments_link(segments[left], segments[right]) {
                union_sets(&mut parents, left, right);
            }
        }
    }

    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut roots = Vec::new();
    for segment_index in 0..segments.len() {
        let root = find_set(&mut parents, segment_index);
        let group_index = roots
            .iter()
            .position(|existing| *existing == root)
            .unwrap_or_else(|| {
                roots.push(root);
                groups.push(Vec::new());
                groups.len() - 1
            });
        groups[group_index].push(segment_index);
    }

    let mut clusters = Vec::new();
    for group in groups {
        let mut represented_samples = vec![false; sample_count];
        for &segment_index in &group {
            represented_samples[segments[segment_index].sample_index] = true;
        }
        let represented = represented_samples
            .into_iter()
            .filter(|value| *value)
            .count();
        if represented * 4 < sample_count {
            continue;
        }
        let endpoint_threshold = ((represented as f32) * 0.02).ceil() as usize;
        let start = select_endpoint(
            group.iter().map(|index| segments[*index].sequence_start),
            endpoint_threshold,
            true,
        );
        let end_inclusive = select_endpoint(
            group.iter().map(|index| segments[*index].sequence_end),
            endpoint_threshold,
            false,
        );
        let model_start = select_endpoint(
            group
                .iter()
                .map(|index| usize::from(segments[*index].model_start)),
            endpoint_threshold,
            true,
        );
        let model_end = select_endpoint(
            group
                .iter()
                .map(|index| usize::from(segments[*index].model_end)),
            endpoint_threshold,
            false,
        );
        if start <= end_inclusive && model_start <= model_end {
            clusters.push(SampleCluster {
                start,
                end: end_inclusive + 1,
                probability: represented as f32 / sample_count as f32,
            });
        }
    }
    clusters.sort_by_key(|cluster| cluster.start);

    let mut dominated = vec![false; clusters.len()];
    for left in 0..clusters.len() {
        for right in left + 1..clusters.len() {
            let overlap = clusters[left]
                .end
                .min(clusters[right].end)
                .saturating_sub(clusters[left].start.max(clusters[right].start));
            let shorter = (clusters[left].end - clusters[left].start)
                .min(clusters[right].end - clusters[right].start);
            if overlap * 5 >= shorter * 4 {
                if clusters[left].probability > clusters[right].probability {
                    dominated[right] = true;
                } else {
                    dominated[left] = true;
                }
            }
        }
    }
    clusters
        .into_iter()
        .enumerate()
        .filter_map(|(index, cluster)| (!dominated[index]).then_some(cluster))
        .collect()
}

fn select_endpoint(
    coordinates: impl Iterator<Item = usize>,
    threshold: usize,
    select_leftmost: bool,
) -> usize {
    let coordinates: Vec<usize> = coordinates.collect();
    let minimum = coordinates.iter().copied().min().unwrap_or(0);
    let maximum = coordinates.iter().copied().max().unwrap_or(minimum);
    let mut counts = vec![0_usize; maximum - minimum + 1];
    for coordinate in coordinates {
        counts[coordinate - minimum] += 1;
    }
    if select_leftmost {
        counts
            .iter()
            .position(|count| *count >= threshold)
            .map_or(minimum, |offset| minimum + offset)
    } else {
        counts
            .iter()
            .rposition(|count| *count >= threshold)
            .map_or(maximum, |offset| minimum + offset)
    }
}

fn sample_segments_link(left: SampleSegment, right: SampleSegment) -> bool {
    let sequence_overlap = left
        .sequence_end
        .min(right.sequence_end)
        .saturating_sub(left.sequence_start.max(right.sequence_start))
        + 1;
    let shorter_sequence = (left.sequence_end - left.sequence_start + 1)
        .min(right.sequence_end - right.sequence_start + 1);
    if sequence_overlap * 5 < shorter_sequence * 4 {
        return false;
    }

    let model_overlap = usize::from(left.model_end.min(right.model_end))
        .saturating_sub(usize::from(left.model_start.max(right.model_start)));
    let shorter_model = usize::from(left.model_end - left.model_start + 1)
        .min(usize::from(right.model_end - right.model_start + 1));
    if model_overlap * 5 < shorter_model * 4 {
        return false;
    }

    let left_start_diagonal = left.sequence_start as isize - left.model_start as isize;
    let right_start_diagonal = right.sequence_start as isize - right.model_start as isize;
    let left_end_diagonal = left.sequence_end as isize - left.model_end as isize;
    let right_end_diagonal = right.sequence_end as isize - right.model_end as isize;
    (left_start_diagonal - right_start_diagonal).abs() <= 4
        || (left_end_diagonal - right_end_diagonal).abs() <= 4
}

fn find_set(parents: &mut [usize], index: usize) -> usize {
    if parents[index] != index {
        parents[index] = find_set(parents, parents[index]);
    }
    parents[index]
}

fn union_sets(parents: &mut [usize], left: usize, right: usize) {
    let left_root = find_set(parents, left);
    let right_root = find_set(parents, right);
    if left_root != right_root {
        parents[right_root] = left_root;
    }
}

struct EaselFastRng {
    state: u32,
}

impl EaselFastRng {
    fn new(seed: u32) -> Self {
        let mut state = mix_three(seed, 87_654_321, 12_345_678);
        if state == 0 {
            state = 42;
        }
        Self { state }
    }

    fn random(&mut self) -> f64 {
        self.state = self.state.wrapping_mul(69_069).wrapping_add(1);
        f64::from(self.state) / 4_294_967_296.0
    }

    fn choose_log_weights(&mut self, weights: &[f32]) -> usize {
        let maximum = weights.iter().copied().fold(f32::NEG_INFINITY, f32::max);
        let mut norm = 0.0_f64;
        for weight in weights {
            norm += f64::from((*weight - maximum).exp());
        }
        let roll = self.random();
        let mut cumulative = 0.0_f64;
        for (index, weight) in weights.iter().enumerate() {
            cumulative += f64::from((*weight - maximum).exp());
            if roll < cumulative / norm {
                return index;
            }
        }
        weights.len().saturating_sub(1)
    }
}

fn mix_three(mut a: u32, mut b: u32, mut c: u32) -> u32 {
    a = a.wrapping_sub(b).wrapping_sub(c) ^ (c >> 13);
    b = b.wrapping_sub(c).wrapping_sub(a) ^ a.wrapping_shl(8);
    c = c.wrapping_sub(a).wrapping_sub(b) ^ (b >> 13);
    a = a.wrapping_sub(b).wrapping_sub(c) ^ (c >> 12);
    b = b.wrapping_sub(c).wrapping_sub(a) ^ a.wrapping_shl(16);
    c = c.wrapping_sub(a).wrapping_sub(b) ^ (b >> 5);
    a = a.wrapping_sub(b).wrapping_sub(c) ^ (c >> 3);
    b = b.wrapping_sub(c).wrapping_sub(a) ^ a.wrapping_shl(10);
    c.wrapping_sub(a).wrapping_sub(b) ^ (b >> 15)
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

/// HMMER3's isolated-domain Forward score, before the null2 composition
/// correction. Domain envelopes are rescored in local, unihit mode while the
/// profile and null length models remain configured for the complete target.
#[cfg(test)]
pub(crate) fn domain_forward_bit_score(
    profile: &Profile<'_>,
    envelope: &[u8],
    target_length: usize,
) -> f32 {
    if envelope.is_empty() || target_length == 0 || envelope.len() > target_length {
        return f32::NEG_INFINITY;
    }
    let model_length = profile.consensus().len();
    let envelope_length = envelope.len();
    let (loop_score, move_score) = unihit_length_scores(target_length);

    let mut previous = vec![f32::NEG_INFINITY; (model_length + 1) * STATE_COUNT];
    let mut current = vec![f32::NEG_INFINITY; (model_length + 1) * STATE_COUNT];
    let mut previous_special = [f32::NEG_INFINITY; SPECIAL_COUNT];
    previous_special[N] = 0.0;
    previous_special[B] = move_score;

    for residue in envelope {
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

        // Isolated envelopes are always rescored in unihit mode: E can only
        // move to C, and J is unreachable.
        special[J] = f32::NEG_INFINITY;
        special[C] = logsum(previous_special[C] + loop_score, special[E]);
        special[N] = previous_special[N] + loop_score;
        special[B] = special[N] + move_score;

        std::mem::swap(&mut previous, &mut current);
        previous_special = special;
    }

    let raw_score = previous_special[C] + move_score;
    let null_score = null_one_score(target_length);
    let outside_envelope = (target_length - envelope_length) as f32
        * (target_length as f32 / (target_length as f32 + 3.0)).ln();
    (raw_score + outside_envelope - null_score) / LN_2
}

/// Rescore one Viterbi envelope using HMMER's isolated-domain workflow:
/// unihit Forward/Backward and null2 by posterior expectation.
pub(crate) fn domain_bit_score(profile: &Profile<'_>, seed: &RawDomain, sequence: &[u8]) -> f32 {
    domain_bit_score_with_null2(profile, seed, sequence, None)
}

pub(crate) fn domain_bit_score_with_null2(
    profile: &Profile<'_>,
    seed: &RawDomain,
    sequence: &[u8],
    trace_null2_bias: Option<f32>,
) -> f32 {
    let envelope = &sequence[seed.start..seed.end];
    let forward = ForwardMatrix::fill(profile, envelope, sequence.len());
    let env_score = forward.score();
    let null2_bias = trace_null2_bias.unwrap_or_else(|| {
        let backward = BackwardMatrix::fill(profile, envelope, sequence.len());
        let posterior = PosteriorMatrix::decode(profile, &forward, backward);
        posterior.null2_bias_bits(profile, envelope)
    });
    let outside_envelope = (sequence.len() - envelope.len()) as f32
        * (sequence.len() as f32 / (sequence.len() as f32 + 3.0)).ln();
    (env_score + outside_envelope - null_one_score(sequence.len())) / LN_2 - null2_bias
}

/// Replace a Viterbi display with HMMER's posterior optimal-accuracy trace.
pub(crate) fn realign_domain(
    profile: &Profile<'_>,
    seed: &RawDomain,
    sequence: &[u8],
) -> RawDomain {
    let model_end = seed
        .steps
        .iter()
        .rev()
        .find(|step| step.sequence_index.is_some())
        .map_or(0, |step| step.model_position);
    let envelope_start = seed.start.saturating_sub(20);
    let envelope_end =
        if (104..120).contains(&model_end) && sequence.len().saturating_sub(seed.end) > 30 {
            // This also gives ANARCI's J-region recovery enough sequence to bridge
            // an unusually long CDR3. Normal domains use a tightly bounded HMMER-
            // like envelope so adjacent domains cannot be conflated.
            sequence.len()
        } else {
            seed.end.saturating_add(20).min(sequence.len())
        };
    let envelope = &sequence[envelope_start..envelope_end];
    realign_envelope_slice(
        profile,
        seed,
        sequence,
        envelope_start,
        envelope_end,
        envelope,
    )
}

pub(crate) fn realign_envelope(
    profile: &Profile<'_>,
    envelope: &RawDomain,
    sequence: &[u8],
) -> RawDomain {
    let envelope_sequence = &sequence[envelope.start..envelope.end];
    realign_envelope_slice(
        profile,
        envelope,
        sequence,
        envelope.start,
        envelope.end,
        envelope_sequence,
    )
}

fn realign_envelope_slice(
    profile: &Profile<'_>,
    fallback: &RawDomain,
    sequence: &[u8],
    envelope_start: usize,
    envelope_end: usize,
    envelope: &[u8],
) -> RawDomain {
    let forward = ForwardMatrix::fill(profile, envelope, sequence.len());
    let backward = BackwardMatrix::fill(profile, envelope, sequence.len());
    let posterior = PosteriorMatrix::decode(profile, &forward, backward);
    posterior
        .optimal_accuracy_trace(profile, envelope_start)
        .unwrap_or_else(|| {
            let mut domain = fallback.clone();
            domain.start = envelope_start;
            domain.end = envelope_end;
            domain
        })
}

/// ANARCI performs a second, permissive profile scan after IMGT 104 when a
/// single-domain alignment leaves a long suffix without a J region. Rebuild
/// the intervening CDR3 and splice the recovered J trace onto the V trace.
pub(crate) fn recover_long_cdr3(
    profiles: &[Profile<'_>],
    domain: &RawDomain,
    sequence: &[u8],
) -> RawDomain {
    let Some(last_emitting) = domain
        .steps
        .iter()
        .rev()
        .find(|step| step.sequence_index.is_some())
    else {
        return domain.clone();
    };
    let Some(last_index) = last_emitting.sequence_index else {
        return domain.clone();
    };
    if last_emitting.model_position >= 120 || last_index.saturating_add(30) >= sequence.len() {
        return domain.clone();
    }
    let Some(cysteine_step) = domain.steps.iter().position(|step| {
        step.state == TraceState::Match
            && step.model_position == 104
            && step.sequence_index.is_some()
    }) else {
        return domain.clone();
    };
    let cysteine_index = domain.steps[cysteine_step]
        .sequence_index
        .unwrap_or(last_index);
    let suffix_offset = cysteine_index + 1;
    let suffix = &sequence[suffix_offset..];

    let mut best: Option<(f32, &str, RawDomain)> = None;
    for profile in profiles {
        for seed in viterbi_domains(profile, suffix) {
            if seed.end <= seed.start {
                continue;
            }
            let score = domain_bit_score(profile, &seed, suffix);
            if score < 10.0 {
                continue;
            }
            let aligned = realign_domain(profile, &seed, suffix);
            let first_model = aligned
                .steps
                .iter()
                .find(|step| step.sequence_index.is_some())
                .map_or(u16::MAX, |step| step.model_position);
            let last_model = aligned
                .steps
                .iter()
                .rev()
                .find(|step| step.sequence_index.is_some())
                .map_or(0, |step| step.model_position);
            if first_model > 117 || last_model < 126 {
                continue;
            }
            let significance = profile.forward_lambda() * (score - profile.forward_tau());
            let replace = best.as_ref().is_none_or(|(current, name, _)| {
                significance > *current || (significance == *current && profile.name() < *name)
            });
            if replace {
                best = Some((significance, profile.name(), aligned));
            }
        }
    }
    let Some((_, _, recovered)) = best else {
        return domain.clone();
    };

    let mut j_region: Vec<TraceStep> = recovered
        .steps
        .into_iter()
        .filter(|step| step.model_position >= 117 && step.sequence_index.is_some())
        .map(|mut step| {
            step.sequence_index = step.sequence_index.map(|index| index + suffix_offset);
            step
        })
        .collect();
    let Some(first_j_index) = j_region.first().and_then(|step| step.sequence_index) else {
        return domain.clone();
    };
    if first_j_index <= cysteine_index {
        return domain.clone();
    }

    let mut steps = domain.steps[..=cysteine_step].to_vec();
    let mut next_model = 105_u16;
    for sequence_index in suffix_offset..first_j_index {
        let (state, model_position) = if next_model >= 116 {
            (TraceState::Insert, 116)
        } else {
            let position = next_model;
            next_model += 1;
            (TraceState::Match, position)
        };
        steps.push(TraceStep {
            state,
            model_position,
            sequence_index: Some(sequence_index),
        });
    }
    steps.append(&mut j_region);
    let start = steps
        .iter()
        .filter_map(|step| step.sequence_index)
        .min()
        .unwrap_or(domain.start);
    let end = steps
        .iter()
        .filter_map(|step| step.sequence_index)
        .max()
        .map_or(domain.end, |index| index + 1);
    RawDomain { steps, start, end }
}

struct ForwardMatrix {
    model_length: usize,
    sequence_length: usize,
    cells: Vec<f32>,
    specials: Vec<f32>,
    loop_score: f32,
    move_score: f32,
    end_loop_score: f32,
    end_move_score: f32,
}

impl ForwardMatrix {
    fn fill(profile: &Profile<'_>, sequence: &[u8], target_length: usize) -> Self {
        Self::fill_with_mode(profile, sequence, target_length, false)
    }

    fn fill_multihit(profile: &Profile<'_>, sequence: &[u8], target_length: usize) -> Self {
        Self::fill_with_mode(profile, sequence, target_length, true)
    }

    fn fill_with_mode(
        profile: &Profile<'_>,
        sequence: &[u8],
        target_length: usize,
        multihit: bool,
    ) -> Self {
        let model_length = profile.consensus().len();
        let sequence_length = sequence.len();
        let logsum = Logsum::shared();
        let row_width = (model_length + 1) * STATE_COUNT;
        let mut cells = vec![f32::NEG_INFINITY; (sequence_length + 1) * row_width];
        let mut specials = vec![f32::NEG_INFINITY; (sequence_length + 1) * SPECIAL_COUNT];
        let (loop_score, move_score) = if multihit {
            length_scores(target_length)
        } else {
            unihit_length_scores(target_length)
        };
        let end_loop_score = if multihit { -LN_2 } else { f32::NEG_INFINITY };
        let end_move_score = if multihit { -LN_2 } else { 0.0 };
        set_special(&mut specials, 0, N, 0.0);
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
            let residue = usize::from(sequence[sequence_position - 1]);
            let begin_state_score = special(&specials, sequence_position - 1, B);
            let mut end_score = f32::NEG_INFINITY;
            for model_position in 1..=model_length {
                let (from_match, from_insert, from_delete) = if model_position > 1 {
                    let previous_transitions = &transition_score_rows[model_position - 2];
                    (
                        previous_row[model_position - 1][MATCH] + previous_transitions[MM],
                        previous_row[model_position - 1][INSERT] + previous_transitions[IM],
                        previous_row[model_position - 1][DELETE] + previous_transitions[DM],
                    )
                } else {
                    (f32::NEG_INFINITY, f32::NEG_INFINITY, f32::NEG_INFINITY)
                };
                let from_begin = begin_state_score + local_entry_scores[model_position - 1];
                let match_score = logsum.four(from_match, from_insert, from_begin, from_delete)
                    + match_score_rows[model_position - 1][residue];
                current_row[model_position][MATCH] = match_score;

                let insert_score = if model_position < model_length {
                    let transitions = &transition_score_rows[model_position - 1];
                    logsum.pair(
                        previous_row[model_position][MATCH] + transitions[MI],
                        previous_row[model_position][INSERT] + transitions[II],
                    )
                } else {
                    f32::NEG_INFINITY
                };
                current_row[model_position][INSERT] = insert_score;

                let delete_score = if model_position > 1 {
                    let previous_transitions = &transition_score_rows[model_position - 2];
                    logsum.pair(
                        current_row[model_position - 1][MATCH] + previous_transitions[MD],
                        current_row[model_position - 1][DELETE] + previous_transitions[DD],
                    )
                } else {
                    f32::NEG_INFINITY
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
            end_loop_score,
            end_move_score,
        }
    }

    fn score(&self) -> f32 {
        special(&self.specials, self.sequence_length, C) + self.move_score
    }

    fn stochastic_trace(
        &self,
        profile: &Profile<'_>,
        sequence_offset: usize,
        rng: &mut EaselFastRng,
        choices: &mut [f32],
        optimized_end_order: bool,
    ) -> Vec<RawDomain> {
        #[derive(Clone, Copy)]
        enum State {
            C,
            End,
            Match,
            Insert,
            Delete,
            Begin,
            Join,
            Prefix,
        }

        let mut sequence_position = self.sequence_length;
        let mut model_position = 0_usize;
        let mut state = State::C;
        let mut reversed_steps = Vec::new();
        let mut reversed_domains = Vec::new();

        loop {
            match state {
                State::C => {
                    let weights = [
                        if sequence_position > 0 {
                            special(&self.specials, sequence_position - 1, C) + self.loop_score
                        } else {
                            f32::NEG_INFINITY
                        },
                        special(&self.specials, sequence_position, E) + self.end_move_score,
                    ];
                    if rng.choose_log_weights(&weights) == 0 {
                        sequence_position = sequence_position.saturating_sub(1);
                    } else {
                        state = State::End;
                    }
                }
                State::End => {
                    if optimized_end_order {
                        // hmmscan's E traceback visits four striped M lanes,
                        // then four D lanes, for each SIMD segment.
                        let segment_count = self.model_length.div_ceil(4);
                        let choice_count = segment_count * 8;
                        for segment in 0..segment_count {
                            for lane in 0..4 {
                                let position = lane * segment_count + segment + 1;
                                choices[segment * 8 + lane] = if position <= self.model_length {
                                    matrix_cell(
                                        &self.cells,
                                        self.model_length,
                                        sequence_position,
                                        position,
                                        MATCH,
                                    )
                                } else {
                                    f32::NEG_INFINITY
                                };
                                choices[segment * 8 + 4 + lane] =
                                    if position > 1 && position <= self.model_length {
                                        matrix_cell(
                                            &self.cells,
                                            self.model_length,
                                            sequence_position,
                                            position,
                                            DELETE,
                                        )
                                    } else {
                                        f32::NEG_INFINITY
                                    };
                            }
                        }
                        let selected = rng.choose_log_weights(&choices[..choice_count]);
                        let segment = selected / 8;
                        let lane = selected % 4;
                        model_position = lane * segment_count + segment + 1;
                        state = if selected % 8 < 4 {
                            State::Match
                        } else {
                            State::Delete
                        };
                    } else {
                        choices.fill(f32::NEG_INFINITY);
                        for (position, choice) in choices
                            .iter_mut()
                            .take(self.model_length + 1)
                            .enumerate()
                            .skip(1)
                        {
                            *choice = matrix_cell(
                                &self.cells,
                                self.model_length,
                                sequence_position,
                                position,
                                MATCH,
                            );
                        }
                        for (choice_index, choice) in
                            choices.iter_mut().enumerate().skip(self.model_length + 2)
                        {
                            let position = choice_index - self.model_length;
                            *choice = matrix_cell(
                                &self.cells,
                                self.model_length,
                                sequence_position,
                                position,
                                DELETE,
                            );
                        }
                        let selected = rng.choose_log_weights(choices);
                        if selected <= self.model_length {
                            model_position = selected;
                            state = State::Match;
                        } else {
                            model_position = selected - self.model_length;
                            state = State::Delete;
                        }
                    }
                }
                State::Match => {
                    if sequence_position == 0 || model_position == 0 {
                        break;
                    }
                    reversed_steps.push(TraceStep {
                        state: TraceState::Match,
                        model_position: model_position as u16,
                        sequence_index: Some(sequence_offset + sequence_position - 1),
                    });
                    let previous_sequence = sequence_position - 1;
                    let previous_model = model_position - 1;
                    let weights = [
                        special(&self.specials, previous_sequence, B)
                            + profile.local_entry_score(model_position),
                        if model_position > 1 {
                            matrix_cell(
                                &self.cells,
                                self.model_length,
                                previous_sequence,
                                previous_model,
                                MATCH,
                            ) + profile.transition_score(previous_model, MM)
                        } else {
                            f32::NEG_INFINITY
                        },
                        if model_position > 1 {
                            matrix_cell(
                                &self.cells,
                                self.model_length,
                                previous_sequence,
                                previous_model,
                                INSERT,
                            ) + profile.transition_score(previous_model, IM)
                        } else {
                            f32::NEG_INFINITY
                        },
                        if model_position > 1 {
                            matrix_cell(
                                &self.cells,
                                self.model_length,
                                previous_sequence,
                                previous_model,
                                DELETE,
                            ) + profile.transition_score(previous_model, DM)
                        } else {
                            f32::NEG_INFINITY
                        },
                    ];
                    state = match rng.choose_log_weights(&weights) {
                        0 => State::Begin,
                        1 => State::Match,
                        2 => State::Insert,
                        _ => State::Delete,
                    };
                    sequence_position = previous_sequence;
                    model_position = previous_model;
                }
                State::Delete => {
                    if model_position <= 1 {
                        break;
                    }
                    reversed_steps.push(TraceStep {
                        state: TraceState::Delete,
                        model_position: model_position as u16,
                        sequence_index: None,
                    });
                    let previous_model = model_position - 1;
                    let weights = [
                        matrix_cell(
                            &self.cells,
                            self.model_length,
                            sequence_position,
                            previous_model,
                            MATCH,
                        ) + profile.transition_score(previous_model, MD),
                        matrix_cell(
                            &self.cells,
                            self.model_length,
                            sequence_position,
                            previous_model,
                            DELETE,
                        ) + profile.transition_score(previous_model, DD),
                    ];
                    state = if rng.choose_log_weights(&weights) == 0 {
                        State::Match
                    } else {
                        State::Delete
                    };
                    model_position = previous_model;
                }
                State::Insert => {
                    if sequence_position == 0 || model_position == 0 {
                        break;
                    }
                    reversed_steps.push(TraceStep {
                        state: TraceState::Insert,
                        model_position: model_position as u16,
                        sequence_index: Some(sequence_offset + sequence_position - 1),
                    });
                    let previous_sequence = sequence_position - 1;
                    let weights = [
                        matrix_cell(
                            &self.cells,
                            self.model_length,
                            previous_sequence,
                            model_position,
                            MATCH,
                        ) + profile.transition_score(model_position, MI),
                        matrix_cell(
                            &self.cells,
                            self.model_length,
                            previous_sequence,
                            model_position,
                            INSERT,
                        ) + profile.transition_score(model_position, II),
                    ];
                    state = if rng.choose_log_weights(&weights) == 0 {
                        State::Match
                    } else {
                        State::Insert
                    };
                    sequence_position = previous_sequence;
                }
                State::Begin => {
                    if !reversed_steps.is_empty() {
                        reversed_steps.reverse();
                        let start = reversed_steps
                            .iter()
                            .filter_map(|step| step.sequence_index)
                            .min()
                            .unwrap_or(sequence_offset);
                        let end = reversed_steps
                            .iter()
                            .filter_map(|step| step.sequence_index)
                            .max()
                            .map_or(start, |index| index + 1);
                        reversed_domains.push(RawDomain {
                            steps: std::mem::take(&mut reversed_steps),
                            start,
                            end,
                        });
                    }
                    let weights = [
                        special(&self.specials, sequence_position, N) + self.move_score,
                        special(&self.specials, sequence_position, J) + self.move_score,
                    ];
                    state = if rng.choose_log_weights(&weights) == 0 {
                        State::Prefix
                    } else {
                        State::Join
                    };
                }
                State::Join => {
                    let weights = [
                        if sequence_position > 0 {
                            special(&self.specials, sequence_position - 1, J) + self.loop_score
                        } else {
                            f32::NEG_INFINITY
                        },
                        special(&self.specials, sequence_position, E) + self.end_loop_score,
                    ];
                    if rng.choose_log_weights(&weights) == 0 {
                        sequence_position = sequence_position.saturating_sub(1);
                    } else {
                        state = State::End;
                    }
                }
                State::Prefix => {
                    if sequence_position == 0 {
                        break;
                    }
                    sequence_position -= 1;
                }
            }
        }
        reversed_domains.reverse();
        reversed_domains
    }
}

struct BackwardMatrix {
    cells: Vec<f32>,
    specials: Vec<f32>,
}

impl BackwardMatrix {
    fn fill(profile: &Profile<'_>, sequence: &[u8], target_length: usize) -> Self {
        Self::fill_with_mode(profile, sequence, target_length, false)
    }

    fn fill_multihit(profile: &Profile<'_>, sequence: &[u8], target_length: usize) -> Self {
        Self::fill_with_mode(profile, sequence, target_length, true)
    }

    fn fill_with_mode(
        profile: &Profile<'_>,
        sequence: &[u8],
        target_length: usize,
        multihit: bool,
    ) -> Self {
        let model_length = profile.consensus().len();
        let sequence_length = sequence.len();
        let logsum = Logsum::shared();
        let row_width = (model_length + 1) * STATE_COUNT;
        let mut cells = vec![f32::NEG_INFINITY; (sequence_length + 1) * row_width];
        let mut specials = vec![f32::NEG_INFINITY; (sequence_length + 1) * SPECIAL_COUNT];
        let (loop_score, move_score) = if multihit {
            length_scores(target_length)
        } else {
            unihit_length_scores(target_length)
        };
        let end_loop_score = if multihit { -LN_2 } else { f32::NEG_INFINITY };
        let end_move_score = if multihit { -LN_2 } else { 0.0 };

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
            let match_score = logsum.pair(terminal_end, next_delete + transitions[MD]);
            let delete_score = logsum.pair(terminal_end, next_delete + transitions[DD]);
            terminal_row[model_position][MATCH] = match_score;
            terminal_row[model_position][DELETE] = delete_score;
        }

        for sequence_position in (1..sequence_length).rev() {
            let next_start = (sequence_position + 1) * state_row_width;
            let (through_current, next_and_later) = state_cells.split_at_mut(next_start);
            let current_start = sequence_position * state_row_width;
            let current_row = &mut through_current[current_start..current_start + state_row_width];
            let next_row = &next_and_later[..state_row_width];
            let next_residue = usize::from(sequence[sequence_position]);
            let mut begin_score = f32::NEG_INFINITY;
            for model_position in 1..=model_length {
                begin_score = logsum.pair(
                    begin_score,
                    next_row[model_position][MATCH]
                        + local_entry_scores[model_position - 1]
                        + match_score_rows[model_position - 1][next_residue],
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
                let next_match = next_row[model_position + 1][MATCH]
                    + match_score_rows[model_position][next_residue];
                let match_score = logsum.four(
                    next_match + transitions[MM],
                    next_row[model_position][INSERT] + transitions[MI],
                    end_score,
                    current_row[model_position + 1][DELETE] + transitions[MD],
                );
                let insert_score = logsum.pair(
                    next_match + transitions[IM],
                    next_row[model_position][INSERT] + transitions[II],
                );
                let delete_score = logsum.three(
                    next_match + transitions[DM],
                    current_row[model_position + 1][DELETE] + transitions[DD],
                    end_score,
                );
                current_row[model_position][MATCH] = match_score;
                current_row[model_position][INSERT] = insert_score;
                current_row[model_position][DELETE] = delete_score;
            }
        }

        if sequence_length > 0 {
            let first_row = &state_cells[state_row_width..state_row_width * 2];
            let first_residue = usize::from(sequence[0]);
            let mut begin_score = f32::NEG_INFINITY;
            for model_position in 1..=model_length {
                begin_score = logsum.pair(
                    begin_score,
                    first_row[model_position][MATCH]
                        + local_entry_scores[model_position - 1]
                        + match_score_rows[model_position - 1][first_residue],
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

struct PosteriorMatrix {
    model_length: usize,
    sequence_length: usize,
    cells: Vec<f32>,
    specials: Vec<f32>,
}

impl PosteriorMatrix {
    fn decode(
        _profile: &Profile<'_>,
        forward: &ForwardMatrix,
        mut backward: BackwardMatrix,
    ) -> Self {
        let overall_score = forward.score();
        let (loop_score, _) = unihit_length_scores_from_move(forward.move_score);
        for sequence_position in 1..=forward.sequence_length {
            let mut denominator = 0.0_f32;
            for model_position in 1..=forward.model_length {
                let match_probability = (matrix_cell(
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
                ) - overall_score)
                    .exp();
                let insert_probability = if model_position < forward.model_length {
                    (matrix_cell(
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
                    ) - overall_score)
                        .exp()
                } else {
                    0.0
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
                    0.0,
                );
                denominator += match_probability + insert_probability;
            }

            for state in [N, J, C] {
                let probability = (special(&forward.specials, sequence_position - 1, state)
                    + special(&backward.specials, sequence_position, state)
                    + loop_score
                    - overall_score)
                    .exp();
                set_special(
                    &mut backward.specials,
                    sequence_position,
                    state,
                    probability,
                );
                denominator += probability;
            }
            set_special(&mut backward.specials, sequence_position, E, 0.0);
            set_special(&mut backward.specials, sequence_position, B, 0.0);

            if denominator > 0.0 && denominator.is_finite() {
                let scale = denominator.recip();
                for model_position in 1..=forward.model_length {
                    for state in [MATCH, INSERT] {
                        let probability = matrix_cell(
                            &backward.cells,
                            forward.model_length,
                            sequence_position,
                            model_position,
                            state,
                        );
                        set_matrix_cell(
                            &mut backward.cells,
                            forward.model_length,
                            sequence_position,
                            model_position,
                            state,
                            probability * scale,
                        );
                    }
                }
                for state in [N, J, C] {
                    let probability = special(&backward.specials, sequence_position, state) * scale;
                    set_special(
                        &mut backward.specials,
                        sequence_position,
                        state,
                        probability,
                    );
                }
            }
        }
        backward.cells[..(forward.model_length + 1) * STATE_COUNT].fill(0.0);
        backward.specials[..SPECIAL_COUNT].fill(0.0);

        Self {
            model_length: forward.model_length,
            sequence_length: forward.sequence_length,
            cells: backward.cells,
            specials: backward.specials,
        }
    }

    fn null2_bias_bits(&self, profile: &Profile<'_>, sequence: &[u8]) -> f32 {
        if self.sequence_length == 0 {
            return 0.0;
        }
        let mut match_usage = vec![0.0_f32; self.model_length + 1];
        let mut insert_usage = 0.0_f32;
        let mut flank_usage = 0.0_f32;
        for sequence_position in 1..=self.sequence_length {
            for (model_position, usage) in match_usage.iter_mut().enumerate().skip(1) {
                *usage += matrix_cell(
                    &self.cells,
                    self.model_length,
                    sequence_position,
                    model_position,
                    MATCH,
                );
                insert_usage += matrix_cell(
                    &self.cells,
                    self.model_length,
                    sequence_position,
                    model_position,
                    INSERT,
                );
            }
            flank_usage += [N, J, C]
                .into_iter()
                .map(|state| special(&self.specials, sequence_position, state))
                .sum::<f32>();
        }
        let normalizer = (self.sequence_length as f32).recip();
        let mut odds = [0.0_f32; 20];
        for (residue, value) in odds.iter_mut().enumerate() {
            let mut weighted = insert_usage + flank_usage;
            for (model_position, usage) in match_usage.iter().enumerate().skip(1) {
                weighted += usage * profile.match_score(model_position, residue).exp();
            }
            *value = weighted * normalizer;
        }
        let correction: f32 = sequence
            .iter()
            .map(|residue| odds[usize::from(*residue)].ln())
            .sum();
        logsum(0.0, -(256.0_f32).ln() + correction) / LN_2
    }

    fn optimal_accuracy_trace(
        &self,
        profile: &Profile<'_>,
        sequence_offset: usize,
    ) -> Option<RawDomain> {
        let mut oa = OaMatrix::fill(profile, self);
        oa.trace(profile, self, sequence_offset)
    }
}

fn unihit_length_scores_from_move(move_score: f32) -> (f32, f32) {
    let move_probability = move_score.exp();
    ((1.0 - move_probability).ln(), move_score)
}

struct OaMatrix {
    model_length: usize,
    sequence_length: usize,
    cells: Vec<f32>,
    specials: Vec<f32>,
}

impl OaMatrix {
    fn fill(profile: &Profile<'_>, posterior: &PosteriorMatrix) -> Self {
        let model_length = posterior.model_length;
        let sequence_length = posterior.sequence_length;
        let row_width = (model_length + 1) * STATE_COUNT;
        let mut cells = vec![f32::NEG_INFINITY; (sequence_length + 1) * row_width];
        let mut specials = vec![f32::NEG_INFINITY; (sequence_length + 1) * SPECIAL_COUNT];
        set_special(&mut specials, 0, N, 0.0);
        set_special(&mut specials, 0, B, 0.0);

        for sequence_position in 1..=sequence_length {
            let mut end_score = f32::NEG_INFINITY;
            for model_position in 1..=model_length {
                let posterior_match = matrix_cell(
                    &posterior.cells,
                    model_length,
                    sequence_position,
                    model_position,
                    MATCH,
                );
                let mut match_score = if profile.local_entry_score(model_position).is_finite() {
                    special(&specials, sequence_position - 1, B)
                } else {
                    f32::NEG_INFINITY
                };
                if model_position > 1 {
                    for (state, transition) in [(MATCH, MM), (INSERT, IM), (DELETE, DM)] {
                        if profile
                            .transition_score(model_position - 1, transition)
                            .is_finite()
                        {
                            match_score = match_score.max(matrix_cell(
                                &cells,
                                model_length,
                                sequence_position - 1,
                                model_position - 1,
                                state,
                            ));
                        }
                    }
                }
                match_score += posterior_match;
                set_matrix_cell(
                    &mut cells,
                    model_length,
                    sequence_position,
                    model_position,
                    MATCH,
                    match_score,
                );

                let insert_score = if model_position < model_length {
                    let posterior_insert = matrix_cell(
                        &posterior.cells,
                        model_length,
                        sequence_position,
                        model_position,
                        INSERT,
                    );
                    let from_match = if profile.transition_score(model_position, MI).is_finite() {
                        matrix_cell(
                            &cells,
                            model_length,
                            sequence_position - 1,
                            model_position,
                            MATCH,
                        )
                    } else {
                        f32::NEG_INFINITY
                    };
                    let from_insert = if profile.transition_score(model_position, II).is_finite() {
                        matrix_cell(
                            &cells,
                            model_length,
                            sequence_position - 1,
                            model_position,
                            INSERT,
                        )
                    } else {
                        f32::NEG_INFINITY
                    };
                    from_match.max(from_insert) + posterior_insert
                } else {
                    f32::NEG_INFINITY
                };
                set_matrix_cell(
                    &mut cells,
                    model_length,
                    sequence_position,
                    model_position,
                    INSERT,
                    insert_score,
                );

                let delete_score = if model_position > 1 {
                    let from_match = if profile.transition_score(model_position - 1, MD).is_finite()
                    {
                        matrix_cell(
                            &cells,
                            model_length,
                            sequence_position,
                            model_position - 1,
                            MATCH,
                        )
                    } else {
                        f32::NEG_INFINITY
                    };
                    let from_delete =
                        if profile.transition_score(model_position - 1, DD).is_finite() {
                            matrix_cell(
                                &cells,
                                model_length,
                                sequence_position,
                                model_position - 1,
                                DELETE,
                            )
                        } else {
                            f32::NEG_INFINITY
                        };
                    from_match.max(from_delete)
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
                // In a local profile, an internal match may exit to E but an
                // internal delete may not. D_M is the sole delete-state exit;
                // this is the slightly asymmetric recurrence used by
                // p7_GOptimalAccuracy().
                end_score = end_score.max(match_score);
                if model_position == model_length {
                    end_score = end_score.max(delete_score);
                }
            }

            set_special(&mut specials, sequence_position, E, end_score);
            set_special(&mut specials, sequence_position, J, f32::NEG_INFINITY);
            let previous_suffix = special(&specials, sequence_position - 1, C)
                + special(&posterior.specials, sequence_position, C);
            set_special(
                &mut specials,
                sequence_position,
                C,
                previous_suffix.max(end_score),
            );
            let prefix = special(&specials, sequence_position - 1, N)
                + special(&posterior.specials, sequence_position, N);
            set_special(&mut specials, sequence_position, N, prefix);
            set_special(&mut specials, sequence_position, B, prefix);
        }

        Self {
            model_length,
            sequence_length,
            cells,
            specials,
        }
    }

    fn trace(
        &mut self,
        profile: &Profile<'_>,
        posterior: &PosteriorMatrix,
        sequence_offset: usize,
    ) -> Option<RawDomain> {
        #[derive(Clone, Copy)]
        enum State {
            C,
            E,
            Match,
            Insert,
            Delete,
            Begin,
            Prefix,
        }

        let mut sequence_position = self.sequence_length;
        let mut model_position = 0_usize;
        let mut state = State::C;
        let mut reversed_steps = Vec::new();

        loop {
            match state {
                State::C => {
                    let from_suffix = if sequence_position > 0 {
                        special(&self.specials, sequence_position - 1, C)
                            + special(&posterior.specials, sequence_position, C)
                    } else {
                        f32::NEG_INFINITY
                    };
                    let from_end = special(&self.specials, sequence_position, E);
                    if from_suffix > from_end {
                        sequence_position = sequence_position.checked_sub(1)?;
                    } else {
                        state = State::E;
                    }
                }
                State::E => {
                    let mut maximum = f32::NEG_INFINITY;
                    let mut selected = None;
                    for position in 1..=self.model_length {
                        let match_score = matrix_cell(
                            &self.cells,
                            self.model_length,
                            sequence_position,
                            position,
                            MATCH,
                        );
                        if match_score >= maximum {
                            maximum = match_score;
                            selected = Some((State::Match, position));
                        }
                        let delete_score = matrix_cell(
                            &self.cells,
                            self.model_length,
                            sequence_position,
                            position,
                            DELETE,
                        );
                        if delete_score > maximum {
                            maximum = delete_score;
                            selected = Some((State::Delete, position));
                        }
                    }
                    let (selected_state, selected_position) = selected?;
                    if !maximum.is_finite() {
                        return None;
                    }
                    model_position = selected_position;
                    state = selected_state;
                }
                State::Match => {
                    if sequence_position == 0 || model_position == 0 {
                        return None;
                    }
                    reversed_steps.push(TraceStep {
                        state: TraceState::Match,
                        model_position: model_position as u16,
                        sequence_index: Some(sequence_offset + sequence_position - 1),
                    });
                    let previous_sequence = sequence_position - 1;
                    let previous_model = model_position - 1;
                    let mut selected = State::Match;
                    let mut maximum = if model_position > 1
                        && profile.transition_score(model_position - 1, MM).is_finite()
                    {
                        matrix_cell(
                            &self.cells,
                            self.model_length,
                            previous_sequence,
                            previous_model,
                            MATCH,
                        )
                    } else {
                        f32::NEG_INFINITY
                    };
                    if model_position > 1
                        && profile.transition_score(model_position - 1, IM).is_finite()
                    {
                        let candidate = matrix_cell(
                            &self.cells,
                            self.model_length,
                            previous_sequence,
                            previous_model,
                            INSERT,
                        );
                        if candidate > maximum {
                            maximum = candidate;
                            selected = State::Insert;
                        }
                    }
                    if model_position > 1
                        && profile.transition_score(model_position - 1, DM).is_finite()
                    {
                        let candidate = matrix_cell(
                            &self.cells,
                            self.model_length,
                            previous_sequence,
                            previous_model,
                            DELETE,
                        );
                        if candidate > maximum {
                            maximum = candidate;
                            selected = State::Delete;
                        }
                    }
                    if profile.local_entry_score(model_position).is_finite() {
                        let candidate = special(&self.specials, previous_sequence, B);
                        if candidate > maximum {
                            selected = State::Begin;
                        }
                    }
                    sequence_position = previous_sequence;
                    model_position = previous_model;
                    state = selected;
                }
                State::Insert => {
                    if sequence_position == 0 || model_position == 0 {
                        return None;
                    }
                    reversed_steps.push(TraceStep {
                        state: TraceState::Insert,
                        model_position: model_position as u16,
                        sequence_index: Some(sequence_offset + sequence_position - 1),
                    });
                    let previous_sequence = sequence_position - 1;
                    let from_match = if profile.transition_score(model_position, MI).is_finite() {
                        matrix_cell(
                            &self.cells,
                            self.model_length,
                            previous_sequence,
                            model_position,
                            MATCH,
                        )
                    } else {
                        f32::NEG_INFINITY
                    };
                    let from_insert = if profile.transition_score(model_position, II).is_finite() {
                        matrix_cell(
                            &self.cells,
                            self.model_length,
                            previous_sequence,
                            model_position,
                            INSERT,
                        )
                    } else {
                        f32::NEG_INFINITY
                    };
                    state = if from_match >= from_insert {
                        State::Match
                    } else {
                        State::Insert
                    };
                    sequence_position = previous_sequence;
                }
                State::Delete => {
                    if model_position <= 1 {
                        return None;
                    }
                    reversed_steps.push(TraceStep {
                        state: TraceState::Delete,
                        model_position: model_position as u16,
                        sequence_index: None,
                    });
                    let previous_model = model_position - 1;
                    let from_match = if profile.transition_score(previous_model, MD).is_finite() {
                        matrix_cell(
                            &self.cells,
                            self.model_length,
                            sequence_position,
                            previous_model,
                            MATCH,
                        )
                    } else {
                        f32::NEG_INFINITY
                    };
                    let from_delete = if profile.transition_score(previous_model, DD).is_finite() {
                        matrix_cell(
                            &self.cells,
                            self.model_length,
                            sequence_position,
                            previous_model,
                            DELETE,
                        )
                    } else {
                        f32::NEG_INFINITY
                    };
                    state = if from_match >= from_delete {
                        State::Match
                    } else {
                        State::Delete
                    };
                    model_position = previous_model;
                }
                State::Begin => state = State::Prefix,
                State::Prefix => {
                    if sequence_position == 0 {
                        break;
                    }
                    sequence_position -= 1;
                }
            }
        }

        reversed_steps.reverse();
        let start = reversed_steps
            .iter()
            .filter_map(|step| step.sequence_index)
            .min()?;
        let end = reversed_steps
            .iter()
            .filter_map(|step| step.sequence_index)
            .max()?
            + 1;
        Some(RawDomain {
            steps: reversed_steps,
            start,
            end,
        })
    }
}

fn null_one_score(length: usize) -> f32 {
    let p1 = length as f32 / (length as f32 + 1.0);
    length as f32 * p1.ln() + (1.0 - p1).ln()
}

fn length_scores(length: usize) -> (f32, f32) {
    let move_probability = 3.0 / (length as f32 + 3.0);
    ((1.0 - move_probability).ln(), move_probability.ln())
}

fn unihit_length_scores(length: usize) -> (f32, f32) {
    let move_probability = 2.0 / (length as f32 + 2.0);
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

#[cfg(test)]
fn cell(row: &[f32], _model_length: usize, model_position: usize, state: usize) -> f32 {
    row[model_position * STATE_COUNT + state]
}

#[cfg(test)]
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

#[derive(Clone, Copy)]
struct Logsum(&'static [f32]);

impl Logsum {
    #[inline]
    fn shared() -> Self {
        Self(logsum_lookup())
    }

    #[inline(always)]
    fn pair(self, a: f32, b: f32) -> f32 {
        // A right-hand -infinity produces an infinite difference and returns
        // `a` through the cutoff below, so only the left sentinel needs an
        // explicit branch. Profile scores and DP state never contain NaN.
        if a == f32::NEG_INFINITY {
            return b;
        }
        let maximum = a.max(b);
        let difference = (a - b).abs();
        if difference >= 15.7 {
            maximum
        } else {
            maximum + self.0[(difference * LOGSUM_SCALE) as usize]
        }
    }

    #[inline(always)]
    fn three(self, a: f32, b: f32, c: f32) -> f32 {
        self.pair(self.pair(a, b), c)
    }

    #[inline(always)]
    fn four(self, a: f32, b: f32, c: f32, d: f32) -> f32 {
        self.pair(self.pair(a, b), self.pair(c, d))
    }
}

#[inline]
fn logsum(a: f32, b: f32) -> f32 {
    Logsum::shared().pair(a, b)
}

#[inline]
fn logsum_lookup() -> &'static [f32] {
    LOGSUM_LOOKUP.get_or_init(|| {
        (0..LOGSUM_TABLE_SIZE)
            .map(|index| (1.0 + (-(index as f64) / f64::from(LOGSUM_SCALE)).exp()).ln() as f32)
            .collect::<Vec<_>>()
            .into_boxed_slice()
    })
}

#[cfg(test)]
fn logsum3(a: f32, b: f32, c: f32) -> f32 {
    Logsum::shared().three(a, b, c)
}

#[cfg(test)]
fn logsum4(a: f32, b: f32, c: f32, d: f32) -> f32 {
    Logsum::shared().four(a, b, c, d)
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
        assert!(domain_forward_bit_score(profile, &sequence, sequence.len()) > 100.0);
    }

    #[test]
    fn logsum_matches_hmmer_table_behavior() {
        assert_eq!(logsum(0.0, f32::NEG_INFINITY), 0.0);
        assert_eq!(logsum(f32::NEG_INFINITY, 0.0), 0.0);
        assert_eq!(
            logsum(f32::NEG_INFINITY, f32::NEG_INFINITY),
            f32::NEG_INFINITY
        );
        assert_eq!(logsum(0.0, -15.7), 0.0);
        assert_eq!(logsum(0.0, 0.0), std::f32::consts::LN_2);

        let left = -0.4_f32;
        let right = -0.500_9_f32;
        let exact = left.max(right) + (-(left - right).abs()).exp().ln_1p();
        assert!((logsum(left, right) - exact).abs() < 0.001);
        assert_eq!(logsum_lookup().len(), LOGSUM_TABLE_SIZE);
    }
}
