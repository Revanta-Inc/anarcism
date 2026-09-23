use std::collections::{BTreeMap, BTreeSet};
#[cfg(not(target_family = "wasm"))]
use std::num::NonZeroUsize;
#[cfg(not(target_family = "wasm"))]
use std::sync::Mutex;

use crate::components::connected_components;
use crate::germlines::assign_closest_germline;
use crate::hmm::{
    DomainAlignment, RawDomain, SEQUENCE_BATCH_LANES, SequenceBatchViterbiWorkspace, define_domain,
    domain_score_components_with_posterior, msv_filter_passes, realign_domain, realign_envelope,
    recover_long_cdr3, trace_emission_score, ungapped_filter_score, viterbi_domains,
};
use crate::numbering::number_imgt;
use crate::sequence::{MAX_SEQUENCE_LENGTH, NormalizedSequence, normalize_sequence};
use crate::{
    ChainType, DomainResult, Error, ErrorCode, NumberingOptions, PairValidationOptions,
    PairValidationResult, Profile, ProfileDatabase, ProfileHit, Result, SequenceInput,
    SequenceResult, embedded_profiles,
};

const EXACT_PROFILE_RESERVE: usize = 7;
const MIN_EXACT_PROFILES: usize = 10;
const BATCH_PROFILE_RESERVE: usize = 4;
const MIN_BATCH_PROFILES: usize = 7;
const BATCH_TASK_SEQUENCE_LIMIT: usize = SEQUENCE_BATCH_LANES * 32;

struct SearchPlan<'a> {
    options: &'a NumberingOptions,
    allowed_profiles: u64,
    exact_profile_budget: usize,
}

impl<'a> SearchPlan<'a> {
    fn new(database: &ProfileDatabase<'_>, options: &'a NumberingOptions) -> Result<Self> {
        validate_options(options)?;
        let maximum_alternatives = database.profiles().len().saturating_sub(1);
        if options.alternative_hit_count > maximum_alternatives {
            return Err(Error::sequence(
                ErrorCode::InvalidOptions,
                format!("alternativeHitCount cannot exceed {maximum_alternatives}"),
                None,
            ));
        }
        let allowed_profiles = database
            .profiles()
            .iter()
            .enumerate()
            .filter(|(_, profile)| profile_is_allowed(profile, options))
            .fold(0_u64, |mask, (index, _)| mask | (1_u64 << index));
        Ok(Self {
            options,
            allowed_profiles,
            exact_profile_budget: options
                .alternative_hit_count
                .saturating_add(EXACT_PROFILE_RESERVE)
                .max(MIN_EXACT_PROFILES),
        })
    }

    fn allows_profile(&self, index: usize) -> bool {
        self.allowed_profiles & (1_u64 << index) != 0
    }

    fn initial_trace_budget(&self, sequence_length: usize, available_profiles: usize) -> usize {
        self.exact_profile_budget
            .saturating_mul(estimated_domain_count(sequence_length))
            .min(available_profiles)
    }

    fn batch_trace_budget(&self, sequence_length: usize, available_profiles: usize) -> usize {
        self.options
            .alternative_hit_count
            .saturating_add(BATCH_PROFILE_RESERVE)
            .max(MIN_BATCH_PROFILES)
            .saturating_mul(estimated_domain_count(sequence_length))
            .min(available_profiles)
    }
}

#[derive(Debug)]
struct Candidate<'a> {
    profile: &'a Profile<'a>,
    domain: RawDomain,
    bit_score: f32,
    bias: f32,
    query_start: usize,
    query_end: usize,
    alignment_source: Option<DomainAlignment>,
    aligned_domain: Option<RawDomain>,
}

impl Candidate<'_> {
    fn ranking_score(&self) -> f32 {
        profile_significance(self.profile, self.bit_score)
    }
}

#[derive(Debug)]
struct PreparedRecord<'a> {
    index: usize,
    id: &'a str,
    sequence: NormalizedSequence,
    profiles: Vec<usize>,
}

#[derive(Debug)]
struct ProfileBatchTask {
    profile_index: usize,
    sequence_indices: Vec<usize>,
    #[cfg(not(target_family = "wasm"))]
    work: usize,
}

pub fn number_sequence(sequence: &str, options: &NumberingOptions) -> Result<SequenceResult> {
    number_sequence_with_id("sequence", sequence, options)
}

pub fn number_sequence_with_id(
    id: &str,
    sequence: &str,
    options: &NumberingOptions,
) -> Result<SequenceResult> {
    let database = embedded_profiles()?;
    let plan = SearchPlan::new(database, options)?;
    let normalized = normalize_sequence(id, sequence)?;
    number_normalized_sequence(id, normalized, &plan, database)
}

pub fn number_sequences(
    inputs: &[SequenceInput],
    options: &NumberingOptions,
) -> Result<Vec<SequenceResult>> {
    let database = embedded_profiles()?;
    let plan = SearchPlan::new(database, options)?;
    if inputs.is_empty() {
        return Ok(Vec::new());
    }

    let normalized: Vec<_> = inputs
        .iter()
        .map(|input| normalize_sequence(&input.id, &input.sequence))
        .collect::<Result<_>>()?;
    let records = prepare_records(inputs, normalized, database, &plan);
    number_prepared_records(records, &plan, database).map(strip_result_indices)
}

/// Number a native batch concurrently with an explicit worker limit.
#[cfg(not(target_family = "wasm"))]
pub fn number_sequences_parallel(
    inputs: &[SequenceInput],
    options: &NumberingOptions,
    worker_count: NonZeroUsize,
) -> Result<Vec<SequenceResult>> {
    let database = embedded_profiles()?;
    let plan = SearchPlan::new(database, options)?;
    if inputs.is_empty() {
        return Ok(Vec::new());
    }

    let normalized: Vec<_> = inputs
        .iter()
        .map(|input| normalize_sequence(&input.id, &input.sequence))
        .collect::<Result<_>>()?;
    let worker_count = worker_count.get().min(inputs.len());
    if worker_count <= 1 {
        let records = prepare_records(inputs, normalized, database, &plan);
        return number_prepared_records(records, &plan, database).map(strip_result_indices);
    }

    let initial_profiles =
        parallel_initial_profile_indices(database, &normalized, &plan, worker_count);
    let mut records = records_with_profiles(inputs, normalized, initial_profiles);
    if records.len() >= SEQUENCE_BATCH_LANES {
        let eligible_profiles =
            parallel_batch_eligible_profile_indices(database, &records, &plan, worker_count);
        replace_record_profiles(&mut records, eligible_profiles);
    }
    parallel_number_prepared_records(records, &plan, database, worker_count)
}

fn prepare_records<'a>(
    inputs: &'a [SequenceInput],
    normalized: Vec<NormalizedSequence>,
    database: &ProfileDatabase<'_>,
    plan: &SearchPlan<'_>,
) -> Vec<PreparedRecord<'a>> {
    let profiles = normalized
        .iter()
        .map(|sequence| initial_profile_indices(database, sequence, plan))
        .collect();
    records_with_profiles(inputs, normalized, profiles)
}

fn records_with_profiles<'a>(
    inputs: &'a [SequenceInput],
    normalized: Vec<NormalizedSequence>,
    profiles: Vec<Vec<usize>>,
) -> Vec<PreparedRecord<'a>> {
    inputs
        .iter()
        .zip(normalized)
        .zip(profiles)
        .enumerate()
        .map(|(index, ((input, sequence), profiles))| PreparedRecord {
            index,
            id: &input.id,
            sequence,
            profiles,
        })
        .collect()
}

fn replace_record_profiles(records: &mut [PreparedRecord<'_>], profiles: Vec<Vec<usize>>) {
    debug_assert_eq!(records.len(), profiles.len());
    for (record, profiles) in records.iter_mut().zip(profiles) {
        record.profiles = profiles;
    }
}

#[cfg(not(target_family = "wasm"))]
fn parallel_initial_profile_indices(
    database: &ProfileDatabase<'_>,
    sequences: &[NormalizedSequence],
    plan: &SearchPlan<'_>,
    worker_count: usize,
) -> Vec<Vec<usize>> {
    let weighted: Vec<_> = sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| (sequence.encoded.len(), index, index))
        .collect();
    let shards = distribute_by_weight(weighted, worker_count);

    std::thread::scope(|scope| {
        let workers: Vec<_> = shards
            .into_iter()
            .filter(|shard| !shard.is_empty())
            .map(|shard| {
                scope.spawn(move || {
                    shard
                        .into_iter()
                        .map(|index| {
                            (
                                index,
                                initial_profile_indices(database, &sequences[index], plan),
                            )
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        let mut initial: Vec<Option<Vec<usize>>> = (0..sequences.len()).map(|_| None).collect();
        for worker in workers {
            let shard = match worker.join() {
                Ok(shard) => shard,
                Err(payload) => std::panic::resume_unwind(payload),
            };
            for (index, profiles) in shard {
                initial[index] = Some(profiles);
            }
        }
        initial
            .into_iter()
            .map(|profiles| profiles.expect("prefilter returns every input exactly once"))
            .collect()
    })
}

#[cfg(not(target_family = "wasm"))]
fn distribute_by_weight<T>(
    mut weighted: Vec<(usize, usize, T)>,
    worker_count: usize,
) -> Vec<Vec<T>> {
    if weighted.is_empty() {
        return Vec::new();
    }
    weighted
        .sort_unstable_by(|left, right| right.0.cmp(&left.0).then_with(|| left.1.cmp(&right.1)));
    let worker_count = worker_count.min(weighted.len());
    let mut shards: Vec<Vec<T>> = (0..worker_count).map(|_| Vec::new()).collect();
    let mut loads = vec![0_usize; worker_count];
    for (work, _, task) in weighted {
        let lightest = (0..worker_count)
            .min_by_key(|worker| (loads[*worker], *worker))
            .expect("worker count is nonzero");
        loads[lightest] = loads[lightest].saturating_add(work);
        shards[lightest].push(task);
    }
    shards
}

#[cfg(not(target_family = "wasm"))]
fn parallel_batch_eligible_profile_indices(
    database: &ProfileDatabase<'_>,
    records: &[PreparedRecord<'_>],
    plan: &SearchPlan<'_>,
    worker_count: usize,
) -> Vec<Vec<usize>> {
    let weighted = profile_batch_tasks(database, records)
        .into_iter()
        .enumerate()
        .map(|(task_index, task)| (task.work, task_index, task))
        .collect();
    let shards = distribute_by_weight(weighted, worker_count);

    let rescored = std::thread::scope(|scope| {
        let workers: Vec<_> = shards
            .into_iter()
            .filter(|shard| !shard.is_empty())
            .map(|shard| {
                scope.spawn(move || {
                    let mut scores = Vec::new();
                    let mut workspace = SequenceBatchViterbiWorkspace::new(database.model_length());
                    for task in shard {
                        score_profile_batch_task(
                            database,
                            records,
                            &task,
                            &mut workspace,
                            &mut scores,
                        );
                    }
                    scores
                })
            })
            .collect();
        let mut rescored = rescore_buckets(records);
        for worker in workers {
            let scores = match worker.join() {
                Ok(scores) => scores,
                Err(payload) => std::panic::resume_unwind(payload),
            };
            for (sequence_index, score, profile_index) in scores {
                rescored[sequence_index].push((score, profile_index));
            }
        }
        rescored
    });
    finish_batch_eligible_profile_indices(database, records, rescored, plan)
}

#[cfg(not(target_family = "wasm"))]
fn profile_batch_work(profile: &Profile<'_>, sequence_length: usize, count: usize) -> usize {
    sequence_length
        .saturating_mul(count.div_ceil(SEQUENCE_BATCH_LANES))
        .saturating_mul(profile.consensus().len())
}

#[cfg(not(target_family = "wasm"))]
fn parallel_number_prepared_records(
    mut records: Vec<PreparedRecord<'_>>,
    plan: &SearchPlan<'_>,
    database: &ProfileDatabase<'_>,
    worker_count: usize,
) -> Result<Vec<SequenceResult>> {
    records.sort_unstable_by_key(|record| {
        (
            record
                .sequence
                .encoded
                .len()
                .saturating_mul(record.profiles.len().max(1)),
            record.index,
        )
    });
    let result_count = records.len();
    let queue = Mutex::new(records);
    std::thread::scope(|scope| {
        let workers: Vec<_> = (0..worker_count.min(result_count))
            .map(|_| {
                scope.spawn(|| -> Result<Vec<(usize, SequenceResult)>> {
                    let mut results = Vec::new();
                    loop {
                        let record = queue.lock().expect("scheduler queue is not poisoned").pop();
                        let Some(record) = record else {
                            break;
                        };
                        results.push(number_prepared_record(record, plan, database)?);
                    }
                    Ok(results)
                })
            })
            .collect();
        let mut results: Vec<Option<SequenceResult>> = (0..result_count).map(|_| None).collect();
        for worker in workers {
            let shard = match worker.join() {
                Ok(shard) => shard?,
                Err(payload) => std::panic::resume_unwind(payload),
            };
            for (index, result) in shard {
                results[index] = Some(result);
            }
        }
        Ok(results
            .into_iter()
            .map(|result| result.expect("scheduler returns every input exactly once"))
            .collect())
    })
}

pub fn number_fasta(fasta: &str, options: &NumberingOptions) -> Result<Vec<SequenceResult>> {
    let inputs = crate::fasta::parse_fasta(fasta)?;
    number_sequences(&inputs, options)
}

pub fn validate_antibody_pair(
    vh: &str,
    vl: &str,
    options: &PairValidationOptions,
) -> Result<PairValidationResult> {
    SearchPlan::new(embedded_profiles()?, &options.numbering)?;
    if options.start_max >= options.end_min || options.end_min > MAX_SEQUENCE_LENGTH {
        return Err(Error::sequence(
            ErrorCode::InvalidOptions,
            "pair validation requires startMax < endMin within the sequence-length limit",
            None,
        ));
    }

    let mut heavy_options = options.numbering.clone();
    heavy_options.allowed_chains =
        intersect_pair_chains(heavy_options.allowed_chains.as_deref(), &[ChainType::H]);
    let mut light_options = options.numbering.clone();
    light_options.allowed_chains = intersect_pair_chains(
        light_options.allowed_chains.as_deref(),
        &[ChainType::K, ChainType::L],
    );

    let heavy = if heavy_options
        .allowed_chains
        .as_ref()
        .is_some_and(Vec::is_empty)
    {
        None
    } else {
        let result = number_sequence_with_id("vh", vh, &heavy_options)?;
        single_pair_domain(result.domains)
    };
    let light = if light_options
        .allowed_chains
        .as_ref()
        .is_some_and(Vec::is_empty)
    {
        None
    } else {
        let result = number_sequence_with_id("vl", vl, &light_options)?;
        single_pair_domain(result.domains)
    };

    let mut errors = Vec::new();
    validate_pair_domain("VH", heavy.as_ref(), options, &mut errors);
    validate_pair_domain("VL", light.as_ref(), options, &mut errors);
    Ok(PairValidationResult {
        ok: errors.is_empty(),
        vh: heavy,
        vl: light,
        errors,
    })
}

fn number_prepared_records(
    mut records: Vec<PreparedRecord<'_>>,
    plan: &SearchPlan<'_>,
    database: &ProfileDatabase<'_>,
) -> Result<Vec<(usize, SequenceResult)>> {
    if records.len() >= SEQUENCE_BATCH_LANES {
        let eligible_profiles = batch_eligible_profile_indices(database, &records, plan);
        replace_record_profiles(&mut records, eligible_profiles);
    }
    records
        .into_iter()
        .map(|record| number_prepared_record(record, plan, database))
        .collect()
}

fn number_prepared_record(
    record: PreparedRecord<'_>,
    plan: &SearchPlan<'_>,
    database: &ProfileDatabase<'_>,
) -> Result<(usize, SequenceResult)> {
    number_normalized_with_profiles(
        record.id,
        record.sequence,
        plan,
        database,
        record
            .profiles
            .into_iter()
            .map(|profile_index| &database.profiles()[profile_index]),
    )
    .map(|result| (record.index, result))
}

fn strip_result_indices(indexed: Vec<(usize, SequenceResult)>) -> Vec<SequenceResult> {
    indexed.into_iter().map(|(_, result)| result).collect()
}

fn intersect_pair_chains(
    configured: Option<&[ChainType]>,
    required: &[ChainType],
) -> Option<Vec<ChainType>> {
    Some(
        required
            .iter()
            .copied()
            .filter(|chain| configured.is_none_or(|chains| chains.contains(chain)))
            .collect(),
    )
}

fn single_pair_domain(domains: Vec<DomainResult>) -> Option<DomainResult> {
    if domains.len() == 1 {
        domains.into_iter().next()
    } else {
        None
    }
}

fn validate_pair_domain(
    label: &str,
    domain: Option<&DomainResult>,
    options: &PairValidationOptions,
    errors: &mut Vec<String>,
) {
    let Some(domain) = domain else {
        errors.push(format!(
            "{label} must contain exactly one qualifying domain"
        ));
        return;
    };
    if domain.start > options.start_max {
        errors.push(format!(
            "{label} domain starts at {}; maximum is {}",
            domain.start, options.start_max
        ));
    }
    if domain.end < options.end_min {
        errors.push(format!(
            "{label} domain ends at {}; minimum is {}",
            domain.end, options.end_min
        ));
    }
}

fn number_normalized_sequence(
    id: &str,
    normalized: NormalizedSequence,
    plan: &SearchPlan<'_>,
    database: &ProfileDatabase<'_>,
) -> Result<SequenceResult> {
    let eligible = initial_profile_indices(database, &normalized, plan)
        .into_iter()
        .map(|profile_index| &database.profiles()[profile_index]);
    number_normalized_with_profiles(id, normalized, plan, database, eligible)
}

fn number_normalized_with_profiles<'a>(
    id: &str,
    normalized: NormalizedSequence,
    plan: &SearchPlan<'_>,
    database: &'a ProfileDatabase<'a>,
    eligible: impl IntoIterator<Item = &'a Profile<'a>>,
) -> Result<SequenceResult> {
    let candidates = collect_candidates(eligible, &normalized.encoded);
    let scored_clusters = score_candidate_clusters(
        candidates,
        &normalized.encoded,
        plan.exact_profile_budget,
        plan.options.min_bit_score,
    );
    let domains = render_domains(database, &normalized, plan.options, scored_clusters)?;

    Ok(SequenceResult {
        id: id.to_owned(),
        normalized_sequence: normalized.text,
        domains,
        warnings: normalized.warnings,
    })
}

fn collect_candidates<'a>(
    eligible: impl IntoIterator<Item = &'a Profile<'a>>,
    sequence: &[u8],
) -> Vec<Candidate<'a>> {
    let mut candidates = Vec::new();
    for profile in eligible {
        for domain in viterbi_domains(profile, sequence) {
            if domain.end <= domain.start {
                continue;
            }
            let query_start = domain.start;
            let query_end = domain.end;
            let bit_score = trace_emission_score(profile, &domain, sequence);
            candidates.push(Candidate {
                profile,
                domain,
                bit_score,
                bias: 0.0,
                query_start,
                query_end,
                alignment_source: None,
                aligned_domain: None,
            });
        }
    }
    candidates
}

fn score_candidate_clusters<'a>(
    candidates: Vec<Candidate<'a>>,
    sequence: &[u8],
    exact_profile_budget: usize,
    min_bit_score: f32,
) -> Vec<Vec<Candidate<'a>>> {
    let exact_candidates = cluster_candidates(candidates)
        .into_iter()
        .flat_map(|mut cluster| {
            // Reserve profiles cannot displace a requested exact result.
            cluster.truncate(exact_profile_budget);
            for candidate in &mut cluster {
                score_defined_domain_candidate(candidate, sequence);
            }
            cluster
        })
        .collect();

    // Several Viterbi seeds can resolve to one posterior envelope.
    let mut clusters = cluster_candidates(exact_candidates);

    let is_multidomain_target = clusters.len() > 1;
    if is_multidomain_target {
        // Use bit score inside the calibrated numerical uncertainty band.
        for cluster in &mut clusters {
            sort_multidomain_candidates(cluster);
        }
    }

    clusters
        .into_iter()
        .filter_map(|mut cluster| {
            cluster.retain(|candidate| candidate.bit_score >= min_bit_score);
            (!cluster.is_empty()).then_some(cluster)
        })
        .collect()
}

fn render_domains(
    database: &ProfileDatabase<'_>,
    normalized: &NormalizedSequence,
    options: &NumberingOptions,
    scored_clusters: Vec<Vec<Candidate<'_>>>,
) -> Result<Vec<DomainResult>> {
    // ANARCI's independent E-value always uses the complete database as Z.
    let e_value_search_space = database.profiles().len();
    let domain_count = scored_clusters.len();
    let mut domains = Vec::with_capacity(domain_count);
    for (domain_index, mut cluster) in scored_clusters.into_iter().enumerate() {
        for candidate in cluster
            .iter_mut()
            .take(options.alternative_hit_count.saturating_add(1))
        {
            let aligned = candidate
                .alignment_source
                .take()
                .and_then(|source| source.align(candidate.profile))
                .unwrap_or_else(|| {
                    realign_envelope(candidate.profile, &candidate.domain, &normalized.encoded)
                });
            candidate.query_start = aligned.start;
            candidate.query_end = aligned.end;
            candidate.aligned_domain = Some(aligned);
        }
        let best = cluster.remove(0);
        let aligned_domain = if domain_count == 1 {
            best.aligned_domain.unwrap_or_else(|| {
                realign_envelope(best.profile, &best.domain, &normalized.encoded)
            })
        } else {
            realign_domain(best.profile, &best.domain, &normalized.encoded)
        };
        let aligned_domain = if domain_count == 1 {
            recover_long_cdr3(database.profiles(), &aligned_domain, &normalized.encoded)
        } else {
            aligned_domain
        };
        let (numbering, padded_imgt_alignment, start, end) = number_imgt(
            normalized.text.as_bytes(),
            &aligned_domain,
            best.profile.chain_type(),
            domain_index == 0,
            domain_count == 1,
        );
        let alternative_hits = cluster
            .iter()
            .take(options.alternative_hit_count)
            .map(|candidate| candidate_hit(candidate, e_value_search_space))
            .collect();
        let germline = if options.assign_germline {
            assign_closest_germline(
                &aligned_domain,
                normalized.text.as_bytes(),
                best.profile.chain_type(),
                options.allowed_species.as_deref(),
            )?
        } else {
            None
        };
        domains.push(DomainResult {
            domain_index,
            receptor_type: best.profile.chain_type().receptor_type(),
            chain_type: best.profile.chain_type(),
            species: best.profile.species().to_owned(),
            start,
            end,
            bit_score: best.bit_score,
            e_value: independent_domain_e_value(best.profile, best.bit_score, e_value_search_space),
            bias: best.bias,
            query_start: best.query_start,
            query_end: best.query_end,
            numbering,
            padded_imgt_alignment,
            alternative_hits,
            germline,
        });
    }
    Ok(domains)
}

fn ranked_profile_indices(
    database: &ProfileDatabase<'_>,
    sequence: &[u8],
    plan: &SearchPlan<'_>,
) -> Vec<(f32, usize)> {
    let mut ranked_profiles: Vec<_> = database
        .profiles()
        .iter()
        .enumerate()
        .filter(|(profile_index, _)| plan.allows_profile(*profile_index))
        .filter(|(_, profile)| msv_filter_passes(profile, sequence))
        .map(|(profile_index, profile)| {
            let filter_score = ungapped_filter_score(profile, sequence);
            (profile_significance(profile, filter_score), profile_index)
        })
        .collect();
    ranked_profiles.sort_by(|left, right| {
        right.0.total_cmp(&left.0).then_with(|| {
            database.profiles()[left.1]
                .name()
                .cmp(database.profiles()[right.1].name())
        })
    });
    ranked_profiles
}

fn estimated_domain_count(sequence_length: usize) -> usize {
    ((sequence_length + 49) / 100).clamp(1, 7)
}

fn initial_profile_indices(
    database: &ProfileDatabase<'_>,
    sequence: &NormalizedSequence,
    plan: &SearchPlan<'_>,
) -> Vec<usize> {
    let ranked = ranked_profile_indices(database, &sequence.encoded, plan);
    let budget = plan.initial_trace_budget(sequence.encoded.len(), ranked.len());
    ranked
        .into_iter()
        .take(budget)
        .map(|(_, profile_index)| profile_index)
        .collect()
}

fn batch_eligible_profile_indices(
    database: &ProfileDatabase<'_>,
    records: &[PreparedRecord<'_>],
    plan: &SearchPlan<'_>,
) -> Vec<Vec<usize>> {
    let mut rescored = rescore_buckets(records);
    let mut workspace = SequenceBatchViterbiWorkspace::new(database.model_length());
    let mut task_scores = Vec::new();
    for task in profile_batch_tasks(database, records) {
        score_profile_batch_task(database, records, &task, &mut workspace, &mut task_scores);
        for (sequence_index, score, profile_index) in task_scores.drain(..) {
            rescored[sequence_index].push((score, profile_index));
        }
    }

    finish_batch_eligible_profile_indices(database, records, rescored, plan)
}

fn rescore_buckets(records: &[PreparedRecord<'_>]) -> Vec<Vec<(f32, usize)>> {
    records
        .iter()
        .map(|record| Vec::with_capacity(record.profiles.len()))
        .collect()
}

fn profile_batch_tasks(
    database: &ProfileDatabase<'_>,
    records: &[PreparedRecord<'_>],
) -> Vec<ProfileBatchTask> {
    #[cfg(target_family = "wasm")]
    let _ = database;
    let mut sequences_by_profile_and_length: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for (sequence_index, record) in records.iter().enumerate() {
        if estimated_domain_count(record.sequence.encoded.len()) != 1 {
            continue;
        }
        for &profile_index in &record.profiles {
            sequences_by_profile_and_length
                .entry((profile_index, record.sequence.encoded.len()))
                .or_default()
                .push(sequence_index);
        }
    }

    sequences_by_profile_and_length
        .into_iter()
        .flat_map(|((profile_index, sequence_length), sequence_indices)| {
            #[cfg(target_family = "wasm")]
            let _ = sequence_length;
            sequence_indices
                .chunks(BATCH_TASK_SEQUENCE_LIMIT)
                .map(move |chunk| ProfileBatchTask {
                    profile_index,
                    sequence_indices: chunk.to_vec(),
                    #[cfg(not(target_family = "wasm"))]
                    work: profile_batch_work(
                        &database.profiles()[profile_index],
                        sequence_length,
                        chunk.len(),
                    ),
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

fn score_profile_batch_task(
    database: &ProfileDatabase<'_>,
    records: &[PreparedRecord<'_>],
    task: &ProfileBatchTask,
    workspace: &mut SequenceBatchViterbiWorkspace,
    output: &mut Vec<(usize, f32, usize)>,
) {
    let profile = &database.profiles()[task.profile_index];
    let encoded: Vec<&[u8]> = task
        .sequence_indices
        .iter()
        .map(|index| records[*index].sequence.encoded.as_slice())
        .collect();
    output.extend(
        task.sequence_indices
            .iter()
            .copied()
            .zip(workspace.score(profile, &encoded).iter().copied())
            .map(|(sequence_index, score)| {
                (
                    sequence_index,
                    profile_significance(profile, score),
                    task.profile_index,
                )
            }),
    );
}

fn finish_batch_eligible_profile_indices(
    database: &ProfileDatabase<'_>,
    records: &[PreparedRecord<'_>],
    rescored: Vec<Vec<(f32, usize)>>,
    plan: &SearchPlan<'_>,
) -> Vec<Vec<usize>> {
    rescored
        .into_iter()
        .zip(records)
        .map(|(mut profiles, record)| {
            if estimated_domain_count(record.sequence.encoded.len()) > 1 {
                return record.profiles.clone();
            }
            profiles.sort_by(|left, right| {
                right.0.total_cmp(&left.0).then_with(|| {
                    database.profiles()[left.1]
                        .name()
                        .cmp(database.profiles()[right.1].name())
                })
            });
            let budget = plan.batch_trace_budget(record.sequence.encoded.len(), profiles.len());
            profiles
                .into_iter()
                .take(budget)
                .map(|(_, profile_index)| profile_index)
                .collect()
        })
        .collect()
}

fn score_defined_domain_candidate(candidate: &mut Candidate<'_>, sequence: &[u8]) {
    let definition = define_domain(candidate.profile, &candidate.domain, sequence);
    candidate.domain.start = definition.start;
    candidate.domain.end = definition.end;
    let (components, alignment) = domain_score_components_with_posterior(
        candidate.profile,
        &candidate.domain,
        sequence,
        definition.trace_null2_bias,
    );
    candidate.bit_score = components.bit_score;
    candidate.bias = components.null2_bias_bits;
    candidate.alignment_source = alignment;
}

fn validate_options(options: &NumberingOptions) -> Result<()> {
    if !options.min_bit_score.is_finite() {
        return Err(Error::sequence(
            ErrorCode::InvalidOptions,
            "minBitScore must be finite",
            None,
        ));
    }
    if options.allowed_chains.as_ref().is_some_and(Vec::is_empty) {
        return Err(Error::sequence(
            ErrorCode::InvalidOptions,
            "allowedChains cannot be empty",
            None,
        ));
    }
    if options.allowed_species.as_ref().is_some_and(Vec::is_empty) {
        return Err(Error::sequence(
            ErrorCode::InvalidOptions,
            "allowedSpecies cannot be empty",
            None,
        ));
    }
    Ok(())
}

fn profile_is_allowed(profile: &Profile<'_>, options: &NumberingOptions) -> bool {
    let chain_allowed = options
        .allowed_chains
        .as_ref()
        .is_none_or(|chains| chains.contains(&profile.chain_type()));
    let species_allowed = options.allowed_species.as_ref().is_none_or(|species| {
        species
            .iter()
            .any(|allowed| allowed.eq_ignore_ascii_case(profile.species()))
    });
    chain_allowed && species_allowed
}

fn cluster_candidates<'a>(mut candidates: Vec<Candidate<'a>>) -> Vec<Vec<Candidate<'a>>> {
    candidates.sort_by(|left, right| {
        left.domain
            .start
            .cmp(&right.domain.start)
            .then_with(|| left.domain.end.cmp(&right.domain.end))
            .then_with(|| candidate_order(left, right))
    });
    let mut clusters = connected_components(candidates, domains_overlap);

    for cluster in &mut clusters {
        cluster.sort_by(candidate_order);
        let mut seen = BTreeSet::new();
        cluster.retain(|candidate| seen.insert(candidate.profile.name()));
    }
    clusters.sort_by_key(|cluster| {
        cluster
            .iter()
            .map(|candidate| candidate.domain.start)
            .min()
            .unwrap_or(usize::MAX)
    });
    clusters
}

fn domains_overlap(left: &Candidate<'_>, right: &Candidate<'_>) -> bool {
    let overlap = left
        .domain
        .end
        .min(right.domain.end)
        .saturating_sub(left.domain.start.max(right.domain.start));
    let shorter = (left.domain.end - left.domain.start).min(right.domain.end - right.domain.start);
    overlap * 2 >= shorter
}

fn candidate_order(left: &Candidate<'_>, right: &Candidate<'_>) -> std::cmp::Ordering {
    right
        .ranking_score()
        .total_cmp(&left.ranking_score())
        .then_with(|| left.profile.name().cmp(right.profile.name()))
}

fn sort_multidomain_candidates(candidates: &mut [Candidate<'_>]) {
    candidates.sort_by(candidate_order);
    let mut cohort_start = 0;
    while cohort_start < candidates.len() {
        let anchor = candidates[cohort_start].ranking_score();
        let cohort_end = candidates[cohort_start + 1..]
            .iter()
            .position(|candidate| anchor - candidate.ranking_score() > 0.01)
            .map_or(candidates.len(), |offset| cohort_start + 1 + offset);
        candidates[cohort_start..cohort_end].sort_by(|left, right| {
            right
                .bit_score
                .total_cmp(&left.bit_score)
                .then_with(|| left.profile.name().cmp(right.profile.name()))
        });
        cohort_start = cohort_end;
    }
}

fn profile_significance(profile: &Profile<'_>, bit_score: f32) -> f32 {
    profile.forward_lambda() * (bit_score - profile.forward_tau())
}

fn independent_domain_e_value(profile: &Profile<'_>, bit_score: f32, search_space: usize) -> f64 {
    let score = f64::from(bit_score);
    let tau = f64::from(profile.forward_tau());
    let lambda = f64::from(profile.forward_lambda());
    let log_survival = if score < tau {
        0.0
    } else {
        -lambda * (score - tau)
    };
    search_space as f64 * log_survival.exp()
}

fn candidate_hit(candidate: &Candidate<'_>, e_value_search_space: usize) -> ProfileHit {
    ProfileHit {
        profile: candidate.profile.name().to_owned(),
        chain_type: candidate.profile.chain_type(),
        species: candidate.profile.species().to_owned(),
        bit_score: candidate.bit_score,
        e_value: independent_domain_e_value(
            candidate.profile,
            candidate.bit_score,
            e_value_search_space,
        ),
        bias: candidate.bias,
        query_start: candidate.query_start,
        query_end: candidate.query_end,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ChainType;

    const VH: &str = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";
    const VL: &str = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR";

    #[test]
    fn reference_vh_and_vl_sequences_are_detected_and_numbered() {
        let vh = number_sequence(VH, &NumberingOptions::default()).unwrap();
        assert_eq!(vh.domains.len(), 1);
        assert_eq!(vh.domains[0].chain_type, ChainType::H);
        assert_eq!(vh.domains[0].alternative_hits[0].profile, "human_H");
        assert!(vh.domains[0].e_value.is_finite());
        assert!(vh.domains[0].e_value > 0.0);
        assert!(vh.domains[0].bias >= 0.0);
        assert!(vh.domains[0].query_start < vh.domains[0].query_end);
        assert!(vh.domains[0].alternative_hits[0].e_value > 0.0);
        assert!(vh.domains[0].alternative_hits[0].bias >= 0.0);
        assert!(
            vh.domains[0].alternative_hits[0].query_start
                < vh.domains[0].alternative_hits[0].query_end
        );
        assert_eq!(vh.domains[0].numbering[0].position, 1);
        assert!(vh.domains[0].padded_imgt_alignment.len() >= 120);

        let vl = number_sequence(VL, &NumberingOptions::default()).unwrap();
        assert_eq!(vl.domains.len(), 1);
        assert!(matches!(
            vl.domains[0].chain_type,
            ChainType::K | ChainType::L
        ));
        assert_eq!(vl.domains[0].alternative_hits[0].profile, "rat_K");
    }

    #[test]
    fn unknown_residue_is_supported_and_preserved() {
        let unknown_index = 60;
        let reference = number_sequence(VH, &NumberingOptions::default()).unwrap();
        let reference_residue = reference.domains[0]
            .numbering
            .iter()
            .find(|residue| residue.sequence_index == unknown_index)
            .unwrap();

        let mut sequence = VH.as_bytes().to_vec();
        sequence[unknown_index] = b'X';
        let sequence = String::from_utf8(sequence).unwrap();
        let result = number_sequence(&sequence, &NumberingOptions::default()).unwrap();

        assert_eq!(result.normalized_sequence, sequence);
        assert_eq!(result.domains.len(), 1);
        let residue = result.domains[0]
            .numbering
            .iter()
            .find(|residue| residue.sequence_index == unknown_index)
            .unwrap();
        assert_eq!(residue.amino_acid, 'X');
        assert_eq!(residue.position, reference_residue.position);
        assert_eq!(residue.insertion_code, reference_residue.insertion_code);
        assert!(result.domains[0].padded_imgt_alignment.contains('X'));
    }

    #[test]
    fn increasing_alternative_count_preserves_existing_hits() {
        let with_three = number_sequence(VH, &NumberingOptions::default()).unwrap();
        let with_all = number_sequence(
            VH,
            &NumberingOptions {
                alternative_hit_count: embedded_profiles().unwrap().profiles().len() - 1,
                ..NumberingOptions::default()
            },
        )
        .unwrap();
        let three = &with_three.domains[0];
        let all = &with_all.domains[0];
        assert_eq!(three.chain_type, all.chain_type);
        assert_eq!(three.species, all.species);
        assert_eq!(
            (three.query_start, three.query_end),
            (all.query_start, all.query_end)
        );
        assert_eq!(three.alternative_hits, all.alternative_hits[..3]);
    }

    #[test]
    fn large_profile_batches_are_split_into_equal_length_tasks() {
        let inputs: Vec<_> = (0..BATCH_TASK_SEQUENCE_LIMIT + 1)
            .map(|index| SequenceInput {
                id: index.to_string(),
                sequence: VH.to_owned(),
            })
            .collect();
        let normalized = inputs
            .iter()
            .map(|input| normalize_sequence(&input.id, &input.sequence).unwrap())
            .collect();
        let profiles = vec![vec![0]; inputs.len()];
        let records = records_with_profiles(&inputs, normalized, profiles);
        let tasks = profile_batch_tasks(embedded_profiles().unwrap(), &records);

        assert_eq!(tasks.len(), 2);
        assert_eq!(tasks[0].sequence_indices.len(), BATCH_TASK_SEQUENCE_LIMIT);
        assert_eq!(tasks[1].sequence_indices.len(), 1);
        assert!(tasks.iter().all(|task| {
            let length = records[task.sequence_indices[0]].sequence.encoded.len();
            task.sequence_indices
                .iter()
                .all(|index| records[*index].sequence.encoded.len() == length)
        }));
    }

    #[test]
    fn allowed_species_is_strict() {
        let options = NumberingOptions {
            allowed_species: Some(vec!["cow".into()]),
            ..NumberingOptions::default()
        };
        let result = number_sequence(VH, &options).unwrap();
        assert!(result.domains.iter().all(|domain| domain.species == "cow"));
    }

    #[test]
    fn antibody_pair_validation_accepts_vh_vl_and_rejects_swapped_chains() {
        let valid = validate_antibody_pair(VH, VL, &PairValidationOptions::default()).unwrap();
        assert!(valid.ok);

        let swapped = validate_antibody_pair(VL, VH, &PairValidationOptions::default()).unwrap();
        assert!(!swapped.ok);
        assert!(swapped.errors.len() >= 2);
    }

    #[test]
    fn closest_germlines_match_the_pinned_reference() {
        let options = NumberingOptions {
            assign_germline: true,
            ..NumberingOptions::default()
        };
        let heavy = number_sequence(VH, &options).unwrap();
        let heavy_assignment = heavy.domains[0].germline.as_ref().unwrap();
        assert_eq!(heavy_assignment.species, "mouse");
        assert_eq!(heavy_assignment.v_gene.as_deref(), Some("IGHV14-4*02"));
        assert_eq!(heavy_assignment.j_gene.as_deref(), Some("IGHJ3*01"));
        assert!((heavy_assignment.v_identity.unwrap() - 0.948_979_6).abs() < 1e-6);

        let light = number_sequence(VL, &options).unwrap();
        let light_assignment = light.domains[0].germline.as_ref().unwrap();
        assert_eq!(light_assignment.species, "mouse");
        assert_eq!(light_assignment.v_gene.as_deref(), Some("IGKV6-13*01"));
        assert_eq!(light_assignment.j_gene.as_deref(), Some("IGKJ5*01"));
        assert_eq!(light_assignment.j_identity, Some(1.0));
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn parallel_batch_matches_serial_order_and_results() {
        let mut inputs: Vec<_> = (0..8)
            .flat_map(|index| {
                [
                    SequenceInput {
                        id: format!("heavy-{index}"),
                        sequence: VH.into(),
                    },
                    SequenceInput {
                        id: format!("light-{index}"),
                        sequence: VL.into(),
                    },
                ]
            })
            .collect();
        inputs[0].sequence.replace_range(60..61, "X");
        let options = NumberingOptions::default();
        let scalar: Vec<_> = inputs
            .iter()
            .map(|input| number_sequence_with_id(&input.id, &input.sequence, &options).unwrap())
            .collect();
        let serial = number_sequences(&inputs, &options).unwrap();
        assert_eq!(serial, scalar);
        for worker_count in [1, 2, 8] {
            let parallel = number_sequences_parallel(
                &inputs,
                &options,
                NonZeroUsize::new(worker_count).unwrap(),
            )
            .unwrap();
            assert_eq!(parallel, scalar);
        }

        let expanded_options = NumberingOptions {
            alternative_hit_count: 28,
            assign_germline: true,
            ..NumberingOptions::default()
        };
        let serial = number_sequences(&inputs, &expanded_options).unwrap();
        let parallel =
            number_sequences_parallel(&inputs, &expanded_options, NonZeroUsize::new(8).unwrap())
                .unwrap();
        assert_eq!(parallel, serial);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn native_batches_are_not_record_count_limited() {
        let inputs: Vec<_> = (0..1_001)
            .map(|index| SequenceInput {
                id: format!("short-{index}"),
                sequence: "A".into(),
            })
            .collect();
        let options = NumberingOptions::default();
        let serial = number_sequences(&inputs, &options).unwrap();
        let parallel =
            number_sequences_parallel(&inputs, &options, NonZeroUsize::new(4).unwrap()).unwrap();
        assert_eq!(serial, parallel);
        assert_eq!(parallel.len(), inputs.len());
        assert_eq!(parallel.last().unwrap().id, "short-1000");
    }

    #[test]
    fn sequence_major_batch_preserves_marginal_and_multidomain_results() {
        let corpus: serde_json::Value =
            serde_json::from_str(include_str!("../../../tests/golden/corpus_v2.json")).unwrap();
        let cases = corpus.as_array().unwrap();
        let selected: Vec<_> = ["cdr1len_vl_00", "scfv_okt3_vl_vh_mouse"]
            .into_iter()
            .cycle()
            .take(8)
            .enumerate()
            .map(|(index, id)| {
                let case = cases
                    .iter()
                    .find(|case| case["id"].as_str() == Some(id))
                    .unwrap();
                SequenceInput {
                    id: format!("{id}-{index}"),
                    sequence: case["seq"].as_str().unwrap().to_owned(),
                }
            })
            .collect();
        let options = NumberingOptions::default();
        let scalar: Vec<_> = selected
            .iter()
            .map(|input| number_sequence_with_id(&input.id, &input.sequence, &options).unwrap())
            .collect();
        let batch = number_sequences(&selected, &options).unwrap();
        assert_eq!(batch, scalar);
    }
}
