use std::collections::{BTreeMap, BTreeSet};
#[cfg(not(target_family = "wasm"))]
use std::num::NonZeroUsize;
#[cfg(not(target_family = "wasm"))]
use std::sync::Mutex;

use crate::germlines::assign_closest_germline;
use crate::hmm::{
    RawDomain, SEQUENCE_BATCH_LANES, define_domain, domain_score_components_with_hit_bounds,
    domain_score_components_with_null2, msv_filter_passes, realign_domain, realign_envelope,
    recover_long_cdr3, trace_emission_score, ungapped_filter_score, viterbi_domains,
    viterbi_filter_bit_scores_batch,
};
use crate::numbering::number_imgt;
use crate::sequence::{MAX_SEQUENCE_LENGTH, NormalizedSequence, normalize_sequence};
use crate::{
    ChainType, DomainResult, Error, ErrorCode, NumberingOptions, PairValidationOptions,
    PairValidationResult, Profile, ProfileDatabase, ProfileHit, Result, SequenceInput,
    SequenceResult, embedded_profiles,
};

const MAX_ALTERNATIVE_HITS: usize = 28;
const EXACT_PROFILE_RESERVE: usize = 7;
const MIN_EXACT_PROFILES: usize = 10;
const BATCH_PROFILE_RESERVE: usize = 4;
const MIN_BATCH_PROFILES: usize = 7;

#[derive(Debug)]
struct Candidate<'a> {
    profile: &'a Profile<'a>,
    domain: RawDomain,
    bit_score: f32,
    bias: f32,
    query_start: usize,
    query_end: usize,
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
    validate_options(options)?;
    let normalized = normalize_sequence(id, sequence)?;
    let database = embedded_profiles()?;
    number_normalized_sequence(id, normalized, options, database)
}

pub fn number_sequences(
    inputs: &[SequenceInput],
    options: &NumberingOptions,
) -> Result<Vec<SequenceResult>> {
    validate_options(options)?;
    if inputs.is_empty() {
        return Ok(Vec::new());
    }

    let database = embedded_profiles()?;
    let normalized: Vec<_> = inputs
        .iter()
        .map(|input| normalize_sequence(&input.id, &input.sequence))
        .collect::<Result<_>>()?;
    let records = prepare_records(inputs, normalized, database, options);
    number_prepared_records(records, options, database).map(strip_result_indices)
}

/// Number a native batch concurrently while preserving input order.
///
/// The caller chooses the maximum worker count explicitly so applications
/// that already own a thread pool can avoid accidental oversubscription.
/// WebAssembly builds omit this API and retain the serial batch path.
#[cfg(not(target_family = "wasm"))]
pub fn number_sequences_parallel(
    inputs: &[SequenceInput],
    options: &NumberingOptions,
    worker_count: NonZeroUsize,
) -> Result<Vec<SequenceResult>> {
    validate_options(options)?;
    if inputs.is_empty() {
        return Ok(Vec::new());
    }

    let database = embedded_profiles()?;
    let normalized: Vec<_> = inputs
        .iter()
        .map(|input| normalize_sequence(&input.id, &input.sequence))
        .collect::<Result<_>>()?;
    let worker_count = worker_count.get().min(inputs.len());
    if worker_count <= 1 {
        let records = prepare_records(inputs, normalized, database, options);
        return number_prepared_records(records, options, database).map(strip_result_indices);
    }

    let initial_profiles =
        parallel_initial_profile_indices(database, &normalized, options, worker_count);
    let mut records = records_with_profiles(inputs, normalized, initial_profiles);
    if records.len() >= SEQUENCE_BATCH_LANES {
        let eligible_profiles =
            parallel_batch_eligible_profile_indices(database, &records, options, worker_count);
        replace_record_profiles(&mut records, eligible_profiles);
    }
    parallel_number_prepared_records(records, options, database, worker_count)
}

fn prepare_records<'a>(
    inputs: &'a [SequenceInput],
    normalized: Vec<NormalizedSequence>,
    database: &ProfileDatabase<'_>,
    options: &NumberingOptions,
) -> Vec<PreparedRecord<'a>> {
    let profiles = normalized
        .iter()
        .map(|sequence| initial_profile_indices(database, sequence, options))
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
    options: &NumberingOptions,
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
                                initial_profile_indices(database, &sequences[index], options),
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
    options: &NumberingOptions,
    worker_count: usize,
) -> Vec<Vec<usize>> {
    let weighted = profile_batch_tasks(database, records)
        .into_iter()
        .map(|task| (task.work, task.profile_index, task))
        .collect();
    let shards = distribute_by_weight(weighted, worker_count);

    let rescored = std::thread::scope(|scope| {
        let workers: Vec<_> = shards
            .into_iter()
            .filter(|shard| !shard.is_empty())
            .map(|shard| {
                scope.spawn(move || {
                    let mut scores = Vec::new();
                    for task in shard {
                        scores.extend(score_profile_batch_task(database, records, &task));
                    }
                    scores
                })
            })
            .collect();
        let mut rescored: Vec<Vec<(f32, usize)>> = (0..records.len()).map(|_| Vec::new()).collect();
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
    finish_batch_eligible_profile_indices(database, records, rescored, options)
}

fn profile_batch_work(
    profile: &Profile<'_>,
    records: &[PreparedRecord<'_>],
    sequence_indices: &[usize],
) -> usize {
    let mut counts = BTreeMap::new();
    for index in sequence_indices {
        *counts
            .entry(records[*index].sequence.encoded.len())
            .or_insert(0_usize) += 1;
    }
    counts.into_iter().fold(0_usize, |work, (length, count)| {
        work.saturating_add(
            length
                .saturating_mul(count.div_ceil(SEQUENCE_BATCH_LANES))
                .saturating_mul(profile.consensus().len()),
        )
    })
}

#[cfg(not(target_family = "wasm"))]
fn parallel_number_prepared_records(
    mut records: Vec<PreparedRecord<'_>>,
    options: &NumberingOptions,
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
                        results.push(number_prepared_record(record, options, database)?);
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
    validate_options(&options.numbering)?;
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
    options: &NumberingOptions,
    database: &ProfileDatabase<'_>,
) -> Result<Vec<(usize, SequenceResult)>> {
    if records.len() >= SEQUENCE_BATCH_LANES {
        let eligible_profiles = batch_eligible_profile_indices(database, &records, options);
        replace_record_profiles(&mut records, eligible_profiles);
    }
    records
        .into_iter()
        .map(|record| number_prepared_record(record, options, database))
        .collect()
}

fn number_prepared_record(
    record: PreparedRecord<'_>,
    options: &NumberingOptions,
    database: &ProfileDatabase<'_>,
) -> Result<(usize, SequenceResult)> {
    number_normalized_with_profiles(
        record.id,
        record.sequence,
        options,
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
    options: &NumberingOptions,
    database: &ProfileDatabase<'_>,
) -> Result<SequenceResult> {
    let eligible = initial_profile_indices(database, &normalized, options)
        .into_iter()
        .map(|profile_index| &database.profiles()[profile_index]);
    number_normalized_with_profiles(id, normalized, options, database, eligible)
}

fn number_normalized_with_profiles<'a>(
    id: &str,
    normalized: NormalizedSequence,
    options: &NumberingOptions,
    database: &'a ProfileDatabase<'a>,
    eligible: impl IntoIterator<Item = &'a Profile<'a>>,
) -> Result<SequenceResult> {
    let candidates = collect_candidates(eligible, &normalized.encoded);
    let scored_clusters = score_candidate_clusters(
        candidates,
        &normalized.encoded,
        exact_profile_budget(options),
        options.alternative_hit_count == MAX_ALTERNATIVE_HITS,
        options.min_bit_score,
    );
    let domains = render_domains(database, &normalized, options, scored_clusters)?;

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
            });
        }
    }
    candidates
}

fn score_candidate_clusters<'a>(
    candidates: Vec<Candidate<'a>>,
    sequence: &[u8],
    exact_profile_budget: usize,
    recover_hit_bounds: bool,
    min_bit_score: f32,
) -> Vec<Vec<Candidate<'a>>> {
    let exact_candidates = cluster_candidates(candidates)
        .into_iter()
        .flat_map(|mut cluster| {
            // Viterbi trace scores eliminate the reserve profile before the
            // costlier Forward/Backward posterior pass. Every candidate that can
            // be returned (winner plus requested alternatives) remains exact.
            cluster.truncate(exact_profile_budget);
            for candidate in &mut cluster {
                score_defined_domain_candidate(candidate, sequence, recover_hit_bounds);
            }
            cluster
        })
        .collect();

    // Domain definition can expand several short Viterbi seeds onto the same
    // posterior envelope. Recluster those exact envelopes so one biological
    // domain is not reported once for every seed that led to it.
    let mut clusters = cluster_candidates(exact_candidates);

    let is_multidomain_target = clusters.len() > 1;
    if is_multidomain_target {
        // Optimized HMMER and the scalar browser implementation can differ by
        // a few thousandths in calibrated significance. Within that narrow
        // numerical uncertainty band, prefer the higher domain bit score;
        // this mirrors the optimized ordering at multidomain boundaries.
        for cluster in &mut clusters {
            cluster.sort_by(multidomain_candidate_order);
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
    // ANARCI runs hmmscan against the complete profile database. HMMER's
    // independent domain E-value therefore uses the full database size as Z,
    // even when caller-side chain/species filters restrict returned hits.
    let e_value_search_space = database.profiles().len();
    let domain_count = scored_clusters.len();
    let mut domains = Vec::with_capacity(domain_count);
    for (domain_index, cluster) in scored_clusters.into_iter().enumerate() {
        let best = &cluster[0];
        let aligned_domain = if domain_count == 1 {
            // Domain definition has already selected the exact simple-region
            // or stochastic-cluster envelope. OA alignment must stay inside
            // it; ANARCI's later J-region recovery handles a genuinely long
            // CDR3 without changing the HMMER domain score.
            realign_envelope(best.profile, &best.domain, &normalized.encoded)
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
            best.profile.species(),
            domain_index == 0,
            domain_count == 1,
        );
        let alternative_hits = cluster
            .iter()
            .skip(1)
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
    options: &NumberingOptions,
) -> Vec<(f32, usize)> {
    let mut ranked_profiles: Vec<_> = database
        .profiles()
        .iter()
        .enumerate()
        .filter(|(_, profile)| profile_is_allowed(profile, options))
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

fn initial_trace_budget(
    sequence_length: usize,
    options: &NumberingOptions,
    available_profiles: usize,
) -> usize {
    // Keep every species profile for the likely chain family in reach of the
    // exact Forward pass. Closely related H/K profiles can rank poorly under
    // an ungapped filter yet win decisively once their indel path is scored.
    exact_profile_budget(options)
        .saturating_mul(estimated_domain_count(sequence_length))
        .min(available_profiles)
}

fn exact_profile_budget(options: &NumberingOptions) -> usize {
    options
        .alternative_hit_count
        .saturating_add(EXACT_PROFILE_RESERVE)
        .max(MIN_EXACT_PROFILES)
}

fn initial_profile_indices(
    database: &ProfileDatabase<'_>,
    sequence: &NormalizedSequence,
    options: &NumberingOptions,
) -> Vec<usize> {
    let ranked = ranked_profile_indices(database, &sequence.encoded, options);
    let budget = initial_trace_budget(sequence.encoded.len(), options, ranked.len());
    ranked
        .into_iter()
        .take(budget)
        .map(|(_, profile_index)| profile_index)
        .collect()
}

fn batch_trace_budget(
    sequence_length: usize,
    options: &NumberingOptions,
    available_profiles: usize,
) -> usize {
    options
        .alternative_hit_count
        .saturating_add(BATCH_PROFILE_RESERVE)
        .max(MIN_BATCH_PROFILES)
        .saturating_mul(estimated_domain_count(sequence_length))
        .min(available_profiles)
}

fn batch_eligible_profile_indices(
    database: &ProfileDatabase<'_>,
    records: &[PreparedRecord<'_>],
    options: &NumberingOptions,
) -> Vec<Vec<usize>> {
    let mut rescored: Vec<Vec<(f32, usize)>> = (0..records.len()).map(|_| Vec::new()).collect();
    for task in profile_batch_tasks(database, records) {
        for (sequence_index, score, profile_index) in
            score_profile_batch_task(database, records, &task)
        {
            rescored[sequence_index].push((score, profile_index));
        }
    }

    finish_batch_eligible_profile_indices(database, records, rescored, options)
}

fn profile_batch_tasks(
    database: &ProfileDatabase<'_>,
    records: &[PreparedRecord<'_>],
) -> Vec<ProfileBatchTask> {
    let mut sequences_by_profile = vec![Vec::new(); database.profiles().len()];
    for (sequence_index, record) in records.iter().enumerate() {
        if estimated_domain_count(record.sequence.encoded.len()) != 1 {
            continue;
        }
        for &profile_index in &record.profiles {
            sequences_by_profile[profile_index].push(sequence_index);
        }
    }

    sequences_by_profile
        .into_iter()
        .enumerate()
        .filter_map(|(profile_index, sequence_indices)| {
            if sequence_indices.is_empty() {
                return None;
            }
            let work = profile_batch_work(
                &database.profiles()[profile_index],
                records,
                &sequence_indices,
            );
            Some(ProfileBatchTask {
                profile_index,
                sequence_indices,
                work,
            })
        })
        .collect()
}

fn score_profile_batch_task(
    database: &ProfileDatabase<'_>,
    records: &[PreparedRecord<'_>],
    task: &ProfileBatchTask,
) -> Vec<(usize, f32, usize)> {
    let profile = &database.profiles()[task.profile_index];
    let encoded: Vec<&[u8]> = task
        .sequence_indices
        .iter()
        .map(|index| records[*index].sequence.encoded.as_slice())
        .collect();
    task.sequence_indices
        .iter()
        .copied()
        .zip(viterbi_filter_bit_scores_batch(profile, &encoded))
        .map(|(sequence_index, score)| {
            (
                sequence_index,
                profile_significance(profile, score),
                task.profile_index,
            )
        })
        .collect()
}

fn finish_batch_eligible_profile_indices(
    database: &ProfileDatabase<'_>,
    records: &[PreparedRecord<'_>],
    rescored: Vec<Vec<(f32, usize)>>,
    options: &NumberingOptions,
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
            let budget = batch_trace_budget(record.sequence.encoded.len(), options, profiles.len());
            profiles
                .into_iter()
                .take(budget)
                .map(|(_, profile_index)| profile_index)
                .collect()
        })
        .collect()
}

fn score_defined_domain_candidate(
    candidate: &mut Candidate<'_>,
    sequence: &[u8],
    recover_hit_bounds: bool,
) {
    let definition = define_domain(candidate.profile, &candidate.domain, sequence);
    candidate.domain.start = definition.start;
    candidate.domain.end = definition.end;
    let (components, hit_bounds) = if recover_hit_bounds {
        domain_score_components_with_hit_bounds(
            candidate.profile,
            &candidate.domain,
            sequence,
            definition.trace_null2_bias,
        )
    } else {
        (
            domain_score_components_with_null2(
                candidate.profile,
                &candidate.domain,
                sequence,
                definition.trace_null2_bias,
            ),
            None,
        )
    };
    candidate.bit_score = components.bit_score;
    candidate.bias = components.null2_bias_bits;
    if let Some((query_start, query_end)) = hit_bounds {
        candidate.query_start = query_start;
        candidate.query_end = query_end;
    }
}

fn validate_options(options: &NumberingOptions) -> Result<()> {
    if !options.min_bit_score.is_finite() {
        return Err(Error::sequence(
            ErrorCode::InvalidOptions,
            "minBitScore must be finite",
            None,
        ));
    }
    if options.alternative_hit_count > MAX_ALTERNATIVE_HITS {
        return Err(Error::sequence(
            ErrorCode::InvalidOptions,
            format!("alternativeHitCount cannot exceed {MAX_ALTERNATIVE_HITS}"),
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
    let mut clusters: Vec<Vec<Candidate<'a>>> = Vec::new();
    for candidate in candidates {
        if let Some(cluster) = clusters.iter_mut().find(|cluster| {
            cluster
                .iter()
                .any(|other| domains_overlap(&candidate, other))
        }) {
            cluster.push(candidate);
        } else {
            clusters.push(vec![candidate]);
        }
    }

    // A profile can occasionally produce duplicate HSP-like paths for one domain.
    // Retain only its best candidate in each physical-domain cluster.
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

fn multidomain_candidate_order(left: &Candidate<'_>, right: &Candidate<'_>) -> std::cmp::Ordering {
    if (left.ranking_score() - right.ranking_score()).abs() <= 0.01 {
        right
            .bit_score
            .total_cmp(&left.bit_score)
            .then_with(|| left.profile.name().cmp(right.profile.name()))
    } else {
        candidate_order(left, right)
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
    fn feasibility_sequences_are_detected_and_numbered() {
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
        let inputs: Vec<_> = (0..8)
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
