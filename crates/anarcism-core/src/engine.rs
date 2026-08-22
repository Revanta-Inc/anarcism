use std::collections::BTreeSet;
#[cfg(not(target_family = "wasm"))]
use std::num::NonZeroUsize;

use crate::germlines::assign_closest_germline;
use crate::hmm::{
    RawDomain, define_domain, domain_score_components_with_hit_bounds,
    domain_score_components_with_null2, msv_filter_passes, realign_domain, realign_envelope,
    recover_long_cdr3, trace_emission_score, ungapped_filter_score, viterbi_domains,
};
use crate::numbering::number_imgt;
use crate::sequence::normalize_sequence;
use crate::{
    ChainType, DomainResult, Error, ErrorCode, NumberingOptions, PairValidationOptions,
    PairValidationResult, Profile, ProfileHit, Result, SequenceInput, SequenceResult,
    ValidationLimits, embedded_profiles,
};

#[derive(Debug)]
struct Candidate<'a> {
    profile: &'a Profile<'a>,
    domain: RawDomain,
    bit_score: f32,
    bias: f32,
    query_start: usize,
    query_end: usize,
    ranking_score: f32,
}

pub fn number_sequence(sequence: &str, options: &NumberingOptions) -> Result<SequenceResult> {
    number_sequence_with_id("sequence", sequence, options)
}

pub fn number_sequence_with_id(
    id: &str,
    sequence: &str,
    options: &NumberingOptions,
) -> Result<SequenceResult> {
    number_sequence_with_limits(id, sequence, options, ValidationLimits::default())
}

pub fn number_sequences(
    inputs: &[SequenceInput],
    options: &NumberingOptions,
) -> Result<Vec<SequenceResult>> {
    let limits = ValidationLimits::default();
    number_sequences_with_limits(inputs, options, limits)
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
    let limits = ValidationLimits::default();
    validate_batch_size(inputs, limits)?;
    let worker_count = worker_count.get().min(inputs.len());
    if worker_count <= 1 {
        return number_sequences_with_limits(inputs, options, limits);
    }

    std::thread::scope(|scope| {
        let workers: Vec<_> = (0..worker_count)
            .map(|worker_index| {
                let start = worker_index * inputs.len() / worker_count;
                let end = (worker_index + 1) * inputs.len() / worker_count;
                let chunk = &inputs[start..end];
                scope.spawn(move || {
                    chunk
                        .iter()
                        .map(|input| {
                            number_sequence_with_limits(&input.id, &input.sequence, options, limits)
                        })
                        .collect::<Result<Vec<_>>>()
                })
            })
            .collect();
        let mut results = Vec::with_capacity(inputs.len());
        for worker in workers {
            let shard = match worker.join() {
                Ok(shard) => shard?,
                Err(payload) => std::panic::resume_unwind(payload),
            };
            results.extend(shard);
        }
        Ok(results)
    })
}

pub fn number_fasta(fasta: &str, options: &NumberingOptions) -> Result<Vec<SequenceResult>> {
    let limits = ValidationLimits::default();
    let inputs = crate::fasta::parse_fasta(fasta, limits)?;
    number_sequences_with_limits(&inputs, options, limits)
}

pub fn validate_antibody_pair(
    vh: &str,
    vl: &str,
    options: &PairValidationOptions,
) -> Result<PairValidationResult> {
    validate_options(&options.numbering)?;
    let limits = ValidationLimits::default();
    if options.start_max >= options.end_min || options.end_min > limits.max_sequence_length {
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
        let result = number_sequence_with_limits("vh", vh, &heavy_options, limits)?;
        single_pair_domain(result.domains)
    };
    let light = if light_options
        .allowed_chains
        .as_ref()
        .is_some_and(Vec::is_empty)
    {
        None
    } else {
        let result = number_sequence_with_limits("vl", vl, &light_options, limits)?;
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

fn number_sequences_with_limits(
    inputs: &[SequenceInput],
    options: &NumberingOptions,
    limits: ValidationLimits,
) -> Result<Vec<SequenceResult>> {
    validate_batch_size(inputs, limits)?;
    inputs
        .iter()
        .map(|input| number_sequence_with_limits(&input.id, &input.sequence, options, limits))
        .collect()
}

fn validate_batch_size(inputs: &[SequenceInput], limits: ValidationLimits) -> Result<()> {
    if inputs.len() > limits.max_batch_size {
        return Err(Error::sequence(
            ErrorCode::BatchTooLarge,
            format!(
                "batch contains {} sequences; configured limit is {}",
                inputs.len(),
                limits.max_batch_size
            ),
            None,
        ));
    }
    Ok(())
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

fn number_sequence_with_limits(
    id: &str,
    sequence: &str,
    options: &NumberingOptions,
    limits: ValidationLimits,
) -> Result<SequenceResult> {
    validate_options(options)?;
    let normalized = normalize_sequence(id, sequence, limits)?;
    let database = embedded_profiles()?;
    // ANARCI runs hmmscan against the complete profile database. HMMER's
    // independent domain E-value therefore uses the full database size as Z,
    // even when caller-side chain/species filters restrict returned hits.
    let e_value_search_space = database.profiles().len();
    let mut ranked_profiles: Vec<(f32, &Profile<'_>)> = database
        .profiles()
        .iter()
        .filter(|profile| profile_is_allowed(profile, options))
        .filter(|profile| msv_filter_passes(profile, &normalized.encoded))
        .map(|profile| {
            let filter_score = ungapped_filter_score(profile, &normalized.encoded);
            (profile_significance(profile, filter_score), profile)
        })
        .collect();
    ranked_profiles.sort_by(|left, right| {
        right
            .0
            .total_cmp(&left.0)
            .then_with(|| left.1.name().cmp(right.1.name()))
    });
    let estimated_domains = ((normalized.encoded.len() + 49) / 100).clamp(1, 7);
    // Keep every species profile for the likely chain family in reach of the
    // exact Forward pass. Closely related H/K profiles can rank poorly under
    // an ungapped filter yet win decisively once their indel path is scored.
    let profiles_to_trace = options
        .alternative_hit_count
        .saturating_add(7)
        .max(10)
        .saturating_mul(estimated_domains)
        .min(ranked_profiles.len());
    let profiles_to_score_per_domain = options.alternative_hit_count.saturating_add(7).max(10);
    let eligible = ranked_profiles
        .into_iter()
        .take(profiles_to_trace)
        .map(|(_, profile)| profile);
    let mut candidates = Vec::new();
    for profile in eligible {
        let profile_domains = viterbi_domains(profile, &normalized.encoded);
        for domain in profile_domains.iter().cloned() {
            if domain.end <= domain.start {
                continue;
            }
            let query_start = domain.start;
            let query_end = domain.end;
            let bit_score = trace_emission_score(profile, &domain, &normalized.encoded);
            candidates.push(Candidate {
                profile,
                domain,
                bit_score,
                bias: 0.0,
                query_start,
                query_end,
                ranking_score: profile_significance(profile, bit_score),
            });
        }
    }

    let clusters = cluster_candidates(candidates);
    let mut raw_scored_clusters = Vec::with_capacity(clusters.len());
    for mut cluster in clusters {
        cluster.sort_by(candidate_order);
        // Viterbi trace scores eliminate the reserve profile before the
        // costlier Forward/Backward posterior pass. Every candidate that can
        // be returned (winner plus requested alternatives) remains exact.
        cluster.truncate(profiles_to_score_per_domain);
        for candidate in &mut cluster {
            score_defined_domain_candidate(
                candidate,
                &normalized.encoded,
                options.alternative_hit_count == 28,
            );
        }
        cluster.sort_by(candidate_order);
        raw_scored_clusters.push(cluster);
    }

    // Domain definition can expand several short Viterbi seeds onto the same
    // posterior envelope. Recluster those exact envelopes so one biological
    // domain is not reported once for every seed that led to it.
    let mut raw_scored_clusters =
        cluster_candidates(raw_scored_clusters.into_iter().flatten().collect());

    let is_multidomain_target = raw_scored_clusters.len() > 1;
    if is_multidomain_target {
        // Optimized HMMER and the scalar browser implementation can differ by
        // a few thousandths in calibrated significance. Within that narrow
        // numerical uncertainty band, prefer the higher domain bit score;
        // this mirrors the optimized ordering at multidomain boundaries.
        for cluster in &mut raw_scored_clusters {
            cluster.sort_by(multidomain_candidate_order);
        }
    }

    let mut scored_clusters = Vec::with_capacity(raw_scored_clusters.len());
    for mut cluster in raw_scored_clusters {
        cluster.retain(|candidate| candidate.bit_score >= options.min_bit_score);
        if cluster.is_empty() {
            continue;
        }
        if is_multidomain_target {
            cluster.sort_by(multidomain_candidate_order);
        } else {
            cluster.sort_by(candidate_order);
        }
        scored_clusters.push(cluster);
    }

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

    Ok(SequenceResult {
        id: id.to_owned(),
        normalized_sequence: normalized.text,
        domains,
        warnings: normalized.warnings,
    })
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
    candidate.ranking_score = profile_significance(candidate.profile, candidate.bit_score);
}

fn validate_options(options: &NumberingOptions) -> Result<()> {
    if !options.min_bit_score.is_finite() {
        return Err(Error::sequence(
            ErrorCode::InvalidOptions,
            "minBitScore must be finite",
            None,
        ));
    }
    if options.alternative_hit_count > 28 {
        return Err(Error::sequence(
            ErrorCode::InvalidOptions,
            "alternativeHitCount cannot exceed 28",
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
        .ranking_score
        .total_cmp(&left.ranking_score)
        .then_with(|| left.profile.name().cmp(right.profile.name()))
}

fn multidomain_candidate_order(left: &Candidate<'_>, right: &Candidate<'_>) -> std::cmp::Ordering {
    if (left.ranking_score - right.ranking_score).abs() <= 0.01 {
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
        let inputs = vec![
            SequenceInput {
                id: "heavy".into(),
                sequence: VH.into(),
            },
            SequenceInput {
                id: "light".into(),
                sequence: VL.into(),
            },
        ];
        let options = NumberingOptions::default();
        let serial = number_sequences(&inputs, &options).unwrap();
        for worker_count in [1, 2, 8] {
            let parallel = number_sequences_parallel(
                &inputs,
                &options,
                NonZeroUsize::new(worker_count).unwrap(),
            )
            .unwrap();
            assert_eq!(parallel, serial);
        }
    }
}
