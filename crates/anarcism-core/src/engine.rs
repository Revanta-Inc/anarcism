use std::collections::BTreeSet;

use crate::germlines::assign_closest_germline;
use crate::hmm::{
    RawDomain, forward_bit_score, trace_emission_score, trace_null2_bias_bits,
    ungapped_filter_score, viterbi_domains,
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
    inputs
        .iter()
        .map(|input| number_sequence_with_limits(&input.id, &input.sequence, options, limits))
        .collect()
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
    inputs
        .iter()
        .map(|input| number_sequence_with_limits(&input.id, &input.sequence, options, limits))
        .collect()
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
    let mut ranked_profiles: Vec<(f32, &Profile<'_>)> = database
        .profiles()
        .iter()
        .filter(|profile| profile_is_allowed(profile, options))
        .map(|profile| (ungapped_filter_score(profile, &normalized.encoded), profile))
        .collect();
    ranked_profiles.sort_by(|left, right| {
        right
            .0
            .total_cmp(&left.0)
            .then_with(|| left.1.name().cmp(right.1.name()))
    });
    let estimated_domains = ((normalized.encoded.len() + 49) / 100).clamp(1, 7);
    let profiles_to_trace = (options.alternative_hit_count + 5)
        .saturating_mul(estimated_domains)
        .min(ranked_profiles.len());
    let eligible = ranked_profiles
        .into_iter()
        .take(profiles_to_trace)
        .map(|(_, profile)| profile);
    let mut candidates = Vec::new();
    for profile in eligible {
        for domain in viterbi_domains(profile, &normalized.encoded) {
            if domain.end <= domain.start {
                continue;
            }
            let bit_score = trace_emission_score(profile, &domain, &normalized.encoded);
            candidates.push(Candidate {
                profile,
                domain,
                bit_score,
            });
        }
    }

    let clusters = cluster_candidates(candidates);
    let mut scored_clusters = Vec::with_capacity(clusters.len());
    for mut cluster in clusters {
        cluster.sort_by(candidate_order);
        for candidate in &mut cluster {
            candidate.bit_score = forward_bit_score(
                candidate.profile,
                &normalized.encoded[candidate.domain.start..candidate.domain.end],
            ) - trace_null2_bias_bits(
                candidate.profile,
                &candidate.domain,
                &normalized.encoded,
            );
        }
        cluster.retain(|candidate| candidate.bit_score >= options.min_bit_score);
        if cluster.is_empty() {
            continue;
        }
        cluster.sort_by(candidate_order);
        scored_clusters.push(cluster);
    }

    let domain_count = scored_clusters.len();
    let mut domains = Vec::with_capacity(domain_count);
    for (domain_index, cluster) in scored_clusters.into_iter().enumerate() {
        let best = &cluster[0];
        let (numbering, padded_imgt_alignment, start, end) = number_imgt(
            normalized.text.as_bytes(),
            &best.domain,
            best.profile.chain_type(),
            domain_index == 0,
            domain_count == 1,
        );
        let alternative_hits = cluster
            .iter()
            .skip(1)
            .take(options.alternative_hit_count)
            .map(candidate_hit)
            .collect();
        let germline = if options.assign_germline {
            assign_closest_germline(
                &best.domain,
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
            // Forward/null1 is accurate, but HMMER's final score and E-value also
            // apply null2 composition correction. Keep the public optional field
            // absent until that last correction is implemented and parity-tested.
            e_value: None,
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
        .bit_score
        .total_cmp(&left.bit_score)
        .then_with(|| left.profile.name().cmp(right.profile.name()))
}

fn candidate_hit(candidate: &Candidate<'_>) -> ProfileHit {
    ProfileHit {
        profile: candidate.profile.name().to_owned(),
        chain_type: candidate.profile.chain_type(),
        species: candidate.profile.species().to_owned(),
        bit_score: candidate.bit_score,
        e_value: None,
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
}
