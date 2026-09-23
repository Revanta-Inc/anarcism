use std::sync::OnceLock;

use crate::sequence::{CANONICAL_RESIDUE_COUNT, RESIDUE_CODE_COUNT, UNKNOWN_RESIDUE_INDEX};
use crate::{ChainType, Error, Result};

const MAGIC: &[u8; 8] = b"ANRCPRF3";
const FORMAT_VERSION: u16 = 3;
const TRANSITION_COUNT: usize = 7;
const MSV_SCORE_SCALE: f32 = 3.0 / std::f32::consts::LN_2;
// HMMER's default amino-acid background, in AMINO_ALPHABET order.
const BACKGROUND: [f32; CANONICAL_RESIDUE_COUNT] = [
    0.078_794_5,
    0.015_160_0,
    0.053_522_2,
    0.066_829_8,
    0.039_706_2,
    0.069_507_1,
    0.022_919_8,
    0.059_009_2,
    0.059_442_2,
    0.096_372_8,
    0.023_771_8,
    0.041_438_6,
    0.048_290_4,
    0.039_563_9,
    0.054_097_8,
    0.068_336_4,
    0.054_068_7,
    0.067_341_7,
    0.011_413_5,
    0.030_413_3,
];
pub(crate) const MAX_MODEL_LENGTH: usize = 256;

/// Return the process-wide decoded profile database.
pub fn embedded_profiles() -> Result<&'static ProfileDatabase<'static>> {
    static DATABASE: OnceLock<Result<ProfileDatabase<'static>>> = OnceLock::new();
    DATABASE
        .get_or_init(|| ProfileDatabase::from_bytes(include_bytes!("../../../assets/profiles.bin")))
        .as_ref()
        .map_err(Clone::clone)
}

#[derive(Clone, Debug)]
pub struct ProfileDatabase<'a> {
    profiles: Vec<Profile<'a>>,
    model_length: usize,
}

impl<'a> ProfileDatabase<'a> {
    pub fn from_bytes(bytes: &'a [u8]) -> Result<Self> {
        let mut reader = Reader::new(bytes);
        if reader.take(MAGIC.len())? != MAGIC {
            return Err(Error::model("invalid profile-data magic"));
        }
        if reader.u16()? != FORMAT_VERSION {
            return Err(Error::model("unsupported profile-data version"));
        }
        let scale = reader.u16()?;
        if scale == 0 {
            return Err(Error::model("profile-data score scale is zero"));
        }
        let profile_count = usize::from(reader.u16()?);
        let model_length = usize::from(reader.u16()?);
        if profile_count == 0 || model_length == 0 {
            return Err(Error::model("profile database is empty"));
        }
        if profile_count > 64 || model_length > MAX_MODEL_LENGTH {
            return Err(Error::model(
                "profile database dimensions exceed format limits",
            ));
        }

        let inverse_scale = f32::from(scale).recip();
        let mut profiles = Vec::with_capacity(profile_count);
        for _ in 0..profile_count {
            let name = reader.short_string()?;
            let species = reader.short_string()?;
            let chain_type = ChainType::from_byte(reader.u8()?)
                .ok_or_else(|| Error::model("invalid profile chain type"))?;
            let encoded_receptor = reader.u8()?;
            if encoded_receptor != u8::from(chain_type.receptor_type() == crate::ReceptorType::TR) {
                return Err(Error::model("profile receptor type does not match chain"));
            }
            let checksum = reader.u32()?;
            let msv_mu = reader.f32()?;
            let msv_lambda = reader.f32()?;
            let forward_tau = reader.f32()?;
            let forward_lambda = reader.f32()?;
            if !msv_mu.is_finite() || !msv_lambda.is_finite() || msv_lambda <= 0.0 {
                return Err(Error::model("invalid profile MSV calibration"));
            }
            if !forward_tau.is_finite() || !forward_lambda.is_finite() || forward_lambda <= 0.0 {
                return Err(Error::model("invalid profile E-value calibration"));
            }
            let consensus = reader.take(model_length)?;
            let msv_bias = reader.u8()?;
            let canonical_msv_match_costs =
                reader.take(checked_size(&[model_length, CANONICAL_RESIDUE_COUNT])?)?;
            let local_entry = decode_scores(
                reader.take(checked_size(&[model_length, 3])?)?,
                inverse_scale,
            );
            let canonical_match_scores = decode_scores(
                reader.take(checked_size(&[model_length, CANONICAL_RESIDUE_COUNT, 3])?)?,
                inverse_scale,
            );
            let match_scores = add_unknown_match_scores(&canonical_match_scores);
            let msv_match_costs = add_unknown_msv_match_costs(
                canonical_msv_match_costs,
                &match_scores,
                model_length,
                msv_bias,
            );
            let match_odds = match_scores.iter().map(|score| score.exp()).collect();
            let transition_scores = decode_scores(
                reader.take(checked_size(&[model_length - 1, TRANSITION_COUNT, 3])?)?,
                inverse_scale,
            );
            profiles.push(Profile {
                name,
                species,
                chain_type,
                checksum,
                msv_mu,
                msv_lambda,
                msv_bias,
                msv_match_costs,
                forward_tau,
                forward_lambda,
                consensus,
                local_entry,
                match_scores,
                match_odds,
                transition_scores,
            });
        }
        if !reader.remaining().is_empty() {
            return Err(Error::model("trailing bytes in profile database"));
        }
        Ok(Self {
            profiles,
            model_length,
        })
    }

    pub fn profiles(&self) -> &[Profile<'a>] {
        &self.profiles
    }

    pub const fn model_length(&self) -> usize {
        self.model_length
    }
}

fn checked_size(factors: &[usize]) -> Result<usize> {
    factors.iter().try_fold(1_usize, |size, factor| {
        size.checked_mul(*factor)
            .ok_or_else(|| Error::model("profile-data field size overflows"))
    })
}

#[derive(Clone, Debug)]
pub struct Profile<'a> {
    name: &'a str,
    species: &'a str,
    chain_type: ChainType,
    checksum: u32,
    msv_mu: f32,
    msv_lambda: f32,
    msv_bias: u8,
    msv_match_costs: Box<[u8]>,
    forward_tau: f32,
    forward_lambda: f32,
    consensus: &'a [u8],
    local_entry: Box<[f32]>,
    match_scores: Box<[f32]>,
    match_odds: Box<[f32]>,
    transition_scores: Box<[f32]>,
}

impl Profile<'_> {
    pub const fn name(&self) -> &str {
        self.name
    }

    pub const fn species(&self) -> &str {
        self.species
    }

    pub const fn chain_type(&self) -> ChainType {
        self.chain_type
    }

    pub const fn checksum(&self) -> u32 {
        self.checksum
    }

    pub const fn msv_mu(&self) -> f32 {
        self.msv_mu
    }

    pub const fn msv_lambda(&self) -> f32 {
        self.msv_lambda
    }

    pub const fn msv_bias(&self) -> u8 {
        self.msv_bias
    }

    #[inline]
    pub(crate) fn msv_match_costs_for_residue(&self, residue_index: usize) -> &[u8] {
        let start = residue_index * self.consensus.len();
        &self.msv_match_costs[start..start + self.consensus.len()]
    }

    pub const fn forward_tau(&self) -> f32 {
        self.forward_tau
    }

    pub const fn forward_lambda(&self) -> f32 {
        self.forward_lambda
    }

    pub const fn consensus(&self) -> &[u8] {
        self.consensus
    }

    #[inline]
    pub(crate) fn local_entry_score(&self, model_position: usize) -> f32 {
        self.local_entry[model_position - 1]
    }

    #[inline]
    pub(crate) fn match_score(&self, model_position: usize, residue_index: usize) -> f32 {
        self.match_scores[(model_position - 1) * RESIDUE_CODE_COUNT + residue_index]
    }

    #[inline]
    pub(crate) fn transition_score(&self, model_position: usize, transition: usize) -> f32 {
        self.transition_scores[(model_position - 1) * TRANSITION_COUNT + transition]
    }

    #[inline]
    pub(crate) fn local_entry_scores(&self) -> &[f32] {
        &self.local_entry
    }

    #[inline]
    pub(crate) fn match_score_rows(&self) -> &[[f32; RESIDUE_CODE_COUNT]] {
        let (rows, remainder) = self.match_scores.as_chunks();
        debug_assert!(remainder.is_empty());
        rows
    }

    #[inline]
    pub(crate) fn match_odds_rows(&self) -> &[[f32; RESIDUE_CODE_COUNT]] {
        let (rows, remainder) = self.match_odds.as_chunks();
        debug_assert!(remainder.is_empty());
        rows
    }

    #[inline]
    pub(crate) fn transition_score_rows(&self) -> &[[f32; TRANSITION_COUNT]] {
        let (rows, remainder) = self.transition_scores.as_chunks();
        debug_assert!(remainder.is_empty());
        rows
    }
}

fn add_unknown_match_scores(canonical_scores: &[f32]) -> Box<[f32]> {
    debug_assert_eq!(canonical_scores.len() % CANONICAL_RESIDUE_COUNT, 0);
    let mut scores =
        Vec::with_capacity(canonical_scores.len() / CANONICAL_RESIDUE_COUNT * RESIDUE_CODE_COUNT);
    for canonical_row in canonical_scores.chunks_exact(CANONICAL_RESIDUE_COUNT) {
        scores.extend_from_slice(canonical_row);
        scores.push(unknown_match_score(canonical_row));
    }
    scores.into_boxed_slice()
}

// HMMER configures the all-degenerate X symbol as the background-weighted
// expected score over the canonical amino acids.
fn unknown_match_score(canonical_scores: &[f32]) -> f32 {
    debug_assert_eq!(canonical_scores.len(), CANONICAL_RESIDUE_COUNT);
    let mut weighted_score = 0.0_f32;
    let mut total_probability = 0.0_f32;
    for (score, probability) in canonical_scores.iter().zip(BACKGROUND) {
        weighted_score += score * probability;
        total_probability += probability;
    }
    weighted_score / total_probability
}

fn add_unknown_msv_match_costs(
    canonical_costs: &[u8],
    match_scores: &[f32],
    model_length: usize,
    bias: u8,
) -> Box<[u8]> {
    debug_assert_eq!(
        canonical_costs.len(),
        model_length * CANONICAL_RESIDUE_COUNT
    );
    debug_assert_eq!(match_scores.len(), model_length * RESIDUE_CODE_COUNT);
    let mut costs = Vec::with_capacity(model_length * RESIDUE_CODE_COUNT);
    costs.extend_from_slice(canonical_costs);
    costs.extend(
        match_scores
            .chunks_exact(RESIDUE_CODE_COUNT)
            .map(|row| msv_biased_byte_cost(row[usize::from(UNKNOWN_RESIDUE_INDEX)], bias)),
    );
    costs.into_boxed_slice()
}

fn msv_biased_byte_cost(score: f32, bias: u8) -> u8 {
    if !score.is_finite() {
        return u8::MAX;
    }
    let cost = -(MSV_SCORE_SCALE * score).round() as i32 + i32::from(bias);
    cost.clamp(0, i32::from(u8::MAX)) as u8
}

fn decode_scores(bytes: &[u8], inverse_scale: f32) -> Box<[f32]> {
    bytes
        .chunks_exact(3)
        .map(|bytes| decode_score(bytes, inverse_scale))
        .collect()
}

fn decode_score(bytes: &[u8], inverse_scale: f32) -> f32 {
    debug_assert_eq!(bytes.len(), 3);
    let encoded = u32::from(bytes[0]) | (u32::from(bytes[1]) << 8) | (u32::from(bytes[2]) << 16);
    let quantized = if encoded & 0x80_0000 == 0 {
        encoded as i32
    } else {
        (encoded | 0xff00_0000) as i32
    };
    if quantized == -(1 << 23) {
        f32::NEG_INFINITY
    } else {
        quantized as f32 * inverse_scale
    }
}

struct Reader<'a> {
    remaining: &'a [u8],
}

impl<'a> Reader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { remaining: bytes }
    }

    const fn remaining(&self) -> &'a [u8] {
        self.remaining
    }

    fn take(&mut self, length: usize) -> Result<&'a [u8]> {
        if length > self.remaining.len() {
            return Err(Error::model("truncated profile database"));
        }
        let (value, remaining) = self.remaining.split_at(length);
        self.remaining = remaining;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16> {
        let bytes = self.take(2)?;
        Ok(u16::from_le_bytes([bytes[0], bytes[1]]))
    }

    fn u32(&mut self) -> Result<u32> {
        let bytes = self.take(4)?;
        Ok(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }

    fn f32(&mut self) -> Result<f32> {
        Ok(f32::from_bits(self.u32()?))
    }

    fn short_string(&mut self) -> Result<&'a str> {
        let length = usize::from(self.u8()?);
        std::str::from_utf8(self.take(length)?)
            .map_err(|_| Error::model("profile database contains invalid UTF-8"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_database_contains_pinned_profile_inventory() {
        let database = embedded_profiles().unwrap();
        assert_eq!(database.model_length(), 128);
        assert_eq!(database.profiles().len(), 29);
        assert_eq!(database.profiles()[0].name(), "human_H");
        assert_eq!(database.profiles()[28].name(), "mouse_D");
    }

    #[test]
    fn embedded_database_is_decoded_once_and_shared() {
        let first = embedded_profiles().unwrap();
        let second = embedded_profiles().unwrap();
        assert!(std::ptr::eq(first, second));
    }

    #[test]
    fn signed_24_bit_scores_decode_without_losing_the_sentinel() {
        let inverse_scale = 32_768.0_f32.recip();
        assert_eq!(decode_score(&[0, 0, 0], inverse_scale), 0.0);
        assert_eq!(
            decode_score(&[0xff, 0xff, 0xff], inverse_scale),
            -inverse_scale
        );
        assert_eq!(
            decode_score(&[0xff, 0xff, 0x7f], inverse_scale),
            ((1 << 23) - 1) as f32 * inverse_scale
        );
        assert_eq!(
            decode_score(&[0, 0, 0x80], inverse_scale),
            f32::NEG_INFINITY
        );
    }

    #[test]
    fn embedded_profiles_expand_x_with_hmmer_degenerate_scores() {
        let profile = &embedded_profiles().unwrap().profiles()[0];
        let unknown = usize::from(UNKNOWN_RESIDUE_INDEX);
        for model_position in [1, 64, 128] {
            let row = &profile.match_score_rows()[model_position - 1];
            let expected = unknown_match_score(&row[..CANONICAL_RESIDUE_COUNT]);
            assert_eq!(row[unknown], expected);
            assert_eq!(profile.match_score(model_position, unknown), expected);
            assert_eq!(
                profile.msv_match_costs_for_residue(unknown)[model_position - 1],
                msv_biased_byte_cost(expected, profile.msv_bias())
            );
        }
    }
}
