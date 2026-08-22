use std::sync::OnceLock;

use crate::{ChainType, Error, Result};

const MAGIC: &[u8; 8] = b"ANRCPRF3";
const FORMAT_VERSION: u16 = 3;
const ALPHABET_SIZE: usize = 20;
const TRANSITION_COUNT: usize = 7;

/// Decodes the embedded profile database once per process and lends it out.
///
/// Numbering reads the database for every sequence, so the decode is cached
/// rather than repeated. The cached value borrows the embedded bytes, which
/// live for the whole program, so the returned reference is `'static`.
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
        if profile_count > 64 || model_length > 256 {
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
            let msv_match_costs = reader.take(checked_size(&[model_length, ALPHABET_SIZE])?)?;
            let local_entry = decode_scores(
                reader.take(checked_size(&[model_length, 3])?)?,
                inverse_scale,
            );
            let match_scores = decode_scores(
                reader.take(checked_size(&[model_length, ALPHABET_SIZE, 3])?)?,
                inverse_scale,
            );
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
    msv_match_costs: &'a [u8],
    forward_tau: f32,
    forward_lambda: f32,
    consensus: &'a [u8],
    local_entry: Box<[f32]>,
    match_scores: Box<[f32]>,
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
    pub fn msv_match_costs_for_residue(&self, residue_index: usize) -> &[u8] {
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
    pub fn local_entry_score(&self, model_position: usize) -> f32 {
        self.local_entry[model_position - 1]
    }

    #[inline]
    pub fn match_score(&self, model_position: usize, residue_index: usize) -> f32 {
        self.match_scores[(model_position - 1) * ALPHABET_SIZE + residue_index]
    }

    #[inline]
    pub fn transition_score(&self, model_position: usize, transition: usize) -> f32 {
        self.transition_scores[(model_position - 1) * TRANSITION_COUNT + transition]
    }

    #[inline]
    pub(crate) fn local_entry_scores(&self) -> &[f32] {
        &self.local_entry
    }

    #[inline]
    pub(crate) fn match_score_rows(&self) -> &[[f32; ALPHABET_SIZE]] {
        let (rows, remainder) = self.match_scores.as_chunks();
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
}
