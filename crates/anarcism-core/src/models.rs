use crate::{ChainType, Error, Result};

const MAGIC: &[u8; 8] = b"ANRCPRF1";
const FORMAT_VERSION: u16 = 1;
const ALPHABET_SIZE: usize = 20;
const TRANSITION_COUNT: usize = 7;

pub fn embedded_profiles() -> Result<ProfileDatabase<'static>> {
    ProfileDatabase::from_bytes(include_bytes!("../../../assets/profiles.bin"))
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
            let forward_tau = reader.f32()?;
            let forward_lambda = reader.f32()?;
            if !forward_tau.is_finite() || !forward_lambda.is_finite() || forward_lambda <= 0.0 {
                return Err(Error::model("invalid profile E-value calibration"));
            }
            let consensus = reader.take(model_length)?;
            let local_entry = reader.take(checked_size(&[model_length, 2])?)?;
            let match_scores = reader.take(checked_size(&[model_length, ALPHABET_SIZE, 2])?)?;
            let transition_scores =
                reader.take(checked_size(&[model_length - 1, TRANSITION_COUNT, 2])?)?;
            profiles.push(Profile {
                name,
                species,
                chain_type,
                checksum,
                forward_tau,
                forward_lambda,
                consensus,
                local_entry,
                match_scores,
                transition_scores,
                scale: f32::from(scale),
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
    forward_tau: f32,
    forward_lambda: f32,
    consensus: &'a [u8],
    local_entry: &'a [u8],
    match_scores: &'a [u8],
    transition_scores: &'a [u8],
    scale: f32,
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

    pub const fn forward_tau(&self) -> f32 {
        self.forward_tau
    }

    pub const fn forward_lambda(&self) -> f32 {
        self.forward_lambda
    }

    pub const fn consensus(&self) -> &[u8] {
        self.consensus
    }

    pub fn local_entry_score(&self, model_position: usize) -> f32 {
        decode_score(self.local_entry, model_position - 1, self.scale)
    }

    pub fn match_score(&self, model_position: usize, residue_index: usize) -> f32 {
        decode_score(
            self.match_scores,
            (model_position - 1) * ALPHABET_SIZE + residue_index,
            self.scale,
        )
    }

    pub fn transition_score(&self, model_position: usize, transition: usize) -> f32 {
        decode_score(
            self.transition_scores,
            (model_position - 1) * TRANSITION_COUNT + transition,
            self.scale,
        )
    }
}

fn decode_score(bytes: &[u8], index: usize, scale: f32) -> f32 {
    let offset = index * 2;
    let quantized = i16::from_le_bytes([bytes[offset], bytes[offset + 1]]);
    if quantized == i16::MIN {
        f32::NEG_INFINITY
    } else {
        f32::from(quantized) / scale
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
}
