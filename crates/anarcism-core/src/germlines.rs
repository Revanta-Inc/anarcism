use crate::hmm::{RawDomain, TraceState};
use crate::{ChainType, Error, GermlineAssignment, Result};

const MAGIC: &[u8; 8] = b"ANRCGER1";
const FORMAT_VERSION: u16 = 1;
const ALIGNMENT_LENGTH: usize = 128;
const PACKED_LENGTH: usize = ALIGNMENT_LENGTH * 5 / 8;
const ALPHABET: &[u8; 21] = b"-ACDEFGHIKLMNPQRSTVWY";

pub fn embedded_germlines() -> Result<GermlineDatabase<'static>> {
    GermlineDatabase::from_bytes(include_bytes!("../../../assets/germlines.bin"))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Segment {
    V,
    J,
}

#[derive(Clone, Debug)]
pub struct GermlineDatabase<'a> {
    germlines: Vec<Germline<'a>>,
}

impl<'a> GermlineDatabase<'a> {
    pub fn from_bytes(bytes: &'a [u8]) -> Result<Self> {
        let mut reader = Reader::new(bytes);
        if reader.take(MAGIC.len())? != MAGIC {
            return Err(Error::model("invalid germline-data magic"));
        }
        if reader.u16()? != FORMAT_VERSION {
            return Err(Error::model("unsupported germline-data version"));
        }
        let species_count = usize::from(reader.u8()?);
        let germline_count = usize::from(reader.u16()?);
        if species_count == 0 || germline_count == 0 {
            return Err(Error::model("germline database is empty"));
        }
        if species_count > 32 || germline_count > 4_096 {
            return Err(Error::model(
                "germline database dimensions exceed format limits",
            ));
        }
        let mut species = Vec::with_capacity(species_count);
        for _ in 0..species_count {
            species.push(reader.short_string()?);
        }

        let mut germlines = Vec::with_capacity(germline_count);
        for _ in 0..germline_count {
            let segment = match reader.u8()? {
                b'V' => Segment::V,
                b'J' => Segment::J,
                _ => return Err(Error::model("invalid germline segment")),
            };
            let chain_type = ChainType::from_byte(reader.u8()?)
                .ok_or_else(|| Error::model("invalid germline chain type"))?;
            let species_index = usize::from(reader.u8()?);
            let species = species
                .get(species_index)
                .copied()
                .ok_or_else(|| Error::model("invalid germline species index"))?;
            let gene = reader.short_string()?;
            let sequence = reader.take(PACKED_LENGTH)?;
            for position in 0..ALIGNMENT_LENGTH {
                if decode_residue(sequence, position)? >= ALPHABET.len() {
                    return Err(Error::model("invalid packed germline residue"));
                }
            }
            germlines.push(Germline {
                segment,
                chain_type,
                species,
                gene,
                sequence,
            });
        }
        if !reader.remaining.is_empty() {
            return Err(Error::model("trailing bytes in germline database"));
        }
        Ok(Self { germlines })
    }

    pub fn len(&self) -> usize {
        self.germlines.len()
    }

    pub fn is_empty(&self) -> bool {
        self.germlines.is_empty()
    }
}

#[derive(Clone, Copy, Debug)]
struct Germline<'a> {
    segment: Segment,
    chain_type: ChainType,
    species: &'a str,
    gene: &'a str,
    sequence: &'a [u8],
}

pub(crate) fn assign_closest_germline(
    domain: &RawDomain,
    sequence: &[u8],
    chain_type: ChainType,
    allowed_species: Option<&[String]>,
) -> Result<Option<GermlineAssignment>> {
    let database = embedded_germlines()?;
    let state_sequence = state_sequence(domain, sequence);
    let best_v = best_match(
        &database,
        Segment::V,
        chain_type,
        allowed_species,
        &state_sequence,
    )?;
    let Some((v, v_identity)) = best_v else {
        return Ok(None);
    };
    let assigned_species = [v.species.to_owned()];
    let best_j = best_match(
        &database,
        Segment::J,
        chain_type,
        Some(&assigned_species),
        &state_sequence,
    )?;

    Ok(Some(GermlineAssignment {
        species: v.species.to_owned(),
        v_gene: Some(v.gene.to_owned()),
        v_identity: Some(v_identity),
        j_gene: best_j.map(|(germline, _)| germline.gene.to_owned()),
        j_identity: best_j.map(|(_, identity)| identity),
    }))
}

fn state_sequence(domain: &RawDomain, sequence: &[u8]) -> [u8; ALIGNMENT_LENGTH] {
    let mut result = [b'-'; ALIGNMENT_LENGTH];
    for step in &domain.steps {
        if step.state != TraceState::Match {
            continue;
        }
        let Some(sequence_index) = step.sequence_index else {
            continue;
        };
        let model_index = usize::from(step.model_position.saturating_sub(1));
        if model_index < ALIGNMENT_LENGTH && sequence_index < sequence.len() {
            result[model_index] = sequence[sequence_index].to_ascii_uppercase();
        }
    }
    result
}

fn best_match<'a>(
    database: &'a GermlineDatabase<'a>,
    segment: Segment,
    chain_type: ChainType,
    allowed_species: Option<&[String]>,
    state_sequence: &[u8; ALIGNMENT_LENGTH],
) -> Result<Option<(&'a Germline<'a>, f32)>> {
    let mut best = None;
    for germline in &database.germlines {
        if germline.segment != segment || germline.chain_type != chain_type {
            continue;
        }
        if allowed_species.is_some_and(|allowed| {
            !allowed
                .iter()
                .any(|species| species.eq_ignore_ascii_case(germline.species))
        }) {
            continue;
        }
        let identity = identity(state_sequence, germline.sequence)?;
        if best.is_none_or(|(_, best_identity)| identity > best_identity) {
            best = Some((germline, identity));
        }
    }
    Ok(best)
}

fn identity(state_sequence: &[u8; ALIGNMENT_LENGTH], packed: &[u8]) -> Result<f32> {
    let mut compared = 0_u16;
    let mut matched = 0_u16;
    for (position, observed) in state_sequence.iter().enumerate() {
        let residue_index = decode_residue(packed, position)?;
        let Some(expected) = ALPHABET.get(residue_index).copied() else {
            return Err(Error::model("invalid packed germline residue"));
        };
        if expected == b'-' {
            continue;
        }
        compared += 1;
        matched += u16::from(*observed == expected);
    }
    if compared == 0 {
        Ok(0.0)
    } else {
        Ok(f32::from(matched) / f32::from(compared))
    }
}

fn decode_residue(packed: &[u8], position: usize) -> Result<usize> {
    let bit_offset = position * 5;
    let byte_offset = bit_offset / 8;
    let shift = bit_offset % 8;
    let first = packed
        .get(byte_offset)
        .copied()
        .ok_or_else(|| Error::model("truncated packed germline sequence"))?;
    let mut value = u16::from(first) >> shift;
    if shift > 3 {
        let second = packed
            .get(byte_offset + 1)
            .copied()
            .ok_or_else(|| Error::model("truncated packed germline sequence"))?;
        value |= u16::from(second) << (8 - shift);
    }
    Ok(usize::from(value & 0x1f))
}

struct Reader<'a> {
    remaining: &'a [u8],
}

impl<'a> Reader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { remaining: bytes }
    }

    fn take(&mut self, length: usize) -> Result<&'a [u8]> {
        if length > self.remaining.len() {
            return Err(Error::model("truncated germline database"));
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

    fn short_string(&mut self) -> Result<&'a str> {
        let length = usize::from(self.u8()?);
        std::str::from_utf8(self.take(length)?)
            .map_err(|_| Error::model("germline database contains invalid UTF-8"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_database_contains_pinned_inventory() {
        let database = embedded_germlines().unwrap();
        assert_eq!(database.len(), 2_389);
    }

    #[test]
    fn decoder_rejects_truncated_and_corrupt_data() {
        assert!(GermlineDatabase::from_bytes(b"short").is_err());
        let mut bytes = include_bytes!("../../../assets/germlines.bin").to_vec();
        bytes.truncate(bytes.len() - 1);
        assert!(GermlineDatabase::from_bytes(&bytes).is_err());
    }
}
