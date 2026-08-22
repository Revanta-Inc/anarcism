use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum ChainType {
    H,
    K,
    L,
    A,
    B,
    G,
    D,
}

impl ChainType {
    pub const ALL: [Self; 7] = [
        Self::H,
        Self::K,
        Self::L,
        Self::A,
        Self::B,
        Self::G,
        Self::D,
    ];

    pub const fn receptor_type(self) -> ReceptorType {
        match self {
            Self::H | Self::K | Self::L => ReceptorType::IG,
            Self::A | Self::B | Self::G | Self::D => ReceptorType::TR,
        }
    }

    pub(crate) fn from_byte(value: u8) -> Option<Self> {
        Some(match value {
            b'H' => Self::H,
            b'K' => Self::K,
            b'L' => Self::L,
            b'A' => Self::A,
            b'B' => Self::B,
            b'G' => Self::G,
            b'D' => Self::D,
            _ => return None,
        })
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum ReceptorType {
    IG,
    TR,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum Region {
    FR1,
    CDR1,
    FR2,
    CDR2,
    FR3,
    CDR3,
    FR4,
}

impl Region {
    pub const fn for_imgt_position(position: u16) -> Self {
        match position {
            1..=26 => Self::FR1,
            27..=38 => Self::CDR1,
            39..=55 => Self::FR2,
            56..=65 => Self::CDR2,
            66..=104 => Self::FR3,
            105..=117 => Self::CDR3,
            _ => Self::FR4,
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(default, rename_all = "camelCase")]
pub struct NumberingOptions {
    pub allowed_chains: Option<Vec<ChainType>>,
    pub allowed_species: Option<Vec<String>>,
    pub min_bit_score: f32,
    pub alternative_hit_count: usize,
    pub assign_germline: bool,
}

impl Default for NumberingOptions {
    fn default() -> Self {
        Self {
            allowed_chains: None,
            allowed_species: None,
            min_bit_score: 80.0,
            alternative_hit_count: 3,
            assign_germline: false,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct NumberedResidue {
    pub sequence_index: usize,
    pub amino_acid: char,
    pub position: u16,
    pub insertion_code: String,
    pub region: Region,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ProfileHit {
    pub profile: String,
    pub chain_type: ChainType,
    pub species: String,
    pub bit_score: f32,
    pub e_value: f64,
    pub bias: f32,
    pub query_start: usize,
    pub query_end: usize,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct GermlineAssignment {
    pub species: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub v_gene: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub v_identity: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub j_gene: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub j_identity: Option<f32>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct DomainResult {
    pub domain_index: usize,
    pub receptor_type: ReceptorType,
    pub chain_type: ChainType,
    pub species: String,
    pub start: usize,
    pub end: usize,
    pub bit_score: f32,
    pub e_value: f64,
    pub bias: f32,
    pub query_start: usize,
    pub query_end: usize,
    pub numbering: Vec<NumberedResidue>,
    pub padded_imgt_alignment: String,
    pub alternative_hits: Vec<ProfileHit>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub germline: Option<GermlineAssignment>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct SequenceInput {
    pub id: String,
    pub sequence: String,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SequenceResult {
    pub id: String,
    pub normalized_sequence: String,
    pub domains: Vec<DomainResult>,
    pub warnings: Vec<String>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(default, rename_all = "camelCase")]
pub struct ValidationLimits {
    pub max_sequence_length: usize,
    pub max_batch_size: usize,
    pub max_fasta_bytes: usize,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(default, rename_all = "camelCase")]
pub struct PairValidationOptions {
    #[serde(flatten)]
    pub numbering: NumberingOptions,
    pub start_max: usize,
    pub end_min: usize,
}

impl Default for PairValidationOptions {
    fn default() -> Self {
        Self {
            numbering: NumberingOptions::default(),
            start_max: 10,
            end_min: 100,
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct PairValidationResult {
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vh: Option<DomainResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vl: Option<DomainResult>,
    pub errors: Vec<String>,
}

impl Default for ValidationLimits {
    fn default() -> Self {
        Self {
            max_sequence_length: 10_000,
            max_batch_size: 1_000,
            max_fasta_bytes: 10 * 1024 * 1024,
        }
    }
}
