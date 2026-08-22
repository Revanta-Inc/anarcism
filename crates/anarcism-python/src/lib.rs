//! PyO3 bindings for the anarcism engine.
//!
//! The Python surface mirrors the browser package but uses snake_case names and
//! keyword arguments. Enumerations cross as plain strings so a result converts
//! to JSON without a custom encoder and matches the JavaScript contract.
//!
//! Every analysis call releases the GIL. `anarcism-core` is pure computation
//! over embedded models with no Python state, so releasing lets a thread pool
//! use more than one core.

use std::num::NonZeroUsize;

use anarcism_core::{
    ChainType, DomainResult as CoreDomain, Error, GermlineAssignment as CoreGermline,
    NumberedResidue as CoreResidue, NumberingOptions, PairValidationOptions,
    PairValidationResult as CorePair, ProfileHit as CoreHit, SequenceInput,
    SequenceResult as CoreSequence, embedded_profiles, number_fasta as core_number_fasta,
    number_sequence_with_id, number_sequences as core_number_sequences,
    number_sequences_parallel as core_number_sequences_parallel,
    validate_antibody_pair as core_validate_pair,
};
use pyo3::create_exception;
use pyo3::exceptions::{PyException, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

create_exception!(
    _anarcism,
    AnarcismError,
    PyException,
    "Raised when anarcism rejects an input or cannot read its embedded models.\n\n\
     Carries a stable `code` string and an optional `input_id` naming the record\n\
     that failed."
);

fn parse_chain(value: &str) -> PyResult<ChainType> {
    match value.to_ascii_uppercase().as_str() {
        "H" => Ok(ChainType::H),
        "K" => Ok(ChainType::K),
        "L" => Ok(ChainType::L),
        "A" => Ok(ChainType::A),
        "B" => Ok(ChainType::B),
        "G" => Ok(ChainType::G),
        "D" => Ok(ChainType::D),
        _ => Err(PyValueError::new_err(format!(
            "unknown chain type {value:?}; expected one of H, K, L, A, B, G, D"
        ))),
    }
}

/// Converts an engine error into `AnarcismError`, preserving the machine-readable code.
fn to_py_error(py: Python<'_>, error: &Error) -> PyErr {
    let raised = AnarcismError::new_err(error.to_string());
    let value = raised.value(py);
    // A freshly built exception instance always accepts attribute assignment.
    let _ = value.setattr("code", error.code.as_str());
    let _ = value.setattr("input_id", error.input_id.clone());
    raised
}

fn numbering_options(
    allowed_chains: Option<Vec<String>>,
    allowed_species: Option<Vec<String>>,
    min_bit_score: Option<f32>,
    alternative_hit_count: Option<usize>,
    assign_germline: bool,
) -> PyResult<NumberingOptions> {
    let defaults = NumberingOptions::default();
    let chains = match allowed_chains {
        None => None,
        Some(values) => Some(
            values
                .iter()
                .map(|value| parse_chain(value))
                .collect::<PyResult<Vec<_>>>()?,
        ),
    };
    Ok(NumberingOptions {
        allowed_chains: chains,
        allowed_species,
        min_bit_score: min_bit_score.unwrap_or(defaults.min_bit_score),
        alternative_hit_count: alternative_hit_count.unwrap_or(defaults.alternative_hit_count),
        assign_germline,
    })
}

#[pyclass(frozen, get_all, module = "anarcism._anarcism")]
#[derive(Debug)]
pub struct NumberedResidue {
    sequence_index: usize,
    amino_acid: String,
    position: u16,
    insertion_code: String,
    region: String,
}

#[pymethods]
impl NumberedResidue {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("sequence_index", self.sequence_index)?;
        dict.set_item("amino_acid", &self.amino_acid)?;
        dict.set_item("position", self.position)?;
        dict.set_item("insertion_code", &self.insertion_code)?;
        dict.set_item("region", &self.region)?;
        Ok(dict)
    }

    fn __repr__(&self) -> String {
        format!(
            "NumberedResidue({}{} {} {})",
            self.position, self.insertion_code, self.amino_acid, self.region
        )
    }
}

#[pyclass(frozen, get_all, module = "anarcism._anarcism")]
#[derive(Debug)]
pub struct ProfileHit {
    profile: String,
    chain_type: String,
    species: String,
    bit_score: f32,
    e_value: f64,
    bias: f32,
    query_start: usize,
    query_end: usize,
}

#[pymethods]
impl ProfileHit {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("profile", &self.profile)?;
        dict.set_item("chain_type", &self.chain_type)?;
        dict.set_item("species", &self.species)?;
        dict.set_item("bit_score", self.bit_score)?;
        dict.set_item("e_value", self.e_value)?;
        dict.set_item("bias", self.bias)?;
        dict.set_item("query_start", self.query_start)?;
        dict.set_item("query_end", self.query_end)?;
        Ok(dict)
    }

    fn __repr__(&self) -> String {
        format!("ProfileHit({} {:.4})", self.profile, self.bit_score)
    }
}

#[pyclass(frozen, get_all, module = "anarcism._anarcism")]
#[derive(Debug)]
pub struct GermlineAssignment {
    species: String,
    v_gene: Option<String>,
    v_identity: Option<f32>,
    j_gene: Option<String>,
    j_identity: Option<f32>,
}

#[pymethods]
impl GermlineAssignment {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("species", &self.species)?;
        dict.set_item("v_gene", self.v_gene.as_deref())?;
        dict.set_item("v_identity", self.v_identity)?;
        dict.set_item("j_gene", self.j_gene.as_deref())?;
        dict.set_item("j_identity", self.j_identity)?;
        Ok(dict)
    }

    fn __repr__(&self) -> String {
        format!(
            "GermlineAssignment({} v={} j={})",
            self.species,
            self.v_gene.as_deref().unwrap_or("-"),
            self.j_gene.as_deref().unwrap_or("-")
        )
    }
}

#[pyclass(frozen, get_all, module = "anarcism._anarcism")]
#[derive(Debug)]
pub struct DomainResult {
    domain_index: usize,
    receptor_type: String,
    chain_type: String,
    species: String,
    start: usize,
    end: usize,
    bit_score: f32,
    e_value: f64,
    bias: f32,
    query_start: usize,
    query_end: usize,
    numbering: Vec<Py<NumberedResidue>>,
    padded_imgt_alignment: String,
    alternative_hits: Vec<Py<ProfileHit>>,
    germline: Option<Py<GermlineAssignment>>,
}

#[pymethods]
impl DomainResult {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("domain_index", self.domain_index)?;
        dict.set_item("receptor_type", &self.receptor_type)?;
        dict.set_item("chain_type", &self.chain_type)?;
        dict.set_item("species", &self.species)?;
        dict.set_item("start", self.start)?;
        dict.set_item("end", self.end)?;
        dict.set_item("bit_score", self.bit_score)?;
        dict.set_item("e_value", self.e_value)?;
        dict.set_item("bias", self.bias)?;
        dict.set_item("query_start", self.query_start)?;
        dict.set_item("query_end", self.query_end)?;
        let numbering = self
            .numbering
            .iter()
            .map(|residue| residue.get().to_dict(py))
            .collect::<PyResult<Vec<_>>>()?;
        dict.set_item("numbering", numbering)?;
        dict.set_item("padded_imgt_alignment", &self.padded_imgt_alignment)?;
        let hits = self
            .alternative_hits
            .iter()
            .map(|hit| hit.get().to_dict(py))
            .collect::<PyResult<Vec<_>>>()?;
        dict.set_item("alternative_hits", hits)?;
        match &self.germline {
            Some(germline) => dict.set_item("germline", germline.get().to_dict(py)?)?,
            None => dict.set_item("germline", py.None())?,
        }
        Ok(dict)
    }

    fn __repr__(&self) -> String {
        format!(
            "DomainResult({} {} {}..{} {:.4} bits)",
            self.chain_type, self.species, self.start, self.end, self.bit_score
        )
    }
}

#[pyclass(frozen, get_all, module = "anarcism._anarcism")]
#[derive(Debug)]
pub struct SequenceResult {
    id: String,
    normalized_sequence: String,
    domains: Vec<Py<DomainResult>>,
    warnings: Vec<String>,
}

#[pymethods]
impl SequenceResult {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("id", &self.id)?;
        dict.set_item("normalized_sequence", &self.normalized_sequence)?;
        let domains = self
            .domains
            .iter()
            .map(|domain| domain.get().to_dict(py))
            .collect::<PyResult<Vec<_>>>()?;
        dict.set_item("domains", domains)?;
        dict.set_item("warnings", self.warnings.clone())?;
        Ok(dict)
    }

    fn __repr__(&self) -> String {
        format!(
            "SequenceResult(id={:?}, domains={})",
            self.id,
            self.domains.len()
        )
    }
}

#[pyclass(frozen, get_all, module = "anarcism._anarcism")]
#[derive(Debug)]
pub struct PairValidationResult {
    ok: bool,
    vh: Option<Py<DomainResult>>,
    vl: Option<Py<DomainResult>>,
    errors: Vec<String>,
}

#[pymethods]
impl PairValidationResult {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("ok", self.ok)?;
        match &self.vh {
            Some(domain) => dict.set_item("vh", domain.get().to_dict(py)?)?,
            None => dict.set_item("vh", py.None())?,
        }
        match &self.vl {
            Some(domain) => dict.set_item("vl", domain.get().to_dict(py)?)?,
            None => dict.set_item("vl", py.None())?,
        }
        dict.set_item("errors", self.errors.clone())?;
        Ok(dict)
    }

    fn __bool__(&self) -> bool {
        self.ok
    }

    fn __repr__(&self) -> String {
        format!(
            "PairValidationResult(ok={}, errors={})",
            if self.ok { "True" } else { "False" },
            self.errors.len()
        )
    }
}

fn residue_to_py(py: Python<'_>, value: &CoreResidue) -> PyResult<Py<NumberedResidue>> {
    Py::new(
        py,
        NumberedResidue {
            sequence_index: value.sequence_index,
            amino_acid: value.amino_acid.to_string(),
            position: value.position,
            insertion_code: value.insertion_code.clone(),
            region: value.region.as_str().to_owned(),
        },
    )
}

fn hit_to_py(py: Python<'_>, value: &CoreHit) -> PyResult<Py<ProfileHit>> {
    Py::new(
        py,
        ProfileHit {
            profile: value.profile.clone(),
            chain_type: value.chain_type.as_str().to_owned(),
            species: value.species.clone(),
            bit_score: value.bit_score,
            e_value: value.e_value,
            bias: value.bias,
            query_start: value.query_start,
            query_end: value.query_end,
        },
    )
}

fn germline_to_py(py: Python<'_>, value: &CoreGermline) -> PyResult<Py<GermlineAssignment>> {
    Py::new(
        py,
        GermlineAssignment {
            species: value.species.clone(),
            v_gene: value.v_gene.clone(),
            v_identity: value.v_identity,
            j_gene: value.j_gene.clone(),
            j_identity: value.j_identity,
        },
    )
}

fn domain_to_py(py: Python<'_>, value: &CoreDomain) -> PyResult<Py<DomainResult>> {
    let numbering = value
        .numbering
        .iter()
        .map(|residue| residue_to_py(py, residue))
        .collect::<PyResult<Vec<_>>>()?;
    let alternative_hits = value
        .alternative_hits
        .iter()
        .map(|hit| hit_to_py(py, hit))
        .collect::<PyResult<Vec<_>>>()?;
    let germline = match &value.germline {
        Some(germline) => Some(germline_to_py(py, germline)?),
        None => None,
    };
    Py::new(
        py,
        DomainResult {
            domain_index: value.domain_index,
            receptor_type: value.receptor_type.as_str().to_owned(),
            chain_type: value.chain_type.as_str().to_owned(),
            species: value.species.clone(),
            start: value.start,
            end: value.end,
            bit_score: value.bit_score,
            e_value: value.e_value,
            bias: value.bias,
            query_start: value.query_start,
            query_end: value.query_end,
            numbering,
            padded_imgt_alignment: value.padded_imgt_alignment.clone(),
            alternative_hits,
            germline,
        },
    )
}

fn sequence_to_py(py: Python<'_>, value: &CoreSequence) -> PyResult<Py<SequenceResult>> {
    let domains = value
        .domains
        .iter()
        .map(|domain| domain_to_py(py, domain))
        .collect::<PyResult<Vec<_>>>()?;
    Py::new(
        py,
        SequenceResult {
            id: value.id.clone(),
            normalized_sequence: value.normalized_sequence.clone(),
            domains,
            warnings: value.warnings.clone(),
        },
    )
}

fn pair_to_py(py: Python<'_>, value: &CorePair) -> PyResult<Py<PairValidationResult>> {
    let vh = match &value.vh {
        Some(domain) => Some(domain_to_py(py, domain)?),
        None => None,
    };
    let vl = match &value.vl {
        Some(domain) => Some(domain_to_py(py, domain)?),
        None => None,
    };
    Py::new(
        py,
        PairValidationResult {
            ok: value.ok,
            vh,
            vl,
            errors: value.errors.clone(),
        },
    )
}

fn sequences_to_py(py: Python<'_>, values: &[CoreSequence]) -> PyResult<Vec<Py<SequenceResult>>> {
    values
        .iter()
        .map(|value| sequence_to_py(py, value))
        .collect()
}

/// Number one sequence and return its detected domains.
#[pyfunction]
#[pyo3(signature = (
    sequence,
    *,
    id = "sequence".to_owned(),
    allowed_chains = None,
    allowed_species = None,
    min_bit_score = None,
    alternative_hit_count = None,
    assign_germline = false,
))]
#[expect(
    clippy::too_many_arguments,
    reason = "keyword options mirror the JS API"
)]
fn number_sequence(
    py: Python<'_>,
    sequence: String,
    id: String,
    allowed_chains: Option<Vec<String>>,
    allowed_species: Option<Vec<String>>,
    min_bit_score: Option<f32>,
    alternative_hit_count: Option<usize>,
    assign_germline: bool,
) -> PyResult<Py<SequenceResult>> {
    let options = numbering_options(
        allowed_chains,
        allowed_species,
        min_bit_score,
        alternative_hit_count,
        assign_germline,
    )?;
    let result = py
        .detach(|| number_sequence_with_id(&id, &sequence, &options))
        .map_err(|error| to_py_error(py, &error))?;
    sequence_to_py(py, &result)
}

/// Number a batch of `(id, sequence)` pairs, preserving input order.
///
/// `workers` above 1 spreads the batch over that many threads and returns
/// results identical to the serial path.
#[pyfunction]
#[pyo3(signature = (
    inputs,
    *,
    workers = None,
    allowed_chains = None,
    allowed_species = None,
    min_bit_score = None,
    alternative_hit_count = None,
    assign_germline = false,
))]
#[expect(
    clippy::too_many_arguments,
    reason = "keyword options mirror the JS API"
)]
fn number_sequences(
    py: Python<'_>,
    inputs: Vec<(String, String)>,
    workers: Option<usize>,
    allowed_chains: Option<Vec<String>>,
    allowed_species: Option<Vec<String>>,
    min_bit_score: Option<f32>,
    alternative_hit_count: Option<usize>,
    assign_germline: bool,
) -> PyResult<Vec<Py<SequenceResult>>> {
    let options = numbering_options(
        allowed_chains,
        allowed_species,
        min_bit_score,
        alternative_hit_count,
        assign_germline,
    )?;
    let records: Vec<SequenceInput> = inputs
        .into_iter()
        .map(|(id, sequence)| SequenceInput { id, sequence })
        .collect();
    let worker_count =
        match workers {
            None => None,
            Some(value) => Some(NonZeroUsize::new(value).ok_or_else(|| {
                PyValueError::new_err("workers must be a positive integer or None")
            })?),
        };
    let results = py
        .detach(|| match worker_count {
            Some(count) => core_number_sequences_parallel(&records, &options, count),
            None => core_number_sequences(&records, &options),
        })
        .map_err(|error| to_py_error(py, &error))?;
    sequences_to_py(py, &results)
}

/// Number every record in a FASTA document.
#[pyfunction]
#[pyo3(signature = (
    fasta,
    *,
    allowed_chains = None,
    allowed_species = None,
    min_bit_score = None,
    alternative_hit_count = None,
    assign_germline = false,
))]
fn number_fasta(
    py: Python<'_>,
    fasta: String,
    allowed_chains: Option<Vec<String>>,
    allowed_species: Option<Vec<String>>,
    min_bit_score: Option<f32>,
    alternative_hit_count: Option<usize>,
    assign_germline: bool,
) -> PyResult<Vec<Py<SequenceResult>>> {
    let options = numbering_options(
        allowed_chains,
        allowed_species,
        min_bit_score,
        alternative_hit_count,
        assign_germline,
    )?;
    let results = py
        .detach(|| core_number_fasta(&fasta, &options))
        .map_err(|error| to_py_error(py, &error))?;
    sequences_to_py(py, &results)
}

/// Check that a heavy/light pair each contain exactly one well-placed domain.
#[pyfunction]
#[pyo3(signature = (
    vh,
    vl,
    *,
    start_max = None,
    end_min = None,
    allowed_chains = None,
    allowed_species = None,
    min_bit_score = None,
    alternative_hit_count = None,
    assign_germline = false,
))]
#[expect(
    clippy::too_many_arguments,
    reason = "keyword options mirror the JS API"
)]
fn validate_antibody_pair(
    py: Python<'_>,
    vh: String,
    vl: String,
    start_max: Option<usize>,
    end_min: Option<usize>,
    allowed_chains: Option<Vec<String>>,
    allowed_species: Option<Vec<String>>,
    min_bit_score: Option<f32>,
    alternative_hit_count: Option<usize>,
    assign_germline: bool,
) -> PyResult<Py<PairValidationResult>> {
    let defaults = PairValidationOptions::default();
    let options = PairValidationOptions {
        numbering: numbering_options(
            allowed_chains,
            allowed_species,
            min_bit_score,
            alternative_hit_count,
            assign_germline,
        )?,
        start_max: start_max.unwrap_or(defaults.start_max),
        end_min: end_min.unwrap_or(defaults.end_min),
    };
    let result = py
        .detach(|| core_validate_pair(&vh, &vl, &options))
        .map_err(|error| to_py_error(py, &error))?;
    pair_to_py(py, &result)
}

/// Chain types present in the embedded profile inventory.
#[pyfunction]
fn chains(py: Python<'_>) -> PyResult<Vec<&'static str>> {
    let database = embedded_profiles().map_err(|error| to_py_error(py, &error))?;
    let mut result = Vec::new();
    for profile in database.profiles() {
        let name = profile.chain_type().as_str();
        if !result.contains(&name) {
            result.push(name);
        }
    }
    Ok(result)
}

/// Species present in the embedded profile inventory.
#[pyfunction]
fn species(py: Python<'_>) -> PyResult<Vec<String>> {
    let database = embedded_profiles().map_err(|error| to_py_error(py, &error))?;
    let mut result: Vec<String> = Vec::new();
    for profile in database.profiles() {
        if !result.iter().any(|value| value == profile.species()) {
            result.push(profile.species().to_owned());
        }
    }
    Ok(result)
}

#[pymodule]
fn _anarcism(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    module.add("AnarcismError", module.py().get_type::<AnarcismError>())?;
    module.add_class::<NumberedResidue>()?;
    module.add_class::<ProfileHit>()?;
    module.add_class::<GermlineAssignment>()?;
    module.add_class::<DomainResult>()?;
    module.add_class::<SequenceResult>()?;
    module.add_class::<PairValidationResult>()?;
    module.add_function(wrap_pyfunction!(number_sequence, module)?)?;
    module.add_function(wrap_pyfunction!(number_sequences, module)?)?;
    module.add_function(wrap_pyfunction!(number_fasta, module)?)?;
    module.add_function(wrap_pyfunction!(validate_antibody_pair, module)?)?;
    module.add_function(wrap_pyfunction!(chains, module)?)?;
    module.add_function(wrap_pyfunction!(species, module)?)?;
    Ok(())
}
