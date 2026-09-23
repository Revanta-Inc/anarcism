//! Minimal JSON-over-C-ABI bridge for JavaScript hosts (browsers and Node.js).
//! The JavaScript wrapper is public; raw exports are internal.

use anarcism_core::{
    ChainType, Error, NumberingOptions, PairValidationOptions, SequenceInput, embedded_profiles,
    number_fasta, number_sequence, number_sequences, validate_antibody_pair,
};
use serde::{Deserialize, Serialize};

const API_VERSION: u32 = 1;
const MAX_REQUEST_BYTES: usize = 12 * 1024 * 1024;

#[derive(Debug, Deserialize)]
#[serde(tag = "method", rename_all = "camelCase")]
enum Request {
    Metadata,
    NumberSequence {
        sequence: String,
        #[serde(default)]
        options: NumberingOptions,
    },
    NumberSequences {
        inputs: Vec<SequenceInput>,
        #[serde(default)]
        options: NumberingOptions,
    },
    NumberFasta {
        fasta: String,
        #[serde(default)]
        options: NumberingOptions,
    },
    ValidateAntibodyPair {
        vh: String,
        vl: String,
        #[serde(default)]
        options: PairValidationOptions,
    },
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ModuleMetadata {
    version: &'static str,
    chains: Vec<ChainType>,
    species: Vec<String>,
}

#[derive(Serialize)]
struct Success<T> {
    ok: bool,
    value: T,
}

#[derive(Serialize)]
struct Failure<E> {
    ok: bool,
    error: E,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct RequestError {
    code: &'static str,
    message: String,
}

#[doc(hidden)]
pub fn handle_request(bytes: &[u8]) -> Vec<u8> {
    if bytes.len() > MAX_REQUEST_BYTES {
        return encode_failure(&RequestError {
            code: "REQUEST_TOO_LARGE",
            message: format!(
                "API request contains {} bytes; limit is {MAX_REQUEST_BYTES}",
                bytes.len()
            ),
        });
    }
    let request = match serde_json::from_slice::<Request>(bytes) {
        Ok(request) => request,
        Err(error) => {
            return encode_failure(&RequestError {
                code: "INVALID_REQUEST",
                message: format!(
                    "invalid API request at line {}, column {}",
                    error.line(),
                    error.column()
                ),
            });
        }
    };

    match request {
        Request::Metadata => encode_result(module_metadata()),
        Request::NumberSequence { sequence, options } => {
            encode_result(number_sequence(&sequence, &options))
        }
        Request::NumberSequences { inputs, options } => {
            encode_result(number_sequences(&inputs, &options))
        }
        Request::NumberFasta { fasta, options } => encode_result(number_fasta(&fasta, &options)),
        Request::ValidateAntibodyPair { vh, vl, options } => {
            encode_result(validate_antibody_pair(&vh, &vl, &options))
        }
    }
}

fn module_metadata() -> Result<ModuleMetadata, Error> {
    let database = embedded_profiles()?;
    let mut chains = Vec::new();
    let mut species = Vec::new();
    for profile in database.profiles() {
        if !chains.contains(&profile.chain_type()) {
            chains.push(profile.chain_type());
        }
        if !species.iter().any(|value| value == profile.species()) {
            species.push(profile.species().to_owned());
        }
    }
    Ok(ModuleMetadata {
        version: env!("CARGO_PKG_VERSION"),
        chains,
        species,
    })
}

fn encode_result<T: Serialize>(result: Result<T, Error>) -> Vec<u8> {
    match result {
        Ok(value) => serde_json::to_vec(&Success { ok: true, value }),
        Err(error) => serde_json::to_vec(&Failure { ok: false, error }),
    }
    .unwrap_or_else(|_| {
        br#"{"ok":false,"error":{"code":"INTERNAL","message":"response serialization failed"}}"#
            .to_vec()
    })
}

fn encode_failure<T: Serialize>(error: &T) -> Vec<u8> {
    serde_json::to_vec(&Failure { ok: false, error }).unwrap_or_else(|_| {
        br#"{"ok":false,"error":{"code":"INTERNAL","message":"response serialization failed"}}"#
            .to_vec()
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn anarcism_api_version() -> u32 {
    API_VERSION
}

#[cfg(target_arch = "wasm32")]
#[unsafe(no_mangle)]
pub extern "C" fn anarcism_alloc(length: u32) -> u32 {
    let Ok(length) = usize::try_from(length) else {
        return 0;
    };
    if length == 0 || length > MAX_REQUEST_BYTES {
        return 0;
    }
    let allocation = vec![0_u8; length].into_boxed_slice();
    Box::into_raw(allocation).cast::<u8>() as u32
}

#[cfg(target_arch = "wasm32")]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn anarcism_free(pointer: u32, length: u32) {
    if pointer == 0 || length == 0 {
        return;
    }
    let raw = std::ptr::slice_from_raw_parts_mut(pointer as *mut u8, length as usize);
    // SAFETY: The JavaScript wrapper passes exactly the pointer and length returned
    // by anarcism_alloc or anarcism_call, and relinquishes each allocation once.
    drop(unsafe { Box::from_raw(raw) });
}

#[cfg(target_arch = "wasm32")]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn anarcism_call(pointer: u32, length: u32) -> u64 {
    if pointer == 0 || length == 0 {
        return pack_response(encode_failure(&RequestError {
            code: "INVALID_REQUEST",
            message: "API request is empty".into(),
        }));
    }
    let length = length as usize;
    if length > MAX_REQUEST_BYTES {
        return pack_response(encode_failure(&RequestError {
            code: "REQUEST_TOO_LARGE",
            message: format!("API request contains {length} bytes; limit is {MAX_REQUEST_BYTES}"),
        }));
    }
    // SAFETY: The JavaScript wrapper creates this region with anarcism_alloc and
    // writes `length` bytes before making the call.
    let request = unsafe { std::slice::from_raw_parts(pointer as *const u8, length) };
    pack_response(handle_request(request))
}

#[cfg(target_arch = "wasm32")]
fn pack_response(response: Vec<u8>) -> u64 {
    let output = response.into_boxed_slice();
    let length = output.len() as u32;
    let pointer = Box::into_raw(output).cast::<u8>() as u32;
    (u64::from(length) << 32) | u64::from(pointer)
}

#[cfg(test)]
mod tests {
    use super::*;

    const VH: &str = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";

    #[test]
    fn json_bridge_returns_structured_success_and_failure() {
        let request = format!(r#"{{"method":"numberSequence","sequence":"{VH}"}}"#);
        let success: serde_json::Value =
            serde_json::from_slice(&handle_request(request.as_bytes()))
                .expect("bridge response is JSON");
        assert_eq!(success["ok"], true);
        assert_eq!(success["value"]["domains"][0]["chainType"], "H");

        let failure: serde_json::Value =
            serde_json::from_slice(&handle_request(b"not json")).expect("bridge error is JSON");
        assert_eq!(failure["ok"], false);
        assert_eq!(failure["error"]["code"], "INVALID_REQUEST");
    }

    #[test]
    fn metadata_comes_from_the_embedded_profile_inventory() {
        let response: serde_json::Value =
            serde_json::from_slice(&handle_request(br#"{"method":"metadata"}"#))
                .expect("bridge response is JSON");
        assert_eq!(response["ok"], true);
        assert_eq!(response["value"]["version"], env!("CARGO_PKG_VERSION"));
        assert_eq!(
            response["value"]["chains"],
            serde_json::json!(["H", "K", "L", "A", "B", "G", "D"])
        );
        assert_eq!(
            response["value"]["species"],
            serde_json::json!([
                "human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca", "cow"
            ])
        );
    }

    #[test]
    fn pair_options_deserialize_with_partial_defaults() {
        let request = format!(
            r#"{{"method":"validateAntibodyPair","vh":"{VH}","vl":"{VH}","options":{{"startMax":5}}}}"#
        );
        let response: serde_json::Value =
            serde_json::from_slice(&handle_request(request.as_bytes()))
                .expect("bridge response is JSON");
        assert_eq!(response["ok"], true);
        assert_eq!(response["value"]["ok"], false);
    }
}
