const encoder = new TextEncoder();
const decoder = new TextDecoder();
let wasm;

export class AnarcismError extends Error {
  constructor(error) {
    super(error?.message ?? "anarcism failed");
    this.name = "AnarcismError";
    this.code = error?.code ?? "INTERNAL";
    this.inputId = error?.inputId;
  }
}

function useInstance(instance) {
  const exports = instance?.exports;
  if (
    !(exports?.memory instanceof WebAssembly.Memory) ||
    typeof exports.anarcism_alloc !== "function" ||
    typeof exports.anarcism_free !== "function" ||
    typeof exports.anarcism_call !== "function" ||
    exports.anarcism_api_version() !== 1
  ) {
    throw new TypeError("incompatible anarcism WebAssembly module");
  }
  wasm = exports;
  return api;
}

export function initSync(source) {
  const module = source instanceof WebAssembly.Module
    ? source
    : new WebAssembly.Module(toBytes(source));
  return useInstance(new WebAssembly.Instance(module, {}));
}

export async function init(source = new URL("./anarcism.wasm", import.meta.url)) {
  if (source instanceof WebAssembly.Module) {
    return useInstance(await WebAssembly.instantiate(source, {}));
  }
  if (source instanceof Response) {
    return instantiateResponse(source);
  }
  if (typeof source === "string" || source instanceof URL || source instanceof Request) {
    return instantiateResponse(await fetch(source));
  }
  const result = await WebAssembly.instantiate(toBytes(source), {});
  return useInstance(result.instance);
}

async function instantiateResponse(response) {
  if (!response.ok) {
    throw new Error(`could not load anarcism WASM: HTTP ${response.status}`);
  }
  // ArrayBuffer instantiation is consistently supported by the three major
  // engines and still compiles asynchronously. At this package size it avoids
  // browser-specific streaming/Response-clone lifecycle edge cases.
  const result = await WebAssembly.instantiate(await response.arrayBuffer(), {});
  return useInstance(result.instance);
}

function toBytes(source) {
  if (source instanceof ArrayBuffer) return source;
  if (ArrayBuffer.isView(source)) {
    return new Uint8Array(source.buffer, source.byteOffset, source.byteLength);
  }
  throw new TypeError("expected WebAssembly bytes, module, response, or URL");
}

function invoke(request) {
  if (!wasm) {
    throw new AnarcismError({
      code: "NOT_INITIALIZED",
      message: "call init() or initSync() before numbering sequences",
    });
  }
  const input = encoder.encode(JSON.stringify(request));
  const inputPointer = wasm.anarcism_alloc(input.byteLength);
  if (!inputPointer) {
    throw new AnarcismError({
      code: "REQUEST_TOO_LARGE",
      message: "request exceeds the WebAssembly input limit",
    });
  }
  let packed;
  try {
    new Uint8Array(wasm.memory.buffer, inputPointer, input.byteLength).set(input);
    packed = wasm.anarcism_call(inputPointer, input.byteLength);
  } finally {
    wasm.anarcism_free(inputPointer, input.byteLength);
  }

  const outputPointer = Number(packed & 0xffff_ffffn);
  const outputLength = Number(packed >> 32n);
  let response;
  try {
    const output = new Uint8Array(wasm.memory.buffer, outputPointer, outputLength);
    response = JSON.parse(decoder.decode(output));
  } finally {
    wasm.anarcism_free(outputPointer, outputLength);
  }
  if (!response.ok) throw new AnarcismError(response.error);
  return response.value;
}

export function numberSequence(sequence, options = {}) {
  return invoke({ method: "numberSequence", sequence, options });
}

export function numberSequences(inputs, options = {}) {
  return invoke({ method: "numberSequences", inputs, options });
}

export function numberFasta(fasta, options = {}) {
  return invoke({ method: "numberFasta", fasta, options });
}

export function validateAntibodyPair(vh, vl, options = {}) {
  return invoke({ method: "validateAntibodyPair", vh, vl, options });
}

const api = Object.freeze({
  numberSequence,
  numberSequences,
  numberFasta,
  validateAntibodyPair,
});

export default init;
