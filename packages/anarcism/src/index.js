const encoder = new TextEncoder();
const decoder = new TextDecoder();

export class AnarcismError extends Error {
	constructor(error) {
		super(error?.message ?? "anarcism failed");
		this.name = "AnarcismError";
		this.code = error?.code ?? "INTERNAL";
		this.inputId = error?.inputId;
	}
}

export class Anarcism {
	#wasm;
	#metadata;

	constructor(instance) {
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
		this.#wasm = exports;
		this.#metadata = this.#invoke({ method: "metadata" });
		Object.freeze(this);
	}

	static async create({ source = new URL("./anarcism.wasm", import.meta.url), signal } = {}) {
		signal?.throwIfAborted();
		if (source instanceof WebAssembly.Module) {
			const instance = await abortable(WebAssembly.instantiate(source, {}), signal);
			return new Anarcism(instance);
		}
		// Node.js cannot fetch file: URLs, which is how the bundled module resolves there.
		const nodeFs = globalThis.process?.getBuiltinModule?.("node:fs/promises");
		if (nodeFs && source instanceof URL && source.protocol === "file:") {
			source = await nodeFs.readFile(source, { signal });
		}
		if (typeof source === "string" || source instanceof URL || source instanceof Request) {
			source = await fetch(source, { signal });
		}
		if (source instanceof Response) {
			if (!source.ok) {
				throw new Error(`could not load anarcism WASM: HTTP ${source.status}`);
			}
			source = await abortable(source.arrayBuffer(), signal);
		}
		const result = await abortable(WebAssembly.instantiate(source, {}), signal);
		return new Anarcism(result.instance);
	}

	static createSync(source) {
		const module = source instanceof WebAssembly.Module ? source : new WebAssembly.Module(source);
		return new Anarcism(new WebAssembly.Instance(module, {}));
	}

	get version() {
		return this.#metadata.version;
	}

	chains() {
		return [...this.#metadata.chains];
	}

	species() {
		return [...this.#metadata.species];
	}

	numberSequence(sequence, options = {}) {
		return this.#invoke({ method: "numberSequence", sequence, options });
	}

	numberSequences(inputs, options = {}) {
		return this.#invoke({ method: "numberSequences", inputs, options });
	}

	numberFasta(fasta, options = {}) {
		return this.#invoke({ method: "numberFasta", fasta, options });
	}

	validateAntibodyPair(vh, vl, options = {}) {
		return this.#invoke({ method: "validateAntibodyPair", vh, vl, options });
	}

	#invoke(request) {
		const wasm = this.#wasm;
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
}

// WebAssembly compilation cannot be cancelled; an abort discards its result.
function abortable(promise, signal) {
	if (!signal) return promise;
	return new Promise((resolve, reject) => {
		const onAbort = () => reject(signal.reason);
		if (signal.aborted) {
			onAbort();
			return;
		}
		signal.addEventListener("abort", onAbort, { once: true });
		promise.then(resolve, reject).finally(() => signal.removeEventListener("abort", onAbort));
	});
}
