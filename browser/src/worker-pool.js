const DEFAULT_MAX_WORKERS = 16;
const DEFAULT_FALLBACK_WORKERS = 4;
const DEFAULT_AUTO_WORKER_CAP = 8;
const MAX_FASTA_BYTES = 10 * 1024 * 1024;
const MAX_BATCH_SIZE = 1_000;
const utf8Encoder = new TextEncoder();

export const MAX_WORKER_COUNT = DEFAULT_MAX_WORKERS;

export function recommendedWorkerCount(
	hardwareConcurrency = globalThis.navigator?.hardwareConcurrency,
) {
	const available = Number.isFinite(Number(hardwareConcurrency))
		? Number(hardwareConcurrency)
		: DEFAULT_FALLBACK_WORKERS;
	return normalizeWorkerCount(Math.min(available, DEFAULT_AUTO_WORKER_CAP), DEFAULT_MAX_WORKERS);
}

export function normalizeWorkerCount(value, maxWorkers = DEFAULT_MAX_WORKERS) {
	const maximum = Math.max(1, Math.trunc(Number(maxWorkers)) || DEFAULT_MAX_WORKERS);
	const count = Math.trunc(Number(value));
	return Math.min(maximum, Math.max(1, Number.isFinite(count) ? count : 1));
}

export class AnalysisWorkerPool {
	#busy = false;
	#maxWorkers;
	#metadata;
	#workerFactory;
	#workers = [];
	#workerUrl;

	constructor(
		workerUrl = new URL("./worker.js", import.meta.url),
		{
			maxWorkers = DEFAULT_MAX_WORKERS,
			workerFactory = (url, options) => new Worker(url, options),
		} = {},
	) {
		this.#workerUrl = workerUrl;
		this.#maxWorkers = normalizeWorkerCount(maxWorkers, Number.MAX_SAFE_INTEGER);
		this.#workerFactory = workerFactory;
	}

	get busy() {
		return this.#busy;
	}

	get size() {
		return this.#workers.length;
	}

	initialize(workerCount = recommendedWorkerCount()) {
		return this.resize(workerCount);
	}

	async resize(value) {
		if (this.#busy) {
			throw makeError("POOL_BUSY", "cannot resize the worker pool during analysis");
		}

		const started = performance.now();
		const targetSize = normalizeWorkerCount(value, this.#maxWorkers);
		if (targetSize < this.#workers.length) {
			for (const worker of this.#workers.splice(targetSize)) worker.terminate();
		}

		const originalSize = this.#workers.length;
		try {
			while (this.#workers.length < targetSize) {
				const index = this.#workers.length;
				this.#workers.push(
					new WorkerClient(this.#workerFactory, this.#workerUrl, `anarcism-${index + 1}`),
				);
			}

			const metadata = await Promise.all(this.#workers.map((worker) => worker.ready));
			assertMatchingMetadata(metadata);
			this.#metadata = metadata[0];
		} catch (error) {
			for (const worker of this.#workers.splice(originalSize)) worker.terminate();
			throw error;
		}

		return {
			...this.#metadata,
			initMs: performance.now() - started,
			workerCount: this.#workers.length,
		};
	}

	async numberSequence(sequence, options = {}) {
		const response = await this.#numberSingle(sequence, options);
		return response.results[0];
	}

	async numberSequences(inputs, options = {}) {
		return (await this.#numberInputs(inputs, options)).results;
	}

	async numberFasta(fasta, options = {}) {
		return (await this.#numberInputs(parseFasta(fasta), options)).results;
	}

	terminate() {
		for (const worker of this.#workers.splice(0)) worker.terminate();
		this.#metadata = undefined;
	}

	async #numberSingle(sequence, options) {
		if (!this.#workers.length) {
			throw makeError("NOT_INITIALIZED", "initialize the worker pool before numbering");
		}
		if (this.#busy) {
			throw makeError("POOL_BUSY", "the worker pool is already numbering a batch");
		}

		this.#busy = true;
		try {
			const response = await this.#workers[0].request({
				type: "number",
				text: sequence,
				options,
			});
			return { ...response, workersUsed: 1 };
		} finally {
			this.#busy = false;
		}
	}

	async #numberInputs(inputs, options) {
		if (!this.#workers.length) {
			throw makeError("NOT_INITIALIZED", "initialize the worker pool before numbering");
		}
		if (this.#busy) {
			throw makeError("POOL_BUSY", "the worker pool is already numbering a batch");
		}
		if (inputs.length > MAX_BATCH_SIZE) {
			throw makeError(
				"BATCH_TOO_LARGE",
				`batch contains ${inputs.length} records; configured limit is ${MAX_BATCH_SIZE}`,
			);
		}
		if (!inputs.length) return { results: [], computeMs: 0, workersUsed: 0 };

		this.#busy = true;
		try {
			const workerCount = Math.min(inputs.length, this.#workers.length);
			const shards = balanceInputs(inputs, workerCount);
			const settledResponses = await Promise.allSettled(
				shards.map((shard, index) =>
					this.#workers[index].request({
						type: "numberBatch",
						inputs: shard.map((item) => item.input),
						options,
					}),
				),
			);
			const failedResponse = settledResponses.find((response) => response.status === "rejected");
			if (failedResponse) throw failedResponse.reason;
			const responses = settledResponses.map((response) => response.value);

			const results = new Array(inputs.length);
			for (let shardIndex = 0; shardIndex < shards.length; shardIndex += 1) {
				const shard = shards[shardIndex];
				const shardResults = responses[shardIndex].results;
				if (shardResults.length !== shard.length) {
					throw makeError("INTERNAL", "an analysis worker returned an incomplete batch");
				}
				for (let resultIndex = 0; resultIndex < shard.length; resultIndex += 1) {
					results[shard[resultIndex].index] = shardResults[resultIndex];
				}
			}

			return {
				results,
				computeMs: Math.max(...responses.map((response) => response.computeMs)),
				workersUsed: workerCount,
			};
		} finally {
			this.#busy = false;
		}
	}
}

export { AnalysisWorkerPool as AnarcismWorkerPool };

export function parseFasta(input) {
	const inputBytes = utf8Encoder.encode(input).byteLength;
	if (inputBytes > MAX_FASTA_BYTES) {
		throw makeError(
			"FASTA_TOO_LARGE",
			`FASTA input contains ${inputBytes} bytes; configured limit is ${MAX_FASTA_BYTES}`,
		);
	}

	const records = [];
	let currentId;
	let currentSequence = [];
	const lines = input.split("\n");

	for (let lineIndex = 0; lineIndex < lines.length; lineIndex += 1) {
		const line = lines[lineIndex].replace(/\r+$/u, "");
		if (line.startsWith(">")) {
			if (currentId !== undefined) {
				pushRecord(records, currentId, currentSequence);
			}
			const id = line.slice(1).trim();
			if (!id) {
				throw makeError("INVALID_FASTA", `FASTA header on line ${lineIndex + 1} is empty`);
			}
			currentId = id;
			currentSequence = [];
		} else if (!line.trim()) {
			continue;
		} else if (currentId === undefined) {
			throw makeError(
				"INVALID_FASTA",
				`FASTA sequence data appears before the first header on line ${lineIndex + 1}`,
			);
		} else {
			currentSequence.push(line);
		}
	}

	if (currentId !== undefined) pushRecord(records, currentId, currentSequence);
	if (!records.length) {
		throw makeError("INVALID_FASTA", "FASTA input contains no records");
	}
	return records;
}

function pushRecord(records, id, sequenceParts) {
	if (records.length === MAX_BATCH_SIZE) {
		throw makeError(
			"BATCH_TOO_LARGE",
			`FASTA contains more than the configured ${MAX_BATCH_SIZE} records`,
		);
	}
	if (!sequenceParts.length) {
		throw makeError("INVALID_FASTA", "FASTA record has an empty sequence", id);
	}
	records.push({ id, sequence: sequenceParts.join("") });
}

function balanceInputs(inputs, workerCount) {
	const shards = Array.from({ length: workerCount }, () => []);
	const loads = new Array(workerCount).fill(0);
	const weightedInputs = inputs
		.map((input, index) => ({ index, input }))
		.sort(
			(left, right) =>
				right.input.sequence.length - left.input.sequence.length || left.index - right.index,
		);

	for (const item of weightedInputs) {
		let lightest = 0;
		for (let index = 1; index < loads.length; index += 1) {
			if (loads[index] < loads[lightest]) lightest = index;
		}
		shards[lightest].push(item);
		loads[lightest] += item.input.sequence.length;
	}

	for (const shard of shards) shard.sort((left, right) => left.index - right.index);
	return shards;
}

class WorkerClient {
	#alive = true;
	#initializationState = "pending";
	#nextRequestId = 1;
	#pending = new Map();
	#ready;
	#rejectReady;
	#resolveReady;
	#terminated = false;
	#worker;

	constructor(workerFactory, workerUrl, name) {
		this.#ready = new Promise((resolve, reject) => {
			this.#resolveReady = resolve;
			this.#rejectReady = reject;
		});
		this.#worker = workerFactory(workerUrl, { type: "module", name });
		this.#worker.addEventListener("message", (event) => this.#handleMessage(event.data));
		this.#worker.addEventListener("error", (event) => {
			this.#fail(makeError("WORKER_FAILED", event.message || "an analysis worker failed"));
		});
		this.#worker.addEventListener("messageerror", () => {
			this.#fail(makeError("WORKER_FAILED", "an analysis worker sent an unreadable response"));
		});
	}

	get ready() {
		return this.#ready;
	}

	async request(message) {
		await this.#ready;
		if (!this.#alive) {
			throw makeError("WORKER_FAILED", "the analysis worker is unavailable");
		}

		const requestId = this.#nextRequestId++;
		return new Promise((resolve, reject) => {
			this.#pending.set(requestId, { resolve, reject });
			try {
				this.#worker.postMessage({ ...message, requestId });
			} catch (error) {
				this.#pending.delete(requestId);
				reject(error);
			}
		});
	}

	terminate() {
		if (this.#terminated) return;
		this.#terminated = true;
		this.#alive = false;
		this.#worker.terminate();
		const error = makeError("WORKER_TERMINATED", "the analysis worker was terminated");
		if (this.#initializationState === "pending") {
			this.#initializationState = "failed";
			this.#rejectReady(error);
		}
		for (const pending of this.#pending.values()) pending.reject(error);
		this.#pending.clear();
	}

	#fail(error) {
		if (this.#terminated) return;
		this.#terminated = true;
		this.#alive = false;
		this.#worker.terminate();
		if (this.#initializationState === "pending") {
			this.#initializationState = "failed";
			this.#rejectReady(error);
		}
		for (const pending of this.#pending.values()) pending.reject(error);
		this.#pending.clear();
	}

	#handleMessage(message) {
		if (message?.type === "ready" && this.#initializationState === "pending") {
			this.#initializationState = "ready";
			this.#resolveReady(message);
			return;
		}
		if (message?.type === "initializationError" && this.#initializationState === "pending") {
			this.#fail(asError(message.error));
			return;
		}

		const pending = this.#pending.get(message?.requestId);
		if (!pending) return;
		this.#pending.delete(message.requestId);
		if (message.type === "result") pending.resolve(message);
		else pending.reject(asError(message.error));
	}
}

function assertMatchingMetadata(metadata) {
	const signature = (value) =>
		JSON.stringify({
			version: value.version,
			chains: value.chains,
			species: value.species,
		});
	const expected = signature(metadata[0]);
	if (metadata.some((value) => signature(value) !== expected)) {
		throw makeError("INTERNAL", "analysis workers loaded inconsistent model metadata");
	}
}

function asError(value) {
	return makeError(value?.code ?? "INTERNAL", value?.message ?? "anarcism failed", value?.inputId);
}

function makeError(code, message, inputId) {
	const error = new Error(message);
	error.name = "AnarcismError";
	error.code = code;
	error.inputId = inputId;
	return error;
}
