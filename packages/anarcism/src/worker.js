import { Anarcism } from "./index.js";

const initialization = initialize();

// Map lookup cannot resolve inherited handler names.
const handlers = new Map([
	[
		"number",
		(engine, request) => ({
			results: [engine.numberSequence(request.text, request.options)],
		}),
	],
	[
		"numberBatch",
		(engine, request) => ({
			results: engine.numberSequences(request.inputs, request.options),
		}),
	],
	[
		"validatePair",
		(engine, request) => ({
			value: engine.validateAntibodyPair(request.vh, request.vl, request.options),
		}),
	],
]);

self.addEventListener("message", async (event) => {
	const request = event.data;
	const handler = handlers.get(request?.type);
	if (!handler) return;

	const engine = await initialization;
	if (!engine) {
		self.postMessage({
			type: "error",
			requestId: request.requestId,
			error: {
				code: "NOT_INITIALIZED",
				message: "the analysis worker could not initialize",
			},
		});
		return;
	}

	const started = performance.now();
	try {
		self.postMessage({
			type: "result",
			requestId: request.requestId,
			...handler(engine, request),
			computeMs: performance.now() - started,
		});
	} catch (error) {
		self.postMessage({
			type: "error",
			requestId: request.requestId,
			error: serializeError(error),
		});
	}
});

async function initialize() {
	const started = performance.now();
	try {
		const engine = await Anarcism.create();
		self.postMessage({
			type: "ready",
			initMs: performance.now() - started,
			version: engine.version,
			chains: engine.chains(),
			species: engine.species(),
		});
		return engine;
	} catch (error) {
		self.postMessage({
			type: "initializationError",
			error: serializeError(error),
		});
		return null;
	}
}

function serializeError(error) {
	return {
		code: error?.code ?? "INTERNAL",
		message: error?.message ?? "anarcism failed",
		inputId: error?.inputId,
	};
}
