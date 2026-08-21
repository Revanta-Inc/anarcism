import init from "./index.js";

const initialization = initialize();

self.addEventListener("message", async (event) => {
  const request = event.data;
  if (request?.type !== "number" && request?.type !== "numberBatch") return;

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
    const results = request.type === "numberBatch"
      ? engine.numberSequences(request.inputs, request.options)
      : [engine.numberSequence(request.text, request.options)];
    self.postMessage({
      type: "result",
      requestId: request.requestId,
      results,
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
    const engine = await init();
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
