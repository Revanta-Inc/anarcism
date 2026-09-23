import { AnalysisWorkerPool } from "@anarcism/worker-pool.js";
import workerUrl from "@anarcism/worker.js?worker&url";

export {
	MAX_WORKER_COUNT,
	normalizeWorkerCount,
	parseFasta,
	recommendedWorkerCount,
} from "@anarcism/worker-pool.js";

export type { WorkerPoolMetadata } from "@anarcism/worker-pool.js";
export type AnalysisPool = AnalysisWorkerPool;
export type {
	ChainType,
	DomainResult,
	NumberingOptions,
	Region,
	SequenceResult,
} from "@anarcism/index.js";

// The published pool resolves `./worker.js` against its own module URL, which
// only holds while the package stays unbundled. Vite emits the worker as a
// separate bundle, so the pool is handed that URL instead.
export function createAnalysisPool(): AnalysisWorkerPool {
	return new AnalysisWorkerPool(workerUrl);
}

export function errorCode(error: unknown): string {
	return typeof error === "object" && error !== null && "code" in error
		? String((error as { code: unknown }).code)
		: "error";
}

export function errorMessage(error: unknown): string {
	return error instanceof Error ? error.message : String(error);
}
