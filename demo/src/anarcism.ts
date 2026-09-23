import { AnarcismWorkerPool } from "@anarcism/worker-pool.js";
import workerUrl from "@anarcism/worker.js?worker&url";

export {
	MAX_WORKER_COUNT,
	normalizeWorkerCount,
	parseFasta,
	recommendedWorkerCount,
} from "@anarcism/worker-pool.js";

export type { WorkerPoolMetadata } from "@anarcism/worker-pool.js";
export type AnalysisPool = AnarcismWorkerPool;
export type {
	ChainType,
	DomainResult,
	NumberingOptions,
	Region,
	SequenceResult,
} from "@anarcism/index.js";

// Vite must supply the worker URL after bundling.
export function createAnalysisPool(workers: number): Promise<AnarcismWorkerPool> {
	return AnarcismWorkerPool.create({ workerUrl, workers });
}

export function errorCode(error: unknown): string {
	return typeof error === "object" && error !== null && "code" in error
		? String((error as { code: unknown }).code)
		: "error";
}

export function errorMessage(error: unknown): string {
	return error instanceof Error ? error.message : String(error);
}
