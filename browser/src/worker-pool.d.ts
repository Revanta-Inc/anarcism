import type { ChainType, NumberingOptions, SequenceResult } from "./index.js";

export const MAX_WORKER_COUNT: number;

export interface WorkerPoolOptions {
	maxWorkers?: number;
	workerFactory?: (url: string | URL, options: WorkerOptions) => Worker;
}

export interface WorkerPoolMetadata {
	initMs: number;
	version: string;
	chains: ChainType[];
	species: string[];
	workerCount: number;
}

export interface SequenceInput {
	id: string;
	sequence: string;
}

export function recommendedWorkerCount(hardwareConcurrency?: number): number;
export function normalizeWorkerCount(value: number, maxWorkers?: number): number;
export function parseFasta(input: string): SequenceInput[];

export class AnalysisWorkerPool {
	constructor(workerUrl?: string | URL, options?: WorkerPoolOptions);
	readonly busy: boolean;
	readonly size: number;
	initialize(workerCount?: number): Promise<WorkerPoolMetadata>;
	resize(workerCount: number): Promise<WorkerPoolMetadata>;
	numberSequence(sequence: string, options?: NumberingOptions): Promise<SequenceResult>;
	numberSequences(inputs: SequenceInput[], options?: NumberingOptions): Promise<SequenceResult[]>;
	numberFasta(fasta: string, options?: NumberingOptions): Promise<SequenceResult[]>;
	terminate(): void;
}

export { AnalysisWorkerPool as AnarcismWorkerPool };
