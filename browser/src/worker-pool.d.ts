import type {
	AbortOptions,
	ChainType,
	NumberingOptions,
	PairValidationOptions,
	PairValidationResult,
	SequenceResult,
} from "./index.js";

export const MAX_WORKER_COUNT: number;

export interface WorkerPoolOptions extends AbortOptions {
	/** Initial worker count; defaults to `recommendedWorkerCount()`. */
	workers?: number;
	/** Worker script URL; defaults to the `worker.js` shipped next to this module. */
	workerUrl?: string | URL;
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

export class AnarcismWorkerPool {
	private constructor();
	/** Starts the workers and resolves once every engine is ready. */
	static create(options?: WorkerPoolOptions): Promise<AnarcismWorkerPool>;
	readonly busy: boolean;
	readonly size: number;
	readonly version: string;
	chains(): ChainType[];
	species(): string[];
	resize(workerCount: number, options?: AbortOptions): Promise<WorkerPoolMetadata>;
	numberSequence(
		sequence: string,
		options?: NumberingOptions & AbortOptions,
	): Promise<SequenceResult>;
	numberSequences(
		inputs: SequenceInput[],
		options?: NumberingOptions & AbortOptions,
	): Promise<SequenceResult[]>;
	numberFasta(fasta: string, options?: NumberingOptions & AbortOptions): Promise<SequenceResult[]>;
	validateAntibodyPair(
		vh: string,
		vl: string,
		options?: PairValidationOptions & AbortOptions,
	): Promise<PairValidationResult>;
	terminate(): void;
}
