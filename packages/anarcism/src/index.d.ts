export type ChainType = "H" | "K" | "L" | "A" | "B" | "G" | "D";
export type ReceptorType = "IG" | "TR";
export type Region = "FR1" | "CDR1" | "FR2" | "CDR2" | "FR3" | "CDR3" | "FR4";

export interface NumberingOptions {
	allowedChains?: ChainType[];
	allowedSpecies?: string[];
	minBitScore?: number;
	alternativeHitCount?: number;
	assignGermline?: boolean;
}

export interface NumberedResidue {
	sequenceIndex: number;
	aminoAcid: string;
	position: number;
	insertionCode: string;
	region: Region;
}

export interface ProfileHit {
	profile: string;
	chainType: ChainType;
	species: string;
	bitScore: number;
	eValue: number;
	bias: number;
	queryStart: number;
	queryEnd: number;
}

export interface GermlineAssignment {
	species: string;
	vGene?: string;
	vIdentity?: number;
	jGene?: string;
	jIdentity?: number;
}

export interface DomainResult {
	domainIndex: number;
	receptorType: ReceptorType;
	chainType: ChainType;
	species: string;
	start: number;
	end: number;
	bitScore: number;
	eValue: number;
	bias: number;
	queryStart: number;
	queryEnd: number;
	numbering: NumberedResidue[];
	paddedImgtAlignment: string;
	alternativeHits: ProfileHit[];
	germline?: GermlineAssignment;
}

export interface SequenceResult {
	id: string;
	normalizedSequence: string;
	domains: DomainResult[];
	warnings: string[];
}

export interface PairValidationOptions extends NumberingOptions {
	startMax?: number;
	endMin?: number;
}

export interface PairValidationResult {
	ok: boolean;
	vh?: DomainResult;
	vl?: DomainResult;
	errors: string[];
}

export class AnarcismError extends Error {
	readonly code: string;
	readonly inputId?: string;
}

export type WasmSource =
	| ArrayBuffer
	| ArrayBufferView
	| WebAssembly.Module
	| Response
	| Request
	| URL
	| string;

export interface AbortOptions {
	/** Rejects the call with `signal.reason` when aborted. */
	signal?: AbortSignal;
}

export interface CreateOptions extends AbortOptions {
	/** Where to load the engine from; defaults to the `anarcism.wasm` shipped next to this module. */
	source?: WasmSource;
}

export class Anarcism {
	private constructor();
	static create(options?: CreateOptions): Promise<Anarcism>;
	static createSync(source: ArrayBuffer | ArrayBufferView | WebAssembly.Module): Anarcism;
	readonly version: string;
	chains(): ChainType[];
	species(): string[];
	numberSequence(sequence: string, options?: NumberingOptions): SequenceResult;
	numberSequences(
		inputs: Array<{ id: string; sequence: string }>,
		options?: NumberingOptions,
	): SequenceResult[];
	numberFasta(fasta: string, options?: NumberingOptions): SequenceResult[];
	validateAntibodyPair(
		vh: string,
		vl: string,
		options?: PairValidationOptions,
	): PairValidationResult;
}
