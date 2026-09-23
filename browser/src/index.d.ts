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
	eValue?: number;
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
	eValue?: number;
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

export interface AnarcismApi {
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

export function init(source?: WasmSource): Promise<AnarcismApi>;
export function initSync(source: ArrayBuffer | ArrayBufferView | WebAssembly.Module): AnarcismApi;
export function numberSequence(sequence: string, options?: NumberingOptions): SequenceResult;
export function numberSequences(
	inputs: Array<{ id: string; sequence: string }>,
	options?: NumberingOptions,
): SequenceResult[];
export function numberFasta(fasta: string, options?: NumberingOptions): SequenceResult[];
export function validateAntibodyPair(
	vh: string,
	vl: string,
	options?: PairValidationOptions,
): PairValidationResult;

export default init;
