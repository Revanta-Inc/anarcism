import { readFile } from "node:fs/promises";
import { availableParallelism } from "node:os";
import { Worker, isMainThread, parentPort, workerData } from "node:worker_threads";

import { Anarcism } from "../packages/anarcism/dist/index.js";

const corpusUrl = new URL("../tests/golden/corpus.json", import.meta.url);
const referenceUrl = new URL("../tests/golden/corpus_reference.jsonl", import.meta.url);
const wasmUrl = new URL("../packages/anarcism/dist/anarcism.wasm", import.meta.url);
const maximumAllowedAbsoluteScoreDifference = 0.2;
const maximumAllowedRelativeEValueDifference = 0.08;
const maximumAllowedAbsoluteBiasDifference = 0.2;

function emptyResult() {
	return {
		domains: 0,
		residues: 0,
		exactDomains: 0,
		exactResidues: 0,
		scoresWithinDisplayedPrecision: 0,
		scoreCount: 0,
		absoluteScoreDifferenceSum: 0,
		maximumScoreDifference: 0,
		maximumBiasDifference: 0,
		eValuesWithinDisplayedPrecision: 0,
		eValueCount: 0,
		relativeEValueDifferenceSum: 0,
		maximumRelativeEValueDifference: 0,
		maximumRelativeEValueDifferenceCase: null,
		failures: [],
	};
}

function mergeResult(target, source) {
	target.domains += source.domains;
	target.residues += source.residues;
	target.exactDomains += source.exactDomains;
	target.exactResidues += source.exactResidues;
	target.scoresWithinDisplayedPrecision += source.scoresWithinDisplayedPrecision;
	target.scoreCount += source.scoreCount;
	target.absoluteScoreDifferenceSum += source.absoluteScoreDifferenceSum;
	target.maximumScoreDifference = Math.max(
		target.maximumScoreDifference,
		source.maximumScoreDifference,
	);
	target.maximumBiasDifference = Math.max(
		target.maximumBiasDifference,
		source.maximumBiasDifference,
	);
	target.eValuesWithinDisplayedPrecision += source.eValuesWithinDisplayedPrecision;
	target.eValueCount += source.eValueCount;
	target.relativeEValueDifferenceSum += source.relativeEValueDifferenceSum;
	if (source.maximumRelativeEValueDifference > target.maximumRelativeEValueDifference) {
		target.maximumRelativeEValueDifference = source.maximumRelativeEValueDifference;
		target.maximumRelativeEValueDifferenceCase = source.maximumRelativeEValueDifferenceCase;
	}
	target.failures.push(...source.failures);
}

function recordEValueParity(result, id, actual, expected, significantDigits) {
	result.eValueCount += 1;
	if (!Number.isFinite(actual) || actual <= 0 || !Number.isFinite(expected) || expected <= 0) {
		result.failures.push(`${id}: invalid E-value`);
		return;
	}
	if (!Number.isSafeInteger(significantDigits) || significantDigits < 1) {
		result.failures.push(`${id}: invalid E-value display precision`);
		return;
	}

	const relativeDifference = Math.abs(actual - expected) / expected;
	result.relativeEValueDifferenceSum += relativeDifference;
	if (relativeDifference > result.maximumRelativeEValueDifference) {
		result.maximumRelativeEValueDifference = relativeDifference;
		result.maximumRelativeEValueDifferenceCase = id;
	}

	const exponent = Math.floor(Math.log10(expected));
	const displayQuantum = 10 ** (exponent - significantDigits + 1);
	const displayTolerance = displayQuantum / 2 + expected * 1e-12;
	if (Math.abs(actual - expected) <= displayTolerance) {
		result.eValuesWithinDisplayedPrecision += 1;
	}
}

// Mirrors tools/generate_golden.py: IMGT labels in residue order, with runs of
// consecutive plain positions written as `a-b`. Returns null when residues are
// not the contiguous sequence span the reference encoding assumes.
function encodeNumbering(sequence, domain) {
	const tokens = [];
	let run = null;
	const flush = () => {
		if (run) tokens.push(run[0] === run[1] ? `${run[0]}` : `${run[0]}-${run[1]}`);
		run = null;
	};
	for (const [offset, residue] of domain.numbering.entries()) {
		if (
			residue.sequenceIndex !== domain.start + offset ||
			residue.aminoAcid !== sequence[residue.sequenceIndex]
		) {
			return null;
		}
		if (!residue.insertionCode && run && residue.position === run[1] + 1) {
			run[1] = residue.position;
			continue;
		}
		flush();
		if (residue.insertionCode) tokens.push(`${residue.position}${residue.insertionCode}`);
		else run = [residue.position, residue.position];
	}
	flush();
	return tokens.join(" ");
}

async function readJson(url) {
	return JSON.parse(await readFile(url, "utf8"));
}

// corpus_reference.jsonl: a header record, then one record per case, sorted by id.
async function readReference(url) {
	const [header, ...cases] = (await readFile(url, "utf8"))
		.trimEnd()
		.split("\n")
		.map((line) => JSON.parse(line));
	return { ...header, cases: new Map(cases.map((entry) => [entry.id, entry])) };
}

async function initializeWasm() {
	return Anarcism.create({ source: await readFile(wasmUrl) });
}

async function checkShard(shardIndex, shardCount) {
	const corpusPromise = readJson(corpusUrl);
	const referencePromise = readReference(referenceUrl);
	const anarcism = await initializeWasm();
	const [corpus, reference] = await Promise.all([corpusPromise, referencePromise]);
	const result = emptyResult();

	for (let caseIndex = shardIndex; caseIndex < corpus.length; caseIndex += shardCount) {
		const testCase = corpus[caseIndex];
		const referenceCase = reference.cases.get(testCase.id);
		if (!referenceCase || testCase.category !== referenceCase.category) {
			result.failures.push(`${testCase.id}: reference identity`);
			continue;
		}
		const observed = anarcism.numberSequence(testCase.seq);
		if (observed.domains.length !== referenceCase.domains.length) {
			result.failures.push(`${testCase.id}: domain count`);
			continue;
		}
		for (let index = 0; index < observed.domains.length; index += 1) {
			result.domains += 1;
			const actual = observed.domains[index];
			const expected = referenceCase.domains[index];
			const numberingExact = encodeNumbering(testCase.seq, actual) === expected.numbering;
			const domainExact =
				`${actual.species}_${actual.chainType}` === expected.profile &&
				actual.start === expected.start &&
				actual.end === expected.end &&
				actual.paddedImgtAlignment === expected.alignment;
			if (domainExact) result.exactDomains += 1;
			else result.failures.push(`${testCase.id}: domain ${index}`);
			result.maximumScoreDifference = Math.max(
				result.maximumScoreDifference,
				Math.abs(actual.bitScore - expected.bitScore),
			);
			const scoreDifference = Math.abs(actual.bitScore - expected.bitScore);
			result.scoreCount += 1;
			result.absoluteScoreDifferenceSum += scoreDifference;
			if (scoreDifference <= reference.format.scoreDisplayPrecisionBits / 2 + 1e-6) {
				result.scoresWithinDisplayedPrecision += 1;
			}
			result.maximumBiasDifference = Math.max(
				result.maximumBiasDifference,
				Math.abs(actual.bias - expected.bias),
			);
			recordEValueParity(
				result,
				`${testCase.id}: domain ${index}`,
				actual.eValue,
				expected.eValue,
				reference.format.eValueDisplaySignificantDigits,
			);
			result.residues += actual.numbering.length;
			if (numberingExact) result.exactResidues += actual.numbering.length;
			else result.failures.push(`${testCase.id}: numbering ${index}`);
		}
	}
	return result;
}

function configuredWorkerCount(caseCount) {
	const configured = process.env.PARITY_WORKERS;
	if (configured !== undefined && !/^\d+$/.test(configured)) {
		throw new Error("PARITY_WORKERS must be a positive integer");
	}
	const requested = configured === undefined ? availableParallelism() : Number(configured);
	if (!Number.isSafeInteger(requested) || requested < 1) {
		throw new Error("PARITY_WORKERS must be a positive integer");
	}
	return Math.min(requested, Math.max(caseCount, 1));
}

function runWorker(shardIndex, shardCount) {
	return new Promise((resolve, reject) => {
		const worker = new Worker(new URL(import.meta.url), {
			workerData: { shardIndex, shardCount },
		});
		let receivedResult = false;
		worker.once("message", (result) => {
			receivedResult = true;
			resolve(result);
		});
		worker.once("error", reject);
		worker.once("exit", (code) => {
			if (code !== 0) reject(new Error(`parity worker ${shardIndex} exited with code ${code}`));
			else if (!receivedResult) reject(new Error(`parity worker ${shardIndex} exited silently`));
		});
	});
}

if (!isMainThread) {
	const result = await checkShard(workerData.shardIndex, workerData.shardCount);
	parentPort.postMessage(result);
} else {
	const [corpus, reference] = await Promise.all([readJson(corpusUrl), readReference(referenceUrl)]);
	const result = emptyResult();
	if (corpus.length !== reference.cases.size) {
		result.failures.push("corpus/reference case count");
	}

	const workerCount = configuredWorkerCount(corpus.length);
	const shardResults =
		corpus.length === reference.cases.size
			? Promise.all(
					Array.from({ length: workerCount }, (_, shardIndex) =>
						runWorker(shardIndex, workerCount),
					),
				)
			: [];
	for (const shardResult of await shardResults) mergeResult(result, shardResult);
	result.failures.sort();

	console.log(
		JSON.stringify(
			{
				reference: reference.reference,
				workers: workerCount,
				cases: corpus.length,
				domains: { exact: result.exactDomains, total: result.domains },
				residues: { exact: result.exactResidues, total: result.residues },
				scores: {
					withinDisplayedPrecision: result.scoresWithinDisplayedPrecision,
					total: result.scoreCount,
					meanAbsoluteDifference:
						result.scoreCount === 0 ? 0 : result.absoluteScoreDifferenceSum / result.scoreCount,
				},
				eValues: {
					withinDisplayedPrecision: result.eValuesWithinDisplayedPrecision,
					total: result.eValueCount,
					meanRelativeDifference:
						result.eValueCount === 0 ? 0 : result.relativeEValueDifferenceSum / result.eValueCount,
					maximumRelativeDifference: result.maximumRelativeEValueDifference,
					maximumRelativeDifferenceCase: result.maximumRelativeEValueDifferenceCase,
					maximumAllowedRelativeDifference: maximumAllowedRelativeEValueDifference,
				},
				maximumAbsoluteBitScoreDifference: result.maximumScoreDifference,
				maximumAllowedAbsoluteBitScoreDifference: maximumAllowedAbsoluteScoreDifference,
				maximumAbsoluteBiasDifference: result.maximumBiasDifference,
				maximumAllowedAbsoluteBiasDifference,
				failures: result.failures,
			},
			null,
			2,
		),
	);

	if (
		result.failures.length > 0 ||
		result.maximumScoreDifference > maximumAllowedAbsoluteScoreDifference ||
		result.maximumBiasDifference > maximumAllowedAbsoluteBiasDifference ||
		result.maximumRelativeEValueDifference > maximumAllowedRelativeEValueDifference
	) {
		process.exitCode = 1;
	}
}
