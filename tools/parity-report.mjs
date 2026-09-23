import { readFile } from "node:fs/promises";
import { availableParallelism } from "node:os";
import { Worker, isMainThread, parentPort, workerData } from "node:worker_threads";

import init, { numberSequence, validateAntibodyPair } from "../browser/dist/index.js";

const originalCorpusUrl = new URL("../tests/golden/corpus.json", import.meta.url);
const expandedCorpusUrl = new URL("../tests/golden/corpus_v2.json", import.meta.url);
const expandedReferenceUrl = new URL("../tests/golden/corpus_v2_reference.json", import.meta.url);
const wasmUrl = new URL("../browser/dist/anarcism.wasm", import.meta.url);
const textEncoder = new TextEncoder();

function emptyResult() {
	return {
		domains: 0,
		residues: 0,
		exactDomains: 0,
		exactResidues: 0,
		maximumScoreDifference: 0,
		failures: [],
	};
}

function mergeResult(target, source) {
	target.domains += source.domains;
	target.residues += source.residues;
	target.exactDomains += source.exactDomains;
	target.exactResidues += source.exactResidues;
	target.maximumScoreDifference = Math.max(
		target.maximumScoreDifference,
		source.maximumScoreDifference,
	);
	target.failures.push(...source.failures);
}

function fnv1a64(value) {
	let hash = 0xcbf29ce484222325n;
	for (const byte of textEncoder.encode(value)) {
		hash ^= BigInt(byte);
		hash = BigInt.asUintN(64, hash * 0x100000001b3n);
	}
	return hash.toString(16).padStart(16, "0");
}

function numberingFingerprint(numbering) {
	return fnv1a64(
		numbering
			.map((residue) =>
				[residue.sequenceIndex, residue.aminoAcid, residue.position, residue.insertionCode].join(
					"|",
				),
			)
			.join("\n") + (numbering.length > 0 ? "\n" : ""),
	);
}

async function readJson(url) {
	return JSON.parse(await readFile(url, "utf8"));
}

async function initializeWasm() {
	await init(await readFile(wasmUrl));
}

async function checkOriginalCorpus(corpus) {
	const result = emptyResult();
	await initializeWasm();

	for (const testCase of corpus.cases) {
		const observed = numberSequence(testCase.sequence);
		if (observed.domains.length !== testCase.referenceDomains.length) {
			result.failures.push(`${testCase.id}: domain count`);
			continue;
		}
		for (let index = 0; index < observed.domains.length; index += 1) {
			result.domains += 1;
			const actual = observed.domains[index];
			const expected = testCase.referenceDomains[index];
			const domainExact =
				actual.chainType === expected.chainType &&
				actual.species === expected.species &&
				actual.start === expected.start &&
				actual.end === expected.end &&
				actual.paddedImgtAlignment === expected.paddedImgtAlignment;
			if (domainExact) result.exactDomains += 1;
			else result.failures.push(`${testCase.id}: domain ${index}`);
			result.maximumScoreDifference = Math.max(
				result.maximumScoreDifference,
				Math.abs(actual.bitScore - expected.bitScore),
			);
			result.residues += expected.numbering.length;
			for (let residueIndex = 0; residueIndex < expected.numbering.length; residueIndex += 1) {
				const left = actual.numbering[residueIndex];
				const right = expected.numbering[residueIndex];
				if (
					left &&
					left.sequenceIndex === right.sequenceIndex &&
					left.aminoAcid === right.aminoAcid &&
					left.position === right.position &&
					left.insertionCode === right.insertionCode
				) {
					result.exactResidues += 1;
				} else {
					result.failures.push(`${testCase.id}: residue ${residueIndex}`);
				}
			}
		}
	}

	for (const pair of corpus.pairs) {
		const observed = validateAntibodyPair(pair.vh, pair.vl);
		if (observed.ok !== pair.ok) result.failures.push(`${pair.id}: pair outcome`);
	}
	return result;
}

async function checkExpandedShard(shardIndex, shardCount) {
	const expandedCorpusPromise = readJson(expandedCorpusUrl);
	const expandedReferencePromise = readJson(expandedReferenceUrl);
	await initializeWasm();
	const [expandedCorpus, expandedReference] = await Promise.all([
		expandedCorpusPromise,
		expandedReferencePromise,
	]);
	const result = emptyResult();

	for (let caseIndex = shardIndex; caseIndex < expandedCorpus.length; caseIndex += shardCount) {
		const testCase = expandedCorpus[caseIndex];
		const referenceCase = expandedReference.cases[caseIndex];
		if (
			!referenceCase ||
			testCase.id !== referenceCase.id ||
			testCase.category !== referenceCase.category
		) {
			result.failures.push(`${testCase.id}: expanded reference identity`);
			continue;
		}
		const observed = numberSequence(testCase.seq);
		if (observed.domains.length !== referenceCase.domains.length) {
			result.failures.push(`${testCase.id}: expanded domain count`);
			continue;
		}
		for (let index = 0; index < observed.domains.length; index += 1) {
			result.domains += 1;
			const actual = observed.domains[index];
			const expected = referenceCase.domains[index];
			const numberingExact =
				actual.numbering.length === expected.numberingLength &&
				numberingFingerprint(actual.numbering) === expected.numberingFnv1a64;
			const domainExact =
				actual.chainType === expected.chainType &&
				actual.species === expected.species &&
				actual.start === expected.start &&
				actual.end === expected.end &&
				fnv1a64(actual.paddedImgtAlignment) === expected.paddedAlignmentFnv1a64;
			if (domainExact) result.exactDomains += 1;
			else result.failures.push(`${testCase.id}: expanded domain ${index}`);
			result.maximumScoreDifference = Math.max(
				result.maximumScoreDifference,
				Math.abs(actual.bitScore - expected.bitScore),
			);
			result.residues += expected.numberingLength;
			if (numberingExact) result.exactResidues += expected.numberingLength;
			else result.failures.push(`${testCase.id}: expanded numbering ${index}`);
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
	const result = await checkExpandedShard(workerData.shardIndex, workerData.shardCount);
	parentPort.postMessage(result);
} else {
	const [corpus, expandedCorpus, expandedReference] = await Promise.all([
		readJson(originalCorpusUrl),
		readJson(expandedCorpusUrl),
		readJson(expandedReferenceUrl),
	]);
	const result = emptyResult();
	if (expandedCorpus.length !== expandedReference.cases.length) {
		result.failures.push("expanded corpus/reference case count");
	}

	const workerCount = configuredWorkerCount(expandedCorpus.length);
	const expandedPromise =
		expandedCorpus.length === expandedReference.cases.length
			? Promise.all(
					Array.from({ length: workerCount }, (_, shardIndex) =>
						runWorker(shardIndex, workerCount),
					),
				)
			: Promise.resolve([]);
	const [originalResult, expandedResults] = await Promise.all([
		checkOriginalCorpus(corpus),
		expandedPromise,
	]);
	mergeResult(result, originalResult);
	for (const shardResult of expandedResults) mergeResult(result, shardResult);
	result.failures.sort();

	console.log(
		JSON.stringify(
			{
				reference: [corpus.reference, expandedReference.reference],
				workers: workerCount,
				cases: {
					original: corpus.cases.length,
					expanded: expandedCorpus.length,
					total: corpus.cases.length + expandedCorpus.length,
				},
				domains: { exact: result.exactDomains, total: result.domains },
				residues: { exact: result.exactResidues, total: result.residues },
				maximumAbsoluteBitScoreDifference: result.maximumScoreDifference,
				pairCases: corpus.pairs.length,
				failures: result.failures,
			},
			null,
			2,
		),
	);

	if (result.failures.length > 0 || result.maximumScoreDifference > 4) process.exitCode = 1;
}
