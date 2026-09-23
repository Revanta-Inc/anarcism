import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

import { Anarcism, AnarcismError } from "../dist/index.js";
import { normalizeWorkerCount, parseFasta, recommendedWorkerCount } from "../dist/worker-pool.js";

const VH =
	"EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";
const VL =
	"DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR";
let wasm;
let anarcism;

test.before(async () => {
	wasm = await readFile(new URL("../dist/anarcism.wasm", import.meta.url));
	anarcism = await Anarcism.create({ source: wasm });
});

test("synchronous initialization accepts bytes and compiled modules", () => {
	const bytesApi = Anarcism.createSync(wasm);
	assert.equal(bytesApi.version, "1.0.0");
	assert.deepEqual(bytesApi.chains(), ["H", "K", "L", "A", "B", "G", "D"]);
	assert.deepEqual(bytesApi.species(), [
		"human",
		"mouse",
		"rat",
		"rabbit",
		"rhesus",
		"pig",
		"alpaca",
		"cow",
	]);
	assert.equal(bytesApi.numberSequence(VH).domains[0].chainType, "H");
	const module = new WebAssembly.Module(wasm);
	assert.equal(Anarcism.createSync(module).numberSequence(VL).domains[0].chainType, "K");
});

test("single, batch, FASTA, and germline APIs execute in WebAssembly", () => {
	const heavy = anarcism.numberSequence(VH, { assignGermline: true });
	assert.equal(heavy.domains[0].chainType, "H");
	assert.equal(heavy.domains[0].end, 120);
	assert.equal(heavy.domains[0].germline.vGene, "IGHV14-4*02");
	assert.ok(Number.isFinite(heavy.domains[0].eValue));
	assert.ok(heavy.domains[0].eValue > 0);
	assert.ok(heavy.domains[0].bias >= 0);
	assert.ok(heavy.domains[0].queryStart < heavy.domains[0].queryEnd);
	assert.ok(heavy.domains[0].alternativeHits[0].eValue > 0);
	assert.ok(heavy.domains[0].alternativeHits[0].bias >= 0);
	assert.ok(
		heavy.domains[0].alternativeHits[0].queryStart < heavy.domains[0].alternativeHits[0].queryEnd,
	);
	assert.equal("confidence" in heavy.domains[0], false);

	const batchInputs = Array.from({ length: 8 }, (_, index) => ({
		id: index % 2 === 0 ? `heavy-${index}` : `light-${index}`,
		sequence: index % 2 === 0 ? VH : VL,
	}));
	const batch = anarcism.numberSequences(batchInputs);
	assert.deepEqual(
		batch.map((result) => result.id),
		batchInputs.map(({ id }) => id),
	);
	assert.deepEqual(
		batch.map((result) => result.domains[0].chainType),
		["H", "K", "H", "K", "H", "K", "H", "K"],
	);
	assert.deepEqual(batch[0].domains, anarcism.numberSequence(VH).domains);

	const fasta = anarcism.numberFasta(`>heavy\n${VH}\n>light\n${VL}\n`);
	assert.deepEqual(
		fasta.map((result) => result.id),
		["heavy", "light"],
	);
});

test("pair validation accepts VH/VL and rejects swapped chains", () => {
	const valid = anarcism.validateAntibodyPair(VH, VL);
	assert.equal(valid.ok, true);
	assert.equal(valid.vh.chainType, "H");
	assert.equal(valid.vl.chainType, "K");

	const swapped = anarcism.validateAntibodyPair(VL, VH);
	assert.equal(swapped.ok, false);
	assert.ok(swapped.errors.length >= 2);
});

test("unknown residues are supported and preserved", () => {
	const unknownIndex = 60;
	const sequence = `${VH.slice(0, unknownIndex)}X${VH.slice(unknownIndex + 1)}`;
	const result = anarcism.numberSequence(sequence);
	assert.equal(result.normalizedSequence, sequence);
	assert.equal(
		result.domains[0].numbering.find((residue) => residue.sequenceIndex === unknownIndex).aminoAcid,
		"X",
	);
});

test("invalid input is a structured JavaScript error", () => {
	assert.throws(
		() => anarcism.numberSequence("ACD!"),
		(error) => error instanceof AnarcismError && error.code === "INVALID_SEQUENCE",
	);
});

test("creation loads the bundled module from the package in Node.js", async () => {
	const bundled = await Anarcism.create();
	assert.equal(bundled.numberSequence(VL).domains[0].chainType, "K");
});

test("creation honours abort signals", async () => {
	const aborted = AbortSignal.abort();
	await assert.rejects(
		Anarcism.create({ source: wasm, signal: aborted }),
		(error) => error === aborted.reason,
	);

	const controller = new AbortController();
	const reason = new Error("stop");
	const pending = Anarcism.create({ source: wasm, signal: controller.signal });
	controller.abort(reason);
	await assert.rejects(pending, (error) => error === reason);

	const created = await Anarcism.create({
		source: wasm,
		signal: new AbortController().signal,
	});
	assert.equal(created.numberSequence(VH).domains[0].chainType, "H");
});

test("worker pool parses FASTA and preserves the core input limits", () => {
	assert.deepEqual(parseFasta("\n>one description\r\nACD\r\nEF\r\n>two\nGH\n"), [
		{ id: "one description", sequence: "ACDEF" },
		{ id: "two", sequence: "GH" },
	]);
	assert.throws(
		() => parseFasta(">empty\n\n>next\nACD\n"),
		(error) => error.code === "INVALID_FASTA" && error.inputId === "empty",
	);
	const tooMany = Array.from({ length: 1_001 }, (_, index) => `>s${index}\nA\n`).join("");
	assert.throws(
		() => parseFasta(tooMany),
		(error) => error.code === "BATCH_TOO_LARGE",
	);
});

test("worker pool uses a bounded hardware-aware default", () => {
	assert.equal(recommendedWorkerCount(64), 8);
	assert.equal(recommendedWorkerCount(2), 2);
	assert.equal(normalizeWorkerCount(0), 1);
	assert.equal(normalizeWorkerCount(100), 16);
});
