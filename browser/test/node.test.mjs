import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

import init, {
  AnarcismError,
  initSync,
  numberFasta,
  numberSequence,
  numberSequences,
  validateAntibodyPair,
} from "../dist/index.js";

const VH = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";
const VL = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR";
let wasm;

test.before(async () => {
  wasm = await readFile(new URL("../dist/anarcism.wasm", import.meta.url));
  await init(wasm);
});

test("synchronous initialization accepts bytes and compiled modules", () => {
  const bytesApi = initSync(wasm);
  assert.equal(bytesApi.version, "0.1.0");
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
  assert.equal(initSync(module).numberSequence(VL).domains[0].chainType, "K");
});

test("single, batch, FASTA, and germline APIs execute in WebAssembly", () => {
  const heavy = numberSequence(VH, { assignGermline: true });
  assert.equal(heavy.domains[0].chainType, "H");
  assert.equal(heavy.domains[0].end, 120);
  assert.equal(heavy.domains[0].germline.vGene, "IGHV14-4*02");
  assert.equal("confidence" in heavy.domains[0], false);

  const batch = numberSequences([
    { id: "heavy", sequence: VH },
    { id: "light", sequence: VL },
  ]);
  assert.deepEqual(batch.map((result) => result.id), ["heavy", "light"]);
  assert.deepEqual(batch.map((result) => result.domains[0].chainType), ["H", "K"]);

  const fasta = numberFasta(`>heavy\n${VH}\n>light\n${VL}\n`);
  assert.deepEqual(fasta.map((result) => result.id), ["heavy", "light"]);
});

test("pair validation accepts VH/VL and rejects swapped chains", () => {
  const valid = validateAntibodyPair(VH, VL);
  assert.equal(valid.ok, true);
  assert.equal(valid.vh.chainType, "H");
  assert.equal(valid.vl.chainType, "K");

  const swapped = validateAntibodyPair(VL, VH);
  assert.equal(swapped.ok, false);
  assert.ok(swapped.errors.length >= 2);
});

test("invalid input is a structured JavaScript error", () => {
  assert.throws(
    () => numberSequence("ACDX"),
    (error) => error instanceof AnarcismError && error.code === "INVALID_SEQUENCE",
  );
});
