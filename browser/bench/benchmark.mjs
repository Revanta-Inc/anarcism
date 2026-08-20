import { readFile } from "node:fs/promises";
import { performance } from "node:perf_hooks";

import init, { numberSequence, numberSequences } from "../dist/index.js";

const VH = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";
const VL = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR";
const wasm = await readFile(new URL("../dist/anarcism.wasm", import.meta.url));

let started = performance.now();
await init(wasm);
const coldInitializationMs = performance.now() - started;

const module = await WebAssembly.compile(wasm);
started = performance.now();
await init(module);
const cachedInitializationMs = performance.now() - started;

for (let index = 0; index < 3; index += 1) {
  numberSequence(VH);
  numberSequence(VL);
}

function distribution(sequence, iterations = 25) {
  const samples = [];
  for (let index = 0; index < iterations; index += 1) {
    const before = performance.now();
    numberSequence(sequence);
    samples.push(performance.now() - before);
  }
  samples.sort((left, right) => left - right);
  return {
    medianMs: samples[Math.floor(samples.length / 2)],
    p95Ms: samples[Math.floor(samples.length * 0.95)],
  };
}

const inputs = Array.from({ length: 100 }, (_, pairIndex) => [
  { id: `vh-${pairIndex}`, sequence: VH },
  { id: `vl-${pairIndex}`, sequence: VL },
]).flat();
started = performance.now();
numberSequences(inputs);
const hundredPairsMs = performance.now() - started;

console.log(JSON.stringify({
  runtime: process.version,
  coldInitializationMs,
  cachedInitializationMs,
  vh: distribution(VH),
  vl: distribution(VL),
  hundredPairsMs,
}, null, 2));
