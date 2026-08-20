import { spawn } from "node:child_process";
import { setTimeout as delay } from "node:timers/promises";

import { chromium } from "playwright";

const origin = "http://127.0.0.1:43991";
const server = spawn(process.execPath, ["test/server.mjs"], { stdio: "ignore" });
let browser;

try {
  for (let attempt = 0; attempt < 50; attempt += 1) {
    try {
      const response = await fetch(origin);
      if (response.ok) break;
    } catch {
      if (attempt === 49) throw new Error("benchmark server did not start");
    }
    await delay(100);
  }

  browser = await chromium.launch();
  const page = await browser.newPage();
  await page.goto(origin);
  await page.waitForFunction(() => globalThis.anarcismStage === "ready");
  const result = await page.evaluate(async ([vh, vl]) => {
    const moduleBytes = await (await fetch("/anarcism.wasm")).arrayBuffer();
    const module = await WebAssembly.compile(moduleBytes);
    let started = performance.now();
    await globalThis.anarcismApi.default(module);
    const cachedInitializationMs = performance.now() - started;

    for (let index = 0; index < 3; index += 1) {
      globalThis.anarcismApi.numberSequence(vh);
      globalThis.anarcismApi.numberSequence(vl);
    }
    const measure = (sequence) => {
      const samples = [];
      for (let index = 0; index < 25; index += 1) {
        started = performance.now();
        globalThis.anarcismApi.numberSequence(sequence);
        samples.push(performance.now() - started);
      }
      samples.sort((left, right) => left - right);
      return {
        medianMs: samples[Math.floor(samples.length / 2)],
        p95Ms: samples[Math.floor(samples.length * 0.95)],
      };
    };

    const inputs = Array.from({ length: 100 }, (_, pairIndex) => [
      { id: `vh-${pairIndex}`, sequence: vh },
      { id: `vl-${pairIndex}`, sequence: vl },
    ]).flat();
    started = performance.now();
    globalThis.anarcismApi.numberSequences(inputs);
    const hundredPairsMs = performance.now() - started;
    return {
      userAgent: navigator.userAgent,
      coldInitializationMs: globalThis.anarcismInitializationMs,
      cachedInitializationMs,
      vh: measure(vh),
      vl: measure(vl),
      hundredPairsMs,
    };
  }, [
    "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA",
    "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR",
  ]);
  console.log(JSON.stringify(result, null, 2));
} finally {
  await browser?.close();
  server.kill();
}
