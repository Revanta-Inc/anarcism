import { spawn } from "node:child_process";
import { setTimeout as delay } from "node:timers/promises";

import { chromium } from "playwright";

const origin = "http://127.0.0.1:43991";
const warmupIterations = setting("ANARCISM_BENCH_WARMUP", 3);
const singleIterations = setting("ANARCISM_BENCH_SINGLE_ITERATIONS", 25);
const pairCount = setting("ANARCISM_BENCH_PAIR_COUNT", 100);
const poolWorkerCount = setting("ANARCISM_BENCH_BROWSER_WORKERS", 4);
const server = spawn(process.execPath, ["test/server.mjs"], { stdio: "ignore" });
let browser;

try {
  let serverReady = false;
  for (let attempt = 0; attempt < 50; attempt += 1) {
    try {
      const response = await fetch(origin);
      if (response.ok) {
        serverReady = true;
        break;
      }
    } catch {}
    await delay(100);
  }
  if (!serverReady) throw new Error("benchmark server did not start");

  browser = await chromium.launch();
  const page = await browser.newPage();
  await page.goto(origin);
  await page.waitForFunction(() => globalThis.anarcismStage === "ready");
  const result = await page.evaluate(async ([vh, vl, config]) => {
    const moduleBytes = await (await fetch("/anarcism.wasm")).arrayBuffer();
    const module = await WebAssembly.compile(moduleBytes);
    let started = performance.now();
    await globalThis.anarcismApi.default(module);
    const cachedInitializationMs = performance.now() - started;

    for (let index = 0; index < config.warmupIterations; index += 1) {
      assertSingleDomain(globalThis.anarcismApi.numberSequence(vh), "H");
      assertSingleDomain(globalThis.anarcismApi.numberSequence(vl), "K");
    }
    const measure = (sequence, expectedChain) => {
      const samples = [];
      for (let index = 0; index < config.singleIterations; index += 1) {
        started = performance.now();
        const sequenceResult = globalThis.anarcismApi.numberSequence(sequence);
        const elapsedMs = performance.now() - started;
        assertSingleDomain(sequenceResult, expectedChain);
        samples.push(elapsedMs);
      }
      samples.sort((left, right) => left - right);
      return {
        medianMs: samples[Math.floor(samples.length / 2)],
        p95Ms: samples[Math.floor(samples.length * 0.95)],
      };
    };

    const inputs = Array.from({ length: config.pairCount }, (_, pairIndex) => [
      { id: `vh-${pairIndex}`, sequence: vh },
      { id: `vl-${pairIndex}`, sequence: vl },
    ]).flat();
    started = performance.now();
    const batch = globalThis.anarcismApi.numberSequences(inputs);
    const batchMs = performance.now() - started;
    assertBatch(batch, inputs.length);

    const { AnarcismWorkerPool } = await import("/browser/dist/worker-pool.js");
    const pool = new AnarcismWorkerPool();
    let workerPool;
    try {
      started = performance.now();
      await pool.initialize(config.poolWorkerCount);
      const poolInitializationMs = performance.now() - started;
      const poolWarmupInputs = inputs.slice(0, Math.min(inputs.length, config.poolWorkerCount * 2));
      started = performance.now();
      await pool.numberSequences(poolWarmupInputs);
      const poolWarmupMs = performance.now() - started;
      started = performance.now();
      const pooledBatch = await pool.numberSequences(inputs);
      const pooledBatchMs = performance.now() - started;
      assertBatch(pooledBatch, inputs.length);
      workerPool = {
        workers: pool.size,
        initializationMs: poolInitializationMs,
        warmupSequenceCount: poolWarmupInputs.length,
        warmupMs: poolWarmupMs,
        totalMs: pooledBatchMs,
        perSequenceMs: pooledBatchMs / inputs.length,
        sequencesPerSecond: inputs.length * 1_000 / pooledBatchMs,
        speedup: batchMs / pooledBatchMs,
      };
    } finally {
      pool.terminate();
    }
    return {
      implementation: "Browser WASM",
      userAgent: navigator.userAgent,
      config: {
        ...config,
        sequenceCount: inputs.length,
        threads: 1,
      },
      coldInitializationMs: globalThis.anarcismInitializationMs,
      cachedInitializationMs,
      vh: measure(vh, "H"),
      vl: measure(vl, "K"),
      batch: {
        totalMs: batchMs,
        perSequenceMs: batchMs / inputs.length,
        sequencesPerSecond: inputs.length * 1_000 / batchMs,
      },
      workerPool,
    };

    function assertSingleDomain(sequenceResult, expectedChain) {
      if (sequenceResult.domains.length !== 1 || sequenceResult.domains[0].chainType !== expectedChain) {
        throw new Error(`WASM returned the wrong domain for ${expectedChain}`);
      }
    }

    function assertBatch(batchResults, expectedLength) {
      if (batchResults.length !== expectedLength) {
        throw new Error("WASM returned an incomplete batch");
      }
      batchResults.forEach((sequenceResult, index) => {
        assertSingleDomain(sequenceResult, index % 2 === 0 ? "H" : "K");
      });
    }
  }, [
    "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA",
    "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR",
    { warmupIterations, singleIterations, pairCount, poolWorkerCount },
  ]);
  console.log(JSON.stringify(result, null, 2));
} finally {
  await browser?.close();
  server.kill();
}

function setting(name, fallback) {
  const value = Number.parseInt(process.env[name] ?? "", 10);
  return Number.isInteger(value) && value > 0 ? value : fallback;
}
