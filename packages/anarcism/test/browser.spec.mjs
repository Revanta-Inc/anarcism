import { readFile } from "node:fs/promises";

import { expect, test } from "@playwright/test";

const VH =
	"EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";
const VL =
	"DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR";
const SCFV = `${VH}GGGGSGGGGSGGGGS${VL}`;
const { version } = JSON.parse(await readFile(new URL("../package.json", import.meta.url), "utf8"));

async function openFixture(page) {
	await page.goto("/");
	await page.waitForFunction(() => Boolean(globalThis.anarcismWorkerPool), null, {
		timeout: 30_000,
	});
}

test("initializes one embedded engine per module worker", async ({ page }) => {
	await openFixture(page);

	const metadata = await page.evaluate(async () => {
		const pool = await globalThis.anarcismWorkerPool.AnarcismWorkerPool.create({ workers: 2 });
		globalThis.pool = pool;
		return {
			workerCount: pool.size,
			version: pool.version,
			chains: pool.chains(),
			species: pool.species(),
		};
	});

	expect(metadata.workerCount).toBe(2);
	expect(metadata.version).toBe(version);
	expect(metadata.chains).toEqual(["H", "K", "L", "A", "B", "G", "D"]);
	expect(metadata.species).toEqual([
		"human",
		"mouse",
		"rat",
		"rabbit",
		"rhesus",
		"pig",
		"alpaca",
		"cow",
	]);

	await expect.poll(() => page.workers().length).toBe(2);
	for (const worker of page.workers()) {
		expect(new URL(worker.url()).pathname).toBe("/worker.js");
	}
});

test("grows, shrinks, and terminates the pool", async ({ page }) => {
	await openFixture(page);
	await page.evaluate(async () => {
		globalThis.pool = await globalThis.anarcismWorkerPool.AnarcismWorkerPool.create({
			workers: 1,
		});
	});
	await expect.poll(() => page.workers().length).toBe(1);

	const grown = await page.evaluate(() => globalThis.pool.resize(3));
	expect(grown.workerCount).toBe(3);
	await expect.poll(() => page.workers().length).toBe(3);

	const shrunk = await page.evaluate(() => globalThis.pool.resize(1));
	expect(shrunk.workerCount).toBe(1);
	await expect.poll(() => page.workers().length).toBe(1);

	await page.evaluate(() => globalThis.pool.terminate());
	await expect.poll(() => page.workers().length).toBe(0);
	expect(await page.evaluate(() => globalThis.pool.size)).toBe(0);
});

test("shards a batch across workers and restores input order", async ({ page }) => {
	await openFixture(page);

	const analysis = await page.evaluate(
		async ([vh, vl, scfv]) => {
			const pool = await globalThis.anarcismWorkerPool.AnarcismWorkerPool.create({ workers: 2 });
			try {
				const batch = await pool.numberSequences([
					{ id: "a", sequence: vl },
					{ id: "b", sequence: scfv },
					{ id: "c", sequence: vh },
					{ id: "d", sequence: vl },
					{ id: "e", sequence: scfv },
				]);
				const single = await pool.numberSequence(vh, { assignGermline: true });
				const unknownIndex = 60;
				const unknownSequence = `${vh.slice(0, unknownIndex)}X${vh.slice(unknownIndex + 1)}`;
				const unknown = await pool.numberSequence(unknownSequence);
				const fasta = await pool.numberFasta(`>heavy\n${vh}\n>light\n${vl}\n`);
				return {
					ids: batch.map((result) => result.id),
					chains: batch.map((result) => result.domains.map((domain) => domain.chainType).join("")),
					single: {
						chainType: single.domains[0].chainType,
						firstPosition: single.domains[0].numbering[0].position,
						vGene: single.domains[0].germline.vGene,
					},
					unknown: {
						normalizedSequence: unknown.normalizedSequence,
						aminoAcid: unknown.domains[0].numbering.find(
							(residue) => residue.sequenceIndex === unknownIndex,
						).aminoAcid,
					},
					fastaIds: fasta.map((result) => result.id),
				};
			} finally {
				pool.terminate();
			}
		},
		[VH, VL, SCFV],
	);

	expect(analysis.ids).toEqual(["a", "b", "c", "d", "e"]);
	expect(analysis.chains).toEqual(["K", "HK", "H", "K", "HK"]);
	expect(analysis.single).toEqual({
		chainType: "H",
		firstPosition: 1,
		vGene: "IGHV14-4*02",
	});
	expect(analysis.unknown.aminoAcid).toBe("X");
	expect(analysis.unknown.normalizedSequence[60]).toBe("X");
	expect(analysis.fastaIds).toEqual(["heavy", "light"]);
});

test("aborts analysis and keeps the pool usable", async ({ page }) => {
	await openFixture(page);

	const outcome = await page.evaluate(
		async ([vh, scfv]) => {
			const { AnarcismWorkerPool } = globalThis.anarcismWorkerPool;
			const capture = async (run) => {
				try {
					await run();
					return "resolved";
				} catch (error) {
					return { name: error.name, message: error.message };
				}
			};

			const preAborted = await capture(() =>
				AnarcismWorkerPool.create({ workers: 1, signal: AbortSignal.abort() }),
			);
			const pool = await AnarcismWorkerPool.create({ workers: 2 });
			try {
				const alreadyAborted = await capture(() =>
					pool.numberSequence(vh, { signal: AbortSignal.abort() }),
				);

				const inputs = Array.from({ length: 1_000 }, (_, index) => ({
					id: `s${index}`,
					sequence: scfv,
				}));
				const controller = new AbortController();
				const started = performance.now();
				const batch = pool.numberSequences(inputs, { signal: controller.signal });
				setTimeout(() => controller.abort(new Error("stopped by test")), 20);
				const midBatch = await capture(() => batch);
				const abortLatencyMs = performance.now() - started;

				const recovered = await pool.numberSequences([
					{ id: "heavy", sequence: vh },
					{ id: "scfv", sequence: scfv },
				]);
				return {
					preAborted,
					alreadyAborted,
					midBatch,
					abortLatencyMs,
					busyAfterAbort: pool.busy,
					size: pool.size,
					recovered: recovered.map((result) => result.domains.length),
				};
			} finally {
				pool.terminate();
			}
		},
		[VH, SCFV],
	);

	expect(outcome.preAborted.name).toBe("AbortError");
	expect(outcome.alreadyAborted.name).toBe("AbortError");
	expect(outcome.midBatch).toEqual({ name: "Error", message: "stopped by test" });
	expect(outcome.abortLatencyMs).toBeLessThan(500);
	expect(outcome.busyAfterAbort).toBe(false);
	expect(outcome.size).toBe(2);
	expect(outcome.recovered).toEqual([1, 2]);
	await expect.poll(() => page.workers().length).toBe(0);
});

test("validates a VH/VL pair in the worker", async ({ page }) => {
	await openFixture(page);

	const validation = await page.evaluate(
		async ([vh, vl]) => {
			const pool = await globalThis.anarcismWorkerPool.AnarcismWorkerPool.create({ workers: 1 });
			try {
				const valid = await pool.validateAntibodyPair(vh, vl);
				const swapped = await pool.validateAntibodyPair(vl, vh);
				return {
					valid: {
						ok: valid.ok,
						vh: valid.vh?.chainType,
						vl: valid.vl?.chainType,
						errors: valid.errors.length,
					},
					swapped: { ok: swapped.ok, errors: swapped.errors.length },
				};
			} finally {
				pool.terminate();
			}
		},
		[VH, VL],
	);

	expect(validation.valid).toEqual({ ok: true, vh: "H", vl: "K", errors: 0 });
	expect(validation.swapped.ok).toBe(false);
	expect(validation.swapped.errors).toBeGreaterThanOrEqual(2);
});

test("surfaces worker failures as structured errors", async ({ page }) => {
	await openFixture(page);

	const failures = await page.evaluate(
		async ([vh]) => {
			const { AnarcismWorkerPool } = globalThis.anarcismWorkerPool;
			const capture = async (run) => {
				try {
					await run();
					return "resolved";
				} catch (error) {
					return { name: error.name, code: error.code };
				}
			};

			const pool = await AnarcismWorkerPool.create({ workers: 1 });
			const invalidSequence = await capture(() => pool.numberSequence("ACD!"));
			const recovered = (await pool.numberSequence(vh)).domains[0].chainType;
			pool.terminate();

			return {
				invalidSequence,
				recovered,
				afterTerminate: await capture(() => pool.numberSequence(vh)),
			};
		},
		[VH],
	);

	expect(failures.invalidSequence).toEqual({
		name: "AnarcismError",
		code: "INVALID_SEQUENCE",
	});
	expect(failures.recovered).toBe("H");
	expect(failures.afterTerminate).toEqual({
		name: "AnarcismError",
		code: "TERMINATED",
	});
});

test("numbers without sending a sequence anywhere", async ({ page }) => {
	const requested = [];
	page.on("request", (request) => requested.push(request.url()));
	await openFixture(page);

	await page.evaluate(
		async ([vh]) => {
			const pool = await globalThis.anarcismWorkerPool.AnarcismWorkerPool.create({ workers: 1 });
			try {
				await pool.numberSequence(vh);
			} finally {
				pool.terminate();
			}
		},
		[VH],
	);

	const origin = new URL(page.url()).origin;
	const paths = requested.map((url) => new URL(url).pathname);
	expect(requested.every((url) => new URL(url).origin === origin)).toBe(true);
	expect(
		paths.every((path) =>
			["/", "/index.js", "/worker-pool.js", "/worker.js", "/anarcism.wasm"].includes(path),
		),
	).toBe(true);
	expect(requested.some((url) => url.includes(VH.slice(0, 16)))).toBe(false);
});
