// The worker pool is the one part of the package that Node cannot cover: it
// needs the DOM `Worker` constructor with `{ type: "module" }`. These tests
// drive it directly and assert nothing about any user interface.
import { expect, test } from "@playwright/test";

const VH =
	"EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";
const VL =
	"DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR";
const SCFV = `${VH}GGGGSGGGGSGGGGS${VL}`;

async function openFixture(page) {
	await page.goto("/");
	await page.waitForFunction(() => Boolean(globalThis.anarcismWorkerPool), null, {
		timeout: 30_000,
	});
}

test("initializes one embedded engine per module worker", async ({ page }) => {
	await openFixture(page);

	const metadata = await page.evaluate(() => {
		globalThis.pool = new globalThis.anarcismWorkerPool.AnalysisWorkerPool();
		return globalThis.pool.initialize(2);
	});

	expect(metadata.workerCount).toBe(2);
	expect(metadata.version).toBe("0.1.0");
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

	// Analysis must run off the page thread, in dedicated workers.
	await expect.poll(() => page.workers().length).toBe(2);
	for (const worker of page.workers()) {
		expect(new URL(worker.url()).pathname).toBe("/worker.js");
	}
});

test("grows, shrinks, and terminates the pool", async ({ page }) => {
	await openFixture(page);
	await page.evaluate(async () => {
		globalThis.pool = new globalThis.anarcismWorkerPool.AnalysisWorkerPool();
		await globalThis.pool.initialize(1);
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

	// Mixed lengths so the length-balancing shards are uneven and the
	// reassembled order is not merely the order each worker happened to finish.
	const analysis = await page.evaluate(
		async ([vh, vl, scfv]) => {
			const pool = new globalThis.anarcismWorkerPool.AnalysisWorkerPool();
			await pool.initialize(2);
			try {
				const batch = await pool.numberSequences([
					{ id: "a", sequence: vl },
					{ id: "b", sequence: scfv },
					{ id: "c", sequence: vh },
					{ id: "d", sequence: vl },
					{ id: "e", sequence: scfv },
				]);
				// `assignGermline` confirms options survive the trip into the worker and
				// that the larger result shape survives structured cloning back out.
				const single = await pool.numberSequence(vh, { assignGermline: true });
				const fasta = await pool.numberFasta(`>heavy\n${vh}\n>light\n${vl}\n`);
				return {
					ids: batch.map((result) => result.id),
					chains: batch.map((result) => result.domains.map((domain) => domain.chainType).join("")),
					single: {
						chainType: single.domains[0].chainType,
						firstPosition: single.domains[0].numbering[0].position,
						vGene: single.domains[0].germline.vGene,
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
	expect(analysis.fastaIds).toEqual(["heavy", "light"]);
});

test("surfaces worker failures as structured errors", async ({ page }) => {
	await openFixture(page);

	const failures = await page.evaluate(
		async ([vh]) => {
			const { AnalysisWorkerPool } = globalThis.anarcismWorkerPool;
			const capture = async (run) => {
				try {
					await run();
					return "resolved";
				} catch (error) {
					return { name: error.name, code: error.code };
				}
			};

			const uninitialized = new AnalysisWorkerPool();
			const pool = new AnalysisWorkerPool();
			await pool.initialize(1);
			const invalidSequence = await capture(() => pool.numberSequence("ACDX"));
			// A rejected request must not poison the worker for the next one.
			const recovered = (await pool.numberSequence(vh)).domains[0].chainType;
			pool.terminate();

			return {
				beforeInitialize: await capture(() => uninitialized.numberSequence(vh)),
				invalidSequence,
				recovered,
				afterTerminate: await capture(() => pool.numberSequence(vh)),
			};
		},
		[VH],
	);

	expect(failures.beforeInitialize).toEqual({
		name: "AnarcismError",
		code: "NOT_INITIALIZED",
	});
	expect(failures.invalidSequence).toEqual({
		name: "AnarcismError",
		code: "INVALID_SEQUENCE",
	});
	expect(failures.recovered).toBe("H");
	expect(failures.afterTerminate).toEqual({
		name: "AnarcismError",
		code: "NOT_INITIALIZED",
	});
});

test("numbers without sending a sequence anywhere", async ({ page }) => {
	const requested = [];
	page.on("request", (request) => requested.push(request.url()));
	await openFixture(page);

	await page.evaluate(
		async ([vh]) => {
			const pool = new globalThis.anarcismWorkerPool.AnalysisWorkerPool();
			await pool.initialize(1);
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
