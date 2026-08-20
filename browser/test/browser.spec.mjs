import { expect, test } from "@playwright/test";

const VH = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";
const VL = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR";

test("initializes the embedded WASM module", async ({ page }) => {
  await page.goto("/");
  await page.waitForFunction(() => globalThis.anarcismStage !== "initializing", null, {
    timeout: 30_000,
  });
  expect(await page.evaluate(() => globalThis.anarcismStage)).toBe("ready");
});

test("numbers a VH/VL pair without network sequence transfer", async ({ page }) => {
  const requests = [];
  page.on("request", (request) => requests.push(request.url()));
  await page.goto("/");
  await page.waitForFunction(() => globalThis.anarcismStage !== "initializing", null, {
    timeout: 30_000,
  });
  expect(await page.evaluate(() => globalThis.anarcismStage)).toBe("ready");
  const heavy = await page.evaluate((vh) => globalThis.anarcismApi.numberSequence(vh), VH);
  const pair = await page.evaluate(
    ([vh, vl]) => globalThis.anarcismApi.validateAntibodyPair(vh, vl),
    [VH, VL],
  );

  expect(heavy.domains[0].chainType).toBe("H");
  expect(heavy.domains[0].numbering[0].position).toBe(1);
  expect(pair.ok).toBe(true);
  expect(requests).toEqual(expect.arrayContaining([
    "http://127.0.0.1:43991/",
    "http://127.0.0.1:43991/index.js",
    "http://127.0.0.1:43991/anarcism.wasm",
  ]));
  expect(requests.some((url) => url.includes(VH.slice(0, 16)))).toBe(false);
});
