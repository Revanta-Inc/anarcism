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

test("demo initializes and renders an IMGT-numbered domain", async ({ page }) => {
  await page.goto("/demo/");

  await expect(page.locator("#run")).toBeEnabled({ timeout: 30_000 });
  expect(page.workers().map((worker) => worker.url())).toEqual([
    "http://127.0.0.1:43991/demo/worker.js",
  ]);
  await expect(page.locator("#ver")).toHaveText("0.1.0");
  await expect(page.locator("#chains input")).toHaveCount(7);
  await expect(page.locator("#species option")).toHaveCount(8);

  await page.locator("#run").click();
  await expect(page.locator(".badge.chain")).toHaveText("H");
  await expect(page.locator(".res .p").first()).toHaveText("1");
  await expect(page.locator("details")).toHaveCount(0);

  await page.locator("#alts").check();
  await page.locator("#run").click();
  await expect(page.locator("details summary")).toHaveText("3 alternative hits");

  await page.locator('[data-ex="lys"]').click();
  await expect(page.locator("#out .none")).toHaveText("No variable domain detected.");
});
