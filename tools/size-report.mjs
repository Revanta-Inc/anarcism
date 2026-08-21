import { readFileSync } from "node:fs";
import { brotliCompressSync, constants, gzipSync } from "node:zlib";
import { resolve } from "node:path";

const root = resolve(import.meta.dirname, "..");
const files = [
  "browser/dist/anarcism.wasm",
  "browser/dist/index.js",
  "browser/dist/index.d.ts",
  "browser/dist/worker-pool.js",
  "browser/dist/worker-pool.d.ts",
  "browser/dist/worker.js",
  "browser/dist/THIRD_PARTY_NOTICES.md",
  "browser/package.json",
];
const limit = 500_000;

function sizes(path) {
  const bytes = readFileSync(resolve(root, path));
  return {
    raw: bytes.byteLength,
    gzip: gzipSync(bytes, { level: 9 }).byteLength,
    brotli: brotliCompressSync(bytes, {
      params: { [constants.BROTLI_PARAM_QUALITY]: 11 },
    }).byteLength,
  };
}

const rows = files.map((path) => ({ path, ...sizes(path) }));
const total = rows.reduce(
  (sum, row) => ({
    raw: sum.raw + row.raw,
    gzip: sum.gzip + row.gzip,
    brotli: sum.brotli + row.brotli,
  }),
  { raw: 0, gzip: 0, brotli: 0 },
);

console.log("file\traw\tgzip\tbrotli");
for (const row of rows) {
  console.log(`${row.path}\t${row.raw}\t${row.gzip}\t${row.brotli}`);
}
console.log(`complete distribution\t${total.raw}\t${total.gzip}\t${total.brotli}`);

for (const asset of ["assets/profiles.bin", "assets/germlines.bin"]) {
  const row = sizes(asset);
  console.log(`${asset} (embedded)\t${row.raw}\t${row.gzip}\t${row.brotli}`);
}

const rustRelease = "target/wasm32-unknown-unknown/release/anarcism_wasm.wasm";
try {
  const row = sizes(rustRelease);
  console.log(`${rustRelease} (Rust release)\t${row.raw}\t${row.gzip}\t${row.brotli}`);
} catch {
  // This file is optional when inspecting an unpacked npm package.
}

if (total.gzip >= limit || total.brotli >= limit) {
  console.error(`compressed distribution exceeds the ${limit}-byte budget`);
  process.exitCode = 1;
}
