import { createReadStream } from "node:fs";
import { stat } from "node:fs/promises";
import { createServer } from "node:http";
import { extname, resolve } from "node:path";

const browserRoot = resolve(import.meta.dirname, "..");
const repositoryRoot = resolve(browserRoot, "..");
const routes = new Map([
  ["/", resolve(import.meta.dirname, "fixture.html")],
  ["/index.js", resolve(browserRoot, "dist/index.js")],
  ["/anarcism.wasm", resolve(browserRoot, "dist/anarcism.wasm")],
  ["/demo/", resolve(repositoryRoot, "demo/index.html")],
  ["/browser/dist/index.js", resolve(browserRoot, "dist/index.js")],
  ["/browser/dist/worker-pool.js", resolve(browserRoot, "dist/worker-pool.js")],
  ["/browser/dist/worker.js", resolve(browserRoot, "dist/worker.js")],
  ["/browser/dist/anarcism.wasm", resolve(browserRoot, "dist/anarcism.wasm")],
]);
const contentTypes = new Map([
  [".html", "text/html; charset=utf-8"],
  [".js", "text/javascript; charset=utf-8"],
  [".wasm", "application/wasm"],
]);

createServer(async (request, response) => {
  const path = routes.get(new URL(request.url, "http://localhost").pathname);
  if (!path) {
    response.writeHead(404).end();
    return;
  }
  try {
    const metadata = await stat(path);
    response.writeHead(200, {
      "content-length": metadata.size,
      "content-type": contentTypes.get(extname(path)) ?? "application/octet-stream",
      "cache-control": "no-store",
    });
    createReadStream(path).pipe(response);
  } catch {
    response.writeHead(500).end();
  }
}).listen(43991, "127.0.0.1");
