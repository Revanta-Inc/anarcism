import { createReadStream } from "node:fs";
import { stat } from "node:fs/promises";
import { createServer } from "node:http";
import { extname, resolve } from "node:path";

const browserRoot = resolve(import.meta.dirname, "..");
// The package is served flat, so `worker-pool.js` resolves `./worker.js`, and
// the worker resolves `./index.js` and `./anarcism.wasm`, exactly as they do
// from an installed `@revanta/anarcism`.
const routes = new Map([
	["/", resolve(import.meta.dirname, "fixture.html")],
	["/index.js", resolve(browserRoot, "dist/index.js")],
	["/worker-pool.js", resolve(browserRoot, "dist/worker-pool.js")],
	["/worker.js", resolve(browserRoot, "dist/worker.js")],
	["/anarcism.wasm", resolve(browserRoot, "dist/anarcism.wasm")],
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
