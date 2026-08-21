import { createReadStream } from "node:fs";
import { stat } from "node:fs/promises";
import { createServer } from "node:http";
import { extname, resolve, sep } from "node:path";

const root = resolve(import.meta.dirname, "..");
const contentTypes = new Map([
  [".html", "text/html; charset=utf-8"],
  [".js", "text/javascript; charset=utf-8"],
  [".css", "text/css; charset=utf-8"],
  [".wasm", "application/wasm"],
]);

createServer(async (request, response) => {
  const url = new URL(request.url, "http://localhost");
  if (url.pathname === "/") {
    response.writeHead(302, { location: "/demo/" }).end();
    return;
  }
  const relative = url.pathname === "/demo/" ? "demo/index.html" : url.pathname.slice(1);
  const path = resolve(root, relative);
  if (!path.startsWith(`${root}${sep}`)) {
    response.writeHead(403).end();
    return;
  }
  try {
    const metadata = await stat(path);
    if (!metadata.isFile()) throw new Error("not a file");
    response.writeHead(200, {
      "content-length": metadata.size,
      "content-type": contentTypes.get(extname(path)) ?? "application/octet-stream",
    });
    createReadStream(path).pipe(response);
  } catch {
    response.writeHead(404).end();
  }
}).listen(8080, "127.0.0.1", () => {
  console.log("anarcism demo: http://127.0.0.1:8080/");
});
