import { fileURLToPath } from "node:url";
import preact from "@preact/preset-vite";
import tailwindcss from "@tailwindcss/vite";
import { defineConfig } from "vite";

// The demo consumes the built package, exactly as an npm dependent would.
// `tools/build-browser.sh` has to run before `vite build`.
const distribution = fileURLToPath(new URL("../browser/dist", import.meta.url));
const root = fileURLToPath(new URL(".", import.meta.url));

export default defineConfig({
  // Relative URLs keep one build working at any GitHub Pages path, so the
  // repository name is not baked into the output.
  base: "./",
  plugins: [preact(), tailwindcss()],
  resolve: {
    alias: { "@anarcism": distribution },
  },
  server: {
    // `browser/dist` sits outside the Vite root.
    fs: { allow: [root, distribution] },
  },
  build: {
    target: "es2022",
  },
});
