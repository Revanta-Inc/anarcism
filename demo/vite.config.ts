import { fileURLToPath } from "node:url";
import preact from "@preact/preset-vite";
import tailwindcss from "@tailwindcss/vite";
import { defineConfig } from "vite";

const distribution = fileURLToPath(new URL("../packages/anarcism/dist", import.meta.url));
const repositoryRoot = fileURLToPath(new URL("..", import.meta.url));

export default defineConfig({
	// Keep the build independent of its deployment path.
	base: "./",
	plugins: [preact(), tailwindcss()],
	resolve: {
		alias: { "@anarcism": distribution },
	},
	server: {
		fs: { allow: [repositoryRoot] },
	},
	build: {
		target: "es2022",
	},
});
