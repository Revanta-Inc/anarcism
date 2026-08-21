import { spawn } from "node:child_process";
import { cpus, platform, release } from "node:os";
import { resolve } from "node:path";

const root = resolve(import.meta.dirname, "..");
const environment = { ...process.env };
const runners = [
  {
    label: "Native Rust (-O3)",
    command: "cargo",
    args: ["run", "--quiet", "--profile", "benchmark", "-p", "anarcism-core", "--example", "runtime_benchmark"],
    cwd: root,
  },
  {
    label: "Browser WASM",
    command: process.execPath,
    args: ["bench/browser-benchmark.mjs"],
    cwd: resolve(root, "browser"),
  },
  {
    label: "ANARCI/HMMER",
    command: resolve(root, ".venv/bin/python"),
    args: ["tools/benchmark-anarci.py"],
    cwd: root,
  },
];

const results = [];
for (const [index, runner] of runners.entries()) {
  process.stderr.write(`[${index + 1}/${runners.length}] ${runner.label}\n`);
  results.push(await run(runner));
}

const config = results[0].config;
for (const result of results.slice(1)) {
  if (Object.keys(config).some((key) => result.config[key] !== config[key])) {
    throw new Error(`benchmark configuration drifted for ${result.implementation}`);
  }
}

const report = {
  machine: {
    platform: `${platform()} ${release()}`,
    architecture: process.arch,
    cpu: cpus()[0]?.model ?? "unknown",
    logicalCpus: cpus().length,
  },
  config,
  results,
};

if (process.argv.includes("--json")) {
  console.log(JSON.stringify(report, null, 2));
} else {
  printMarkdown(report);
}

function run({ command, args, cwd, label }) {
  return new Promise((resolveRun, rejectRun) => {
    const child = spawn(command, args, {
      cwd,
      env: environment,
      stdio: ["ignore", "pipe", "pipe"],
    });
    let stdout = "";
    let stderr = "";
    child.stdout.setEncoding("utf8").on("data", (chunk) => stdout += chunk);
    child.stderr.setEncoding("utf8").on("data", (chunk) => stderr += chunk);
    child.on("error", rejectRun);
    child.on("close", (code) => {
      if (code !== 0) {
        rejectRun(new Error(`${label} exited with ${code}\n${stderr}`));
        return;
      }
      try {
        resolveRun(JSON.parse(stdout));
      } catch (error) {
        rejectRun(new Error(`${label} returned invalid JSON: ${error.message}\n${stdout}\n${stderr}`));
      }
    });
  });
}

function printMarkdown(report) {
  const nativeBatch = report.results[0].batch.totalMs;
  console.log(`# Cross-runtime IMGT benchmark\n`);
  console.log(`${report.machine.cpu} (${report.machine.architecture}), ${report.machine.logicalCpus} logical CPUs`);
  console.log(`${report.config.singleIterations} hot samples per single sequence; ${report.config.pairCount} VH/VL pairs in one batch; one worker per runtime.\n`);
  console.log("| Runtime | VH median / p95 | VL median / p95 | Batch total | ms/sequence | sequences/s | Batch vs Rust |");
  console.log("|---|---:|---:|---:|---:|---:|---:|");
  for (const result of report.results) {
    console.log(
      `| ${result.implementation} | ${durationPair(result.vh)} | ${durationPair(result.vl)} | ` +
      `${fixed(result.batch.totalMs)} ms | ${fixed(result.batch.perSequenceMs)} | ` +
      `${fixed(result.batch.sequencesPerSecond)} | ${(result.batch.totalMs / nativeBatch).toFixed(2)}× |`,
    );
  }
}

function durationPair(distribution) {
  return `${fixed(distribution.medianMs)} / ${fixed(distribution.p95Ms)} ms`;
}

function fixed(value) {
  return value.toFixed(value < 10 ? 3 : 1);
}
