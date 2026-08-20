import init, { numberSequence } from "../browser/dist/index.js";

const button = document.querySelector("#number");
const sequence = document.querySelector("#sequence");
const germline = document.querySelector("#germline");
const status = document.querySelector("#status");
const output = document.querySelector("#result");

try {
  await init("../browser/dist/anarcism.wasm");
  button.disabled = false;
  button.textContent = "Number sequence";
  status.textContent = "Ready";
} catch (error) {
  status.textContent = `Initialization failed: ${error.message}`;
}

button.addEventListener("click", () => {
  const started = performance.now();
  try {
    const result = numberSequence(sequence.value, {
      assignGermline: germline.checked,
    });
    status.textContent = `Completed locally in ${(performance.now() - started).toFixed(1)} ms`;
    output.textContent = JSON.stringify(result, null, 2);
  } catch (error) {
    status.textContent = `${error.code ?? "ERROR"}: ${error.message}`;
    output.textContent = "";
  }
});
