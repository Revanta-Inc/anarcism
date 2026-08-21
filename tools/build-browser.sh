#!/usr/bin/env bash
set -euo pipefail

repository_root="$(cd "$(dirname "$0")/.." && pwd)"
target_wasm="$repository_root/target/wasm32-unknown-unknown/release/anarcism_wasm.wasm"
distribution="$repository_root/browser/dist"

cd "$repository_root"
cargo build --locked --release --target wasm32-unknown-unknown -p anarcism-wasm
mkdir -p "$distribution"
if command -v wasm-opt >/dev/null 2>&1; then
  wasm-opt -Oz --strip-debug \
    --enable-simd \
    --enable-sign-ext \
    --enable-nontrapping-float-to-int \
    --enable-bulk-memory \
    --enable-bulk-memory-opt \
    "$target_wasm" -o "$distribution/anarcism.wasm"
else
  cp "$target_wasm" "$distribution/anarcism.wasm"
  echo "warning: wasm-opt was not found; copied the Rust-optimized artifact" >&2
fi
cp browser/src/index.js browser/src/index.d.ts "$distribution/"
cp THIRD_PARTY_NOTICES.md "$distribution/"
node tools/size-report.mjs
