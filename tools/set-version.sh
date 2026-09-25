#!/usr/bin/env bash
set -euo pipefail

# Sets the release version everywhere it is recorded. The Python package and
# the WebAssembly module read theirs from the Cargo workspace.
if [[ $# -ne 1 || ! "$1" =~ ^[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.-]+)?$ ]]; then
  echo "usage: $0 <major.minor.patch[-prerelease]>" >&2
  exit 2
fi
export VERSION="$1"

repository_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repository_root"

perl -pi -e 's/^version = ".*"$/version = "$ENV{VERSION}"/ if /^\[workspace\.package\]$/ ... /^\[/' Cargo.toml
perl -pi -e 's/^version: .*$/version: $ENV{VERSION}/' CITATION.cff
npm version "$VERSION" --workspaces --include-workspace-root --no-git-tag-version --allow-same-version > /dev/null
cargo update --workspace --offline --quiet
cargo update --workspace --offline --quiet --manifest-path fuzz/Cargo.toml

git status --short
