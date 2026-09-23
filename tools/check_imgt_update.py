#!/usr/bin/env python3
"""Fail when the live IMGT/GENE-DB release differs from the pinned manifest."""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
import tomllib
import urllib.error
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import TypeVar


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "assets" / "MANIFEST.toml"
RELEASE_URL = "https://www.imgt.org/download/GENE-DB/RELEASE"
PROGRAM_VERSIONS_URL = "https://www.imgt.org/IMGTgenedbdoc/programversions.html"
REFERENCE_URL = (
    "https://www.imgt.org/download/GENE-DB/"
    "IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP"
)
USER_AGENT = "anarcism-imgt-update-check/1 (+https://github.com/revanta/anarcism)"
CHUNK_SIZE = 128 * 1024

T = TypeVar("T")


@dataclass(frozen=True)
class DatabaseState:
    release: str
    program_version: str
    reference_sha256: str


def retry(operation: Callable[[], T], *, attempts: int, label: str) -> T:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            return operation()
        except (OSError, TimeoutError, urllib.error.URLError) as error:
            last_error = error
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RuntimeError(
        f"failed to fetch {label} after {attempts} attempts"
    ) from last_error


def request(url: str) -> urllib.request.Request:
    return urllib.request.Request(
        url,
        headers={
            "Accept-Encoding": "identity",
            "User-Agent": USER_AGENT,
        },
    )


def read_url(url: str, *, timeout: float, attempts: int, label: str) -> bytes:
    def read() -> bytes:
        with urllib.request.urlopen(request(url), timeout=timeout) as response:
            return response.read()

    return retry(read, attempts=attempts, label=label)


def hash_url(url: str, *, timeout: float, attempts: int) -> str:
    def digest() -> str:
        checksum = sha256()
        with urllib.request.urlopen(request(url), timeout=timeout) as response:
            while chunk := response.read(CHUNK_SIZE):
                checksum.update(chunk)
        return checksum.hexdigest()

    return retry(digest, attempts=attempts, label="reference export")


def parse_release(content: bytes) -> str:
    release = content.decode("ascii").strip()
    if not re.fullmatch(r"\d{6}-[1-7]", release):
        raise ValueError(f"unexpected IMGT/GENE-DB release marker {release!r}")
    return release


def parse_latest_program_version(content: bytes) -> str:
    # Changelog entries are the authoritative version list.
    versions = set(re.findall(r"\bVersion\s+(\d+(?:\.\d+)+)", content.decode("utf-8")))
    if not versions:
        raise ValueError("no IMGT/GENE-DB program version found")
    return max(versions, key=lambda value: tuple(map(int, value.split("."))))


def load_expected(path: Path) -> DatabaseState:
    with path.open("rb") as manifest_file:
        reference = tomllib.load(manifest_file)["reference"]
    try:
        return DatabaseState(
            release=reference["imgt_genedb_release"],
            program_version=reference["imgt_genedb_program_version"],
            reference_sha256=reference["imgt_reference_export_sha256"],
        )
    except KeyError as error:
        raise ValueError(f"missing {error.args[0]!r} in {path}") from error


def observe(args: argparse.Namespace) -> DatabaseState:
    release = parse_release(
        read_url(
            args.release_url,
            timeout=args.timeout,
            attempts=args.retries,
            label="release marker",
        )
    )
    program_version = parse_latest_program_version(
        read_url(
            args.program_versions_url,
            timeout=args.timeout,
            attempts=args.retries,
            label="program-version page",
        )
    )
    reference_sha256 = hash_url(
        args.reference_url,
        timeout=args.timeout,
        attempts=args.retries,
    )
    return DatabaseState(release, program_version, reference_sha256)


def differences(expected: DatabaseState, observed: DatabaseState) -> list[str]:
    changes = []
    if observed.release != expected.release:
        changes.append(f"release changed: {expected.release} -> {observed.release}")
    if observed.program_version != expected.program_version:
        changes.append(
            "program version changed: "
            f"{expected.program_version} -> {observed.program_version}"
        )
    if observed.reference_sha256 != expected.reference_sha256:
        changes.append(
            "reference export SHA-256 changed: "
            f"{expected.reference_sha256} -> {observed.reference_sha256}"
        )
    return changes


def github_error(message: str) -> None:
    if os.environ.get("GITHUB_ACTIONS") == "true":
        escaped = message.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")
        print(f"::error title=IMGT/GENE-DB update detected::{escaped}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--release-url", default=RELEASE_URL)
    parser.add_argument("--program-versions-url", default=PROGRAM_VERSIONS_URL)
    parser.add_argument("--reference-url", default=REFERENCE_URL)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args(argv)
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if args.retries <= 0:
        parser.error("--retries must be positive")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        expected = load_expected(args.manifest)
        observed = observe(args)
    except (KeyError, OSError, RuntimeError, UnicodeError, ValueError) as error:
        print(f"IMGT/GENE-DB update check could not run: {error}", file=sys.stderr)
        return 2

    changes = differences(expected, observed)
    if changes:
        print("IMGT/GENE-DB update detected:", file=sys.stderr)
        for change in changes:
            print(f"- {change}", file=sys.stderr)
        github_error("; ".join(changes))
        return 1

    print(
        "IMGT/GENE-DB is current: "
        f"release {observed.release}, program {observed.program_version}, "
        f"reference SHA-256 {observed.reference_sha256}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
