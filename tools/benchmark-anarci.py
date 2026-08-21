#!/usr/bin/env python3
"""Benchmark the pinned native ANARCI/HMMER reference on the shared workload."""

from __future__ import annotations

import json
import os
import platform
import subprocess
import time
from importlib.metadata import version
from typing import Any

from anarci import anarci

VH = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA"
VL = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR"
PINNED_ANARCI = "2026.2.13.2"


def setting(name: str, default: int) -> int:
    try:
        value = int(os.environ.get(name, default))
    except ValueError:
        return default
    return value if value > 0 else default


def number(sequences: list[tuple[str, str]]) -> tuple[Any, Any, Any]:
    return anarci(
        sequences,
        scheme="imgt",
        allowed_species=None,
        assign_germline=False,
        bit_score_threshold=80,
        ncpu=1,
    )


def validate(sequence: str, expected_chain: str) -> None:
    numbered, details, _ = number([("sequence", sequence)])
    validate_result(numbered[0], details[0], expected_chain)


def validate_result(numbered: Any, details: Any, expected_chain: str) -> None:
    if not numbered or len(numbered) != 1:
        raise RuntimeError("ANARCI did not return exactly one domain")
    if details[0]["chain_type"] != expected_chain:
        raise RuntimeError(f"ANARCI returned the wrong chain for {expected_chain}")


def distribution(sequence: str, expected_chain: str, iterations: int) -> dict[str, float]:
    samples = []
    for _ in range(iterations):
        started = time.perf_counter()
        numbered, details, _ = number([("sequence", sequence)])
        elapsed_ms = (time.perf_counter() - started) * 1_000
        validate_result(numbered[0], details[0], expected_chain)
        samples.append(elapsed_ms)
    samples.sort()
    return {
        "medianMs": percentile(samples, 0.5),
        "p95Ms": percentile(samples, 0.95),
    }


def percentile(samples: list[float], fraction: float) -> float:
    return samples[min(int(len(samples) * fraction), len(samples) - 1)]


def batch_metrics(total_ms: float, sequence_count: int) -> dict[str, float]:
    return {
        "totalMs": total_ms,
        "perSequenceMs": total_ms / sequence_count,
        "sequencesPerSecond": sequence_count * 1_000 / total_ms,
    }


def main() -> None:
    anarci_version = version("anarci")
    if anarci_version != PINNED_ANARCI:
        raise RuntimeError(
            f"benchmark requires ANARCI {PINNED_ANARCI}, found {anarci_version}"
        )
    hmmer = subprocess.run(
        ["hmmscan", "-h"], check=True, capture_output=True, text=True
    ).stdout.splitlines()[1].removeprefix("# ")
    if not hmmer.startswith("HMMER 3.4 "):
        raise RuntimeError(f"benchmark requires HMMER 3.4, found {hmmer}")

    warmup_iterations = setting("ANARCISM_BENCH_WARMUP", 3)
    single_iterations = setting("ANARCISM_BENCH_SINGLE_ITERATIONS", 25)
    pair_count = setting("ANARCISM_BENCH_PAIR_COUNT", 100)

    for _ in range(warmup_iterations):
        validate(VH, "H")
        validate(VL, "K")

    vh = distribution(VH, "H", single_iterations)
    vl = distribution(VL, "K", single_iterations)
    inputs = [
        item
        for pair_index in range(pair_count)
        for item in ((f"vh-{pair_index}", VH), (f"vl-{pair_index}", VL))
    ]
    started = time.perf_counter()
    numbered, details, _ = number(inputs)
    batch_ms = (time.perf_counter() - started) * 1_000
    if len(numbered) != len(inputs) or len(details) != len(inputs):
        raise RuntimeError("ANARCI returned an incomplete batch")
    for index, (sequence_numbering, sequence_details) in enumerate(
        zip(numbered, details, strict=True)
    ):
        validate_result(
            sequence_numbering,
            sequence_details,
            "H" if index % 2 == 0 else "K",
        )

    print(
        json.dumps(
            {
                "implementation": "ANARCI/HMMER",
                "runtime": {
                    "python": platform.python_version(),
                    "anarci": anarci_version,
                    "hmmer": hmmer,
                },
                "config": {
                    "warmupIterations": warmup_iterations,
                    "singleIterations": single_iterations,
                    "pairCount": pair_count,
                    "sequenceCount": len(inputs),
                    "threads": 1,
                },
                "vh": vh,
                "vl": vl,
                "batch": batch_metrics(batch_ms, len(inputs)),
            }
        )
    )


if __name__ == "__main__":
    main()
