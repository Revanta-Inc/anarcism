use std::env;
use std::hint::black_box;
use std::num::NonZeroUsize;
use std::time::Instant;

use anarcism_core::{
    ChainType, NumberingOptions, SequenceInput, SequenceResult, number_sequences,
    number_sequences_parallel,
};
use serde_json::json;

const VH: &str = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";
const VL: &str = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR";

fn main() {
    let pair_count = setting("ANARCISM_BENCH_PAIR_COUNT", 100);
    let iterations = setting("ANARCISM_BENCH_BATCH_ITERATIONS", 3);
    let available_workers = std::thread::available_parallelism().map_or(1, NonZeroUsize::get);
    let options = NumberingOptions::default();
    let inputs: Vec<_> = (0..pair_count)
        .flat_map(|pair_index| {
            [
                SequenceInput {
                    id: format!("vh-{pair_index}"),
                    sequence: VH.to_owned(),
                },
                SequenceInput {
                    id: format!("vl-{pair_index}"),
                    sequence: VL.to_owned(),
                },
            ]
        })
        .collect();

    let worker_counts = worker_counts(available_workers, inputs.len());
    let mut results = Vec::with_capacity(worker_counts.len());
    for worker_count in worker_counts {
        let workers = NonZeroUsize::new(worker_count).expect("worker count is positive");
        let warmup = run_batch(
            &inputs[..inputs.len().min(worker_count * 2)],
            &options,
            workers,
        );
        validate(&warmup);

        let mut samples = Vec::with_capacity(iterations);
        for _ in 0..iterations {
            let started = Instant::now();
            let batch = run_batch(&inputs, &options, workers);
            let elapsed_ms = started.elapsed().as_secs_f64() * 1_000.0;
            validate(&batch);
            black_box(batch);
            samples.push(elapsed_ms);
        }
        samples.sort_by(f64::total_cmp);
        let median_ms = samples[samples.len() / 2];
        results.push(json!({
            "workers": worker_count,
            "medianMs": median_ms,
            "perSequenceMs": median_ms / inputs.len() as f64,
            "sequencesPerSecond": inputs.len() as f64 * 1_000.0 / median_ms,
        }));
    }

    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "implementation": "Native Rust batch scaling (-O3)",
            "availableWorkers": available_workers,
            "pairCount": pair_count,
            "sequenceCount": inputs.len(),
            "iterations": iterations,
            "results": results,
        }))
        .expect("benchmark result serializes")
    );
}

fn run_batch(
    inputs: &[SequenceInput],
    options: &NumberingOptions,
    workers: NonZeroUsize,
) -> Vec<SequenceResult> {
    if workers.get() == 1 {
        number_sequences(inputs, options).expect("serial batch numbering succeeds")
    } else {
        number_sequences_parallel(inputs, options, workers)
            .expect("parallel batch numbering succeeds")
    }
}

fn validate(results: &[SequenceResult]) {
    for (index, result) in results.iter().enumerate() {
        assert_eq!(result.domains.len(), 1);
        assert_eq!(
            result.domains[0].chain_type,
            if index % 2 == 0 {
                ChainType::H
            } else {
                ChainType::K
            }
        );
    }
}

fn worker_counts(available: usize, input_count: usize) -> Vec<usize> {
    let maximum = available.min(input_count.max(1));
    let mut counts = vec![1];
    for count in [2, 4, 8, maximum] {
        if count <= maximum && counts.last() != Some(&count) {
            counts.push(count);
        }
    }
    counts
}

fn setting(name: &str, default: usize) -> usize {
    env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .filter(|value| *value > 0)
        .unwrap_or(default)
}
