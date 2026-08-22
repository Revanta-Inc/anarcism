use std::env;
use std::hint::black_box;
use std::time::Instant;

use anarcism_core::{
    ChainType, NumberingOptions, SequenceInput, number_sequence, number_sequence_with_id,
    number_sequences,
};
use serde_json::json;

const VH: &str = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";
const VL: &str = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR";

fn main() {
    let warmup_iterations = setting("ANARCISM_BENCH_WARMUP", 3);
    let single_iterations = setting("ANARCISM_BENCH_SINGLE_ITERATIONS", 25);
    let pair_count = setting("ANARCISM_BENCH_PAIR_COUNT", 100);
    let options = NumberingOptions::default();

    for _ in 0..warmup_iterations {
        validate(VH, ChainType::H, &options);
        validate(VL, ChainType::K, &options);
    }

    let vh = distribution(VH, ChainType::H, &options, single_iterations);
    let vl = distribution(VL, ChainType::K, &options, single_iterations);
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
    let scalar_started = Instant::now();
    let scalar_batch: Vec<_> = inputs
        .iter()
        .map(|input| {
            number_sequence_with_id(&input.id, &input.sequence, &options)
                .expect("scalar batch numbering succeeds")
        })
        .collect();
    let scalar_batch_ms = scalar_started.elapsed().as_secs_f64() * 1_000.0;
    let started = Instant::now();
    let batch = number_sequences(&inputs, &options).expect("native batch numbering succeeds");
    let batch_ms = started.elapsed().as_secs_f64() * 1_000.0;
    assert_eq!(batch, scalar_batch);
    assert_eq!(batch.len(), inputs.len());
    for (index, result) in batch.iter().enumerate() {
        let expected_chain = if index % 2 == 0 {
            ChainType::H
        } else {
            ChainType::K
        };
        assert_single_domain(result, expected_chain);
    }
    black_box(batch);

    println!(
        "{}",
        serde_json::to_string(&json!({
            "implementation": "Native Rust (-O3)",
            "config": {
                "warmupIterations": warmup_iterations,
                "singleIterations": single_iterations,
                "pairCount": pair_count,
                "sequenceCount": inputs.len(),
                "threads": 1,
            },
            "vh": vh,
            "vl": vl,
            "scalarBatch": batch_metrics(scalar_batch_ms, inputs.len()),
            "batch": batch_metrics(batch_ms, inputs.len()),
            "batchSpeedup": scalar_batch_ms / batch_ms,
        }))
        .expect("benchmark result serializes")
    );
}

fn distribution(
    sequence: &str,
    expected_chain: ChainType,
    options: &NumberingOptions,
    iterations: usize,
) -> serde_json::Value {
    let mut samples = Vec::with_capacity(iterations);
    for _ in 0..iterations {
        let started = Instant::now();
        let result = number_sequence(sequence, options).expect("native numbering succeeds");
        let elapsed_ms = started.elapsed().as_secs_f64() * 1_000.0;
        assert_single_domain(&result, expected_chain);
        black_box(result);
        samples.push(elapsed_ms);
    }
    samples.sort_by(f64::total_cmp);
    json!({
        "medianMs": percentile(&samples, 0.5),
        "p95Ms": percentile(&samples, 0.95),
    })
}

fn validate(sequence: &str, expected_chain: ChainType, options: &NumberingOptions) {
    let result = number_sequence(sequence, options).expect("native numbering succeeds");
    assert_single_domain(&result, expected_chain);
    black_box(result);
}

fn assert_single_domain(result: &anarcism_core::SequenceResult, expected_chain: ChainType) {
    assert_eq!(result.domains.len(), 1);
    assert_eq!(result.domains[0].chain_type, expected_chain);
}

fn percentile(samples: &[f64], percentile: f64) -> f64 {
    let index = ((samples.len() as f64 * percentile).floor() as usize).min(samples.len() - 1);
    samples[index]
}

fn batch_metrics(total_ms: f64, sequence_count: usize) -> serde_json::Value {
    json!({
        "totalMs": total_ms,
        "perSequenceMs": total_ms / sequence_count as f64,
        "sequencesPerSecond": sequence_count as f64 * 1_000.0 / total_ms,
    })
}

fn setting(name: &str, default: usize) -> usize {
    env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .filter(|value| *value > 0)
        .unwrap_or(default)
}
