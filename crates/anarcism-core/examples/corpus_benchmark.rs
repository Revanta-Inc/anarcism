use std::collections::HashMap;
use std::env;
use std::hint::black_box;
use std::num::NonZeroUsize;
use std::time::Instant;

use anarcism_core::{
    NumberingOptions, SequenceInput, number_sequence, number_sequences, number_sequences_parallel,
};
use serde::Deserialize;
use serde_json::{Value, json};

const WORKLOAD_IDS: [&str; 12] = [
    "trastuzumab_vh",
    "trastuzumab_vl",
    "human_trav12_2_traj33",
    "human_trbv19_trbj2_7",
    "human_trgv9_trgjp",
    "human_trdv2_trdj1",
    "longcdr3_cow_ultralong",
    "scfv_trastuzumab_vh_vl",
    "md_bite_pembro_okt3",
    "pdb_1hzh_h",
    "beta2microglobulin_human",
    "myoglobin_tandem_x8",
];

#[derive(Clone, Deserialize)]
struct CorpusCase {
    id: String,
    seq: String,
    category: String,
}

#[derive(Deserialize)]
struct ReferenceCase {
    id: String,
    domains: Vec<Value>,
}

fn main() {
    let warmups = setting("ANARCISM_BENCH_WARMUP", 3);
    let iterations = setting("ANARCISM_BENCH_WORKLOAD_ITERATIONS", 21);
    let corpus_iterations = setting("ANARCISM_BENCH_CORPUS_ITERATIONS", 7);
    let dense_iterations = setting("ANARCISM_BENCH_DENSE_ITERATIONS", 11);
    let dense_pairs = setting("ANARCISM_BENCH_DENSE_PAIRS", 500);
    let corpus: Vec<CorpusCase> =
        serde_json::from_str(include_str!("../../../tests/golden/corpus.json"))
            .expect("corpus is valid");
    let reference_cases: HashMap<String, ReferenceCase> =
        include_str!("../../../tests/golden/corpus_reference.jsonl")
            .lines()
            .skip(1)
            .map(|line| {
                let case: ReferenceCase =
                    serde_json::from_str(line).expect("reference record is valid");
                (case.id.clone(), case)
            })
            .collect();
    assert_eq!(corpus.len(), reference_cases.len());
    for case in &corpus {
        assert!(
            reference_cases.contains_key(&case.id),
            "{} has no reference",
            case.id
        );
    }

    let options = NumberingOptions::default();
    let workloads: Vec<_> = WORKLOAD_IDS
        .iter()
        .map(|id| {
            let index = corpus
                .iter()
                .position(|case| case.id == *id)
                .unwrap_or_else(|| panic!("missing workload {id}"));
            (&corpus[index], reference_cases[*id].domains.len())
        })
        .collect();

    for _ in 0..warmups {
        for (case, expected_domains) in &workloads {
            validate(case, *expected_domains, &options);
        }
    }
    let workload_results: Vec<_> = workloads
        .iter()
        .map(|(case, expected_domains)| {
            let samples = measure(iterations, || validate(case, *expected_domains, &options));
            json!({
                "id": case.id,
                "category": case.category,
                "residues": case.seq.len(),
                "domains": expected_domains,
                "medianMs": percentile(&samples, 0.5),
                "p95Ms": percentile(&samples, 0.95),
                "samplesMs": samples,
            })
        })
        .collect();

    let inputs: Vec<_> = corpus
        .iter()
        .map(|case| SequenceInput {
            id: case.id.clone(),
            sequence: case.seq.clone(),
        })
        .collect();
    let expected_domains: usize = reference_cases
        .values()
        .map(|case| case.domains.len())
        .sum();
    let heavy = workloads[0].0;
    let light = workloads[1].0;
    let dense_inputs: Vec<_> = (0..dense_pairs)
        .flat_map(|pair_index| {
            [
                SequenceInput {
                    id: format!("dense-vh-{pair_index}"),
                    sequence: heavy.seq.clone(),
                },
                SequenceInput {
                    id: format!("dense-vl-{pair_index}"),
                    sequence: light.seq.clone(),
                },
            ]
        })
        .collect();
    let logical_cpus = std::thread::available_parallelism().map_or(1, NonZeroUsize::get);
    let worker_counts: Vec<_> = [1, 2, 4, 8]
        .into_iter()
        .filter(|workers| *workers <= logical_cpus)
        .collect();
    let corpus_scaling: Vec<_> = worker_counts
        .iter()
        .map(|&workers| {
            let samples = measure(corpus_iterations, || {
                validate_batch(&inputs, &options, workers, expected_domains, "corpus")
            });
            throughput(workers, samples, inputs.len())
        })
        .collect();
    let dense_scaling: Vec<_> = worker_counts
        .iter()
        .map(|&workers| {
            let samples = measure(dense_iterations, || {
                validate_batch(
                    &dense_inputs,
                    &options,
                    workers,
                    dense_inputs.len(),
                    "dense batch",
                )
            });
            throughput(workers, samples, dense_inputs.len())
        })
        .collect();

    println!(
        "{}",
        serde_json::to_string(&json!({
            "implementation": "Native Rust (-O3)",
            "runtime": { "rust": env!("CARGO_PKG_RUST_VERSION") },
            "config": {
                "warmupIterations": warmups,
                "workloadIterations": iterations,
                "corpusIterations": corpus_iterations,
                "denseIterations": dense_iterations,
                "densePairs": dense_pairs,
                "denseSequences": dense_inputs.len(),
                "corpusCases": inputs.len(),
                "corpusDomains": expected_domains,
            },
            "workloads": workload_results,
            "corpus": corpus_scaling.first(),
            "parallelCorpus": corpus_scaling.last(),
            "corpusScaling": corpus_scaling,
            "dense": dense_scaling.first(),
            "parallelDense": dense_scaling.last(),
            "denseScaling": dense_scaling,
        }))
        .expect("benchmark result serializes")
    );
}

fn validate(case: &CorpusCase, expected_domains: usize, options: &NumberingOptions) {
    let result = number_sequence(&case.seq, options).expect("numbering succeeds");
    assert_eq!(result.domains.len(), expected_domains, "{}", case.id);
    black_box(result);
}

fn validate_batch(
    inputs: &[SequenceInput],
    options: &NumberingOptions,
    workers: usize,
    expected_domains: usize,
    label: &str,
) {
    let results = if workers == 1 {
        number_sequences(inputs, options)
    } else {
        number_sequences_parallel(
            inputs,
            options,
            NonZeroUsize::new(workers).expect("worker count is nonzero"),
        )
    }
    .unwrap_or_else(|error| panic!("{label} numbering succeeds: {error}"));
    let observed_domains = results
        .iter()
        .map(|result| result.domains.len())
        .sum::<usize>();
    assert_eq!(observed_domains, expected_domains, "{label}");
    black_box(results);
}

fn measure(mut iterations: usize, mut operation: impl FnMut()) -> Vec<f64> {
    let mut samples = Vec::with_capacity(iterations);
    while iterations > 0 {
        let started = Instant::now();
        operation();
        samples.push(started.elapsed().as_secs_f64() * 1_000.0);
        iterations -= 1;
    }
    samples.sort_by(f64::total_cmp);
    samples
}

fn percentile(samples: &[f64], fraction: f64) -> f64 {
    samples[((samples.len() as f64 * fraction).floor() as usize).min(samples.len() - 1)]
}

fn throughput(workers: usize, samples: Vec<f64>, sequences: usize) -> Value {
    let total_ms = percentile(&samples, 0.5);
    json!({
        "workers": workers,
        "totalMs": total_ms,
        "p95TotalMs": percentile(&samples, 0.95),
        "samplesMs": samples,
        "sequences": sequences,
        "perSequenceMs": total_ms / sequences as f64,
        "sequencesPerSecond": sequences as f64 * 1_000.0 / total_ms,
    })
}

fn setting(name: &str, default: usize) -> usize {
    env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .filter(|value| *value > 0)
        .unwrap_or(default)
}
