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
struct Reference {
    cases: Vec<ReferenceCase>,
}

#[derive(Deserialize)]
struct ReferenceCase {
    id: String,
    domains: Vec<Value>,
}

fn main() {
    let warmups = setting("ANARCISM_BENCH_WARMUP", 2);
    let iterations = setting("ANARCISM_BENCH_WORKLOAD_ITERATIONS", 9);
    let corpus_iterations = setting("ANARCISM_BENCH_CORPUS_ITERATIONS", 3);
    let scaling_repetitions = setting("ANARCISM_BENCH_SCALING_REPETITIONS", 24);
    let corpus: Vec<CorpusCase> =
        serde_json::from_str(include_str!("../../../tests/golden/corpus_v2.json"))
            .expect("corpus is valid");
    let reference: Reference = serde_json::from_str(include_str!(
        "../../../tests/golden/corpus_v2_reference.json"
    ))
    .expect("reference is valid");
    assert_eq!(corpus.len(), reference.cases.len());
    for (case, expected) in corpus.iter().zip(&reference.cases) {
        assert_eq!(case.id, expected.id);
    }

    let options = NumberingOptions::default();
    let workloads: Vec<_> = WORKLOAD_IDS
        .iter()
        .map(|id| {
            let index = corpus
                .iter()
                .position(|case| case.id == *id)
                .unwrap_or_else(|| panic!("missing workload {id}"));
            (&corpus[index], reference.cases[index].domains.len())
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
    let expected_domains: usize = reference.cases.iter().map(|case| case.domains.len()).sum();
    let corpus_samples = measure(corpus_iterations, || {
        let mut observed_domains = 0;
        for chunk in inputs.chunks(1_000) {
            let results = number_sequences(chunk, &options).expect("corpus numbering succeeds");
            observed_domains += results
                .iter()
                .map(|result| result.domains.len())
                .sum::<usize>();
            black_box(results);
        }
        assert_eq!(observed_domains, expected_domains);
    });

    let scaling_inputs: Vec<_> = (0..scaling_repetitions)
        .flat_map(|repetition| {
            workloads.iter().map(move |(case, _)| SequenceInput {
                id: format!("{}-{repetition}", case.id),
                sequence: case.seq.clone(),
            })
        })
        .collect();
    let logical_cpus = std::thread::available_parallelism().map_or(1, NonZeroUsize::get);
    let scaling: Vec<_> = [1, 2, 4, 8]
        .into_iter()
        .filter(|workers| *workers <= logical_cpus)
        .map(|workers| {
            let samples = measure(3, || {
                let results = number_sequences_parallel(
                    &scaling_inputs,
                    &options,
                    NonZeroUsize::new(workers).expect("worker count is nonzero"),
                )
                .expect("parallel numbering succeeds");
                black_box(results);
            });
            let total_ms = percentile(&samples, 0.5);
            throughput(workers, total_ms, scaling_inputs.len())
        })
        .collect();
    let corpus_ms = percentile(&corpus_samples, 0.5);

    println!(
        "{}",
        serde_json::to_string(&json!({
            "implementation": "Native Rust (-O3)",
            "runtime": { "rust": env!("CARGO_PKG_RUST_VERSION") },
            "config": {
                "warmupIterations": warmups,
                "workloadIterations": iterations,
                "corpusIterations": corpus_iterations,
                "corpusCases": inputs.len(),
                "corpusDomains": expected_domains,
                "scalingSequences": scaling_inputs.len(),
            },
            "workloads": workload_results,
            "corpus": throughput(1, corpus_ms, inputs.len()),
            "scaling": scaling,
        }))
        .expect("benchmark result serializes")
    );
}

fn validate(case: &CorpusCase, expected_domains: usize, options: &NumberingOptions) {
    let result = number_sequence(&case.seq, options).expect("numbering succeeds");
    assert_eq!(result.domains.len(), expected_domains, "{}", case.id);
    black_box(result);
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

fn throughput(workers: usize, total_ms: f64, sequences: usize) -> Value {
    json!({
        "workers": workers,
        "totalMs": total_ms,
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
