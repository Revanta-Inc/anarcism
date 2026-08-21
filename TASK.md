You are a principal bioinformatics and Rust/WASM engineer. Design and implement a slim, browser-native replacement for the subset of ANARCI and HMMER needed to validate and number antibody variable-domain sequences.

## Objective

Build a Rust implementation compiled to WebAssembly that runs entirely in modern browsers and provides ANARCI-compatible antibody-domain recognition, chain classification, and IMGT numbering.

The complete distributable must target less than 500 kB compressed.

“Complete distributable” includes:

- The `.wasm` binary
- JavaScript/TypeScript glue
- Embedded profile-HMM/model data required at runtime
- Any lookup tables or other runtime assets

Do not exclude model assets from size accounting or download them after initialization to claim compliance.

This is not a request to port all of HMMER. Implement only the minimal profile-HMM algorithms and ANARCI behavior required for antibody VH/VL validation and IMGT numbering.

## Required biological scope

1. Detect and IMGT-number antibody variable domains:
   - Heavy (H)
   - Kappa light (K)
   - Lambda light (L)

2. Detect and IMGT-number T-cell receptor variable domains:
   - Alpha (A)
   - Beta (B)
   - Gamma (G)
   - Delta (D)

3. Bundle all applicable H/K/L/A/B/G/D profiles available in the ANARCI version
   pinned by this repository, including its human, mouse, rat, rabbit, rhesus,
   pig, alpaca, and cow coverage. Not every species necessarily has every chain
   profile.

4. Detect multiple non-overlapping variable domains within one input sequence,
   including scFv and other engineered multidomain sequences. Return domains in
   sequence order.

5. For every detected domain, return:
   - Chain type and receptor class (IG or TR)
   - Best matching species/profile
   - Zero-based input start and end coordinates
   - IMGT residue numbering with insertion codes
   - Bit score
   - Configurable number of alternative profile hits
   - HMMER-compatible E-value when it can be implemented accurately
   - Padded IMGT alignment
   - CDR1, CDR2, CDR3 and FR1–FR4 annotations

6. Support configurable:
   - Allowed chain types
   - Allowed species
   - Minimum bit-score threshold
   - Alternative-hit count
   - Germline assignment

7. Support optional closest V/J germline assignment, returning:
   - Assigned species
   - V gene and sequence identity
   - J gene and sequence identity

8. Support:
   - Single-sequence numbering
   - Batch numbering
   - FASTA input
   - VH/VL pair validation
   - Deterministic output

Public browser API

type ChainType = 'H' | 'K' | 'L' | 'A' | 'B' | 'G' | 'D';
type ReceptorType = 'IG' | 'TR';

interface NumberingOptions {
  allowedChains?: ChainType[];
  allowedSpecies?: string[];
  minBitScore?: number;
  alternativeHitCount?: number;
  assignGermline?: boolean;
}

interface NumberedResidue {
  sequenceIndex: number;
  aminoAcid: string;
  position: number;
  insertionCode: string;
  region: 'FR1' | 'CDR1' | 'FR2' | 'CDR2' | 'FR3' | 'CDR3' | 'FR4';
}

interface ProfileHit {
  profile: string;
  chainType: ChainType;
  species: string;
  bitScore: number;
  eValue?: number;
}

interface GermlineAssignment {
  species: string;
  vGene?: string;
  vIdentity?: number;
  jGene?: string;
  jIdentity?: number;
}

interface DomainResult {
  domainIndex: number;
  receptorType: ReceptorType;
  chainType: ChainType;
  species: string;
  start: number;
  end: number;
  bitScore: number;
  eValue?: number;
  numbering: NumberedResidue[];
  paddedImgtAlignment: string;
  alternativeHits: ProfileHit[];
  germline?: GermlineAssignment;
}

interface SequenceResult {
  id: string;
  normalizedSequence: string;
  domains: DomainResult[];
  warnings: string[];
}

numberSequence(
  sequence: string,
  options?: NumberingOptions
): SequenceResult;

numberSequences(
  inputs: Array<{ id: string; sequence: string }>,
  options?: NumberingOptions
): SequenceResult[];

numberFasta(
  fasta: string,
  options?: NumberingOptions
): SequenceResult[];

validateAntibodyPair(
  vh: string,
  vl: string,
  options?: NumberingOptions & {
    startMax?: number;
    endMin?: number;
  }
): {
  ok: boolean;
  vh?: DomainResult;
  vl?: DomainResult;
  errors: string[];
};

Provide both synchronous and asynchronous browser-friendly APIs if initialization requires loading or compiling WASM.

## Algorithmic scope

Determine the smallest HMMER-compatible subset necessary for ANARCI-style numbering. Investigate and document whether this requires:

- MSV filtering
- Viterbi filtering/alignment
- Forward scoring
- Profile-HMM local alignment
- Traceback
- E-value calibration
- ANARCI-specific state-to-IMGT mappings

Avoid implementing algorithms that do not materially affect the required outputs.

Do not simply compile the complete HMMER C codebase to WASM. Implement a size-focused Rust subset with only the necessary algorithms, data structures, and profiles.

Use compact binary representations for profile-HMM data. Explore:

- Quantized emission and transition scores
- Shared alphabets and transition tables
- Delta encoding
- Removing unused profiles
- Build-time generation of compact embedded assets
- Brotli/gzip-friendly data layouts
- Lazy decoding without network access
- Feature-gated diagnostic information

## Compatibility requirements

Treat a pinned ANARCI version as the reference implementation.

Create a representative golden corpus containing:

- Known human and non-human VH domains
- Kappa and lambda VL domains
- Truncated domains
- Invalid/non-antibody proteins
- Swapped VH/VL pairs
- Sequences with insertions and deletions
- Multiple domains in one sequence where supported
- Edge cases near the configured IMGT start/end boundaries

For accepted domains, measure parity for:

- Domain detection
- H/K/L classification
- Input-domain boundaries
- IMGT residue positions
- Insertion codes

Clearly define any intentional differences from ANARCI. Do not claim compatibility based only on chain classification.

## Size and performance budgets

Hard target:

- Complete compressed browser distribution: `<500 kB` Smaller is better
- Report both gzip and Brotli sizes
- Report raw WASM, optimized WASM, glue, and embedded-profile sizes separately

Performance targets on a representative laptop browser:

- Initialization: under 100 ms when cached
- Single VH/VL validation: under 50 ms
- Batch of 100 VH/VL pairs: under 1 second

Treat these as targets, measure them, and report actual results. Do not hide initialization or decompression time.

Use size-oriented Rust/WASM configuration, including where appropriate:

- `opt-level = "z"`
- LTO
- A single codegen unit
- `panic = "abort"`
- Stripped symbols
- `wasm-opt -Oz`
- Minimal allocator and formatting machinery
- Carefully selected dependencies with documented size costs

The release implementation may require WASM SIMD128, but must not require WASM threads, shared memory, filesystem access, server calls, Python, or browser extensions.

## Security and privacy

- Never transmit sequences.
- Do not log sequence prefixes, suffixes, or complete molecular payloads.
- Apply explicit sequence-length and batch-size limits.
- Avoid unbounded allocation based on user input.
- Return structured errors rather than panicking across the WASM boundary.
- Include fuzz tests for parsers, profile decoding, and sequence inputs.

## Licensing and provenance

Before implementation, audit and document:

- ANARCI license
- HMMER license
- Licenses and redistribution terms for profile-HMM/model assets
- Provenance of every embedded profile
- Whether modified or quantized profiles may legally be redistributed

Do not copy incompatible code or data. Preserve required notices.

## Required process

Start with a feasibility gate before writing the complete implementation:

1. Identify the exact ANARCI and HMMER behaviors required.
2. Inventory the required profile data and its compressed size.
3. Produce a byte-level size budget.
4. Determine whether `<500 kB` is technically and legally achievable.
5. Build a minimal proof of concept that recognizes and numbers at least one VH and one VL sequence.
6. Measure its compressed size and parity.
7. Continue to the full implementation if the evidence supports the target.

Do not stop at a high-level design. Implement and test the solution.

If the complete requirement cannot fit within 500 kB, provide evidence rather than silently weakening compatibility. Report:

- The measured minimum size
- Which assets or algorithms dominate it
- The smallest viable reduced scope
- A tiered proposal, such as:
  - `<500 kB` common human H/K/L profiles
  - Larger full-species compatibility build
  - Optional diagnostic/E-value build

## Deliverables

Produce:

1. A Rust workspace with:
   - Core sequence and profile-HMM library
   - Compact profile-data generator
   - WASM bindings
   - Browser package
2. TypeScript declarations and ergonomic browser API.
3. Golden compatibility tests against pinned ANARCI output.
4. Rust unit and property tests.
5. Browser tests in Chromium, Firefox, and WebKit.
6. Fuzz tests for untrusted inputs and profile decoding.
7. Performance benchmarks.
8. Automated size reporting and a CI failure threshold at 500 kB.
9. A minimal browser demo.
10. Documentation covering:
    - Algorithm
    - Supported scope
    - Compatibility limitations
    - Model-data provenance
    - Licensing
    - Build and release procedure
    - Size breakdown
    - Benchmark results

## Acceptance criteria

The work is complete when:

- Validation and numbering run entirely in the browser.
- No Python, native HMMER binary, backend request, or runtime model download is needed.
- Golden parity is measured and documented.
- Both single and batch APIs work.
- The browser build is deterministic and reproducible.
- The complete Brotli- or gzip-compressed distributable is below 500 kB, or a documented feasibility report proves why that target cannot be met.
- CI enforces correctness, browser compatibility, and the size budget.
```
