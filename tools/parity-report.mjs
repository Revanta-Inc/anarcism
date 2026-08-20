import { readFile } from "node:fs/promises";

import init, { numberSequence, validateAntibodyPair } from "../browser/dist/index.js";

const corpus = JSON.parse(await readFile(new URL("../tests/golden/corpus.json", import.meta.url)));
const wasm = await readFile(new URL("../browser/dist/anarcism.wasm", import.meta.url));
await init(wasm);

let domains = 0;
let residues = 0;
let exactDomains = 0;
let exactResidues = 0;
let maximumScoreDifference = 0;
const failures = [];

for (const testCase of corpus.cases) {
  const observed = numberSequence(testCase.sequence);
  if (observed.domains.length !== testCase.referenceDomains.length) {
    failures.push(`${testCase.id}: domain count`);
    continue;
  }
  for (let index = 0; index < observed.domains.length; index += 1) {
    domains += 1;
    const actual = observed.domains[index];
    const expected = testCase.referenceDomains[index];
    const domainExact = actual.chainType === expected.chainType
      && actual.species === expected.species
      && actual.start === expected.start
      && actual.end === expected.end
      && actual.paddedImgtAlignment === expected.paddedImgtAlignment;
    if (domainExact) exactDomains += 1;
    else failures.push(`${testCase.id}: domain ${index}`);
    maximumScoreDifference = Math.max(
      maximumScoreDifference,
      Math.abs(actual.bitScore - expected.bitScore),
    );
    residues += expected.numbering.length;
    for (let residueIndex = 0; residueIndex < expected.numbering.length; residueIndex += 1) {
      const left = actual.numbering[residueIndex];
      const right = expected.numbering[residueIndex];
      if (left
        && left.sequenceIndex === right.sequenceIndex
        && left.aminoAcid === right.aminoAcid
        && left.position === right.position
        && left.insertionCode === right.insertionCode) {
        exactResidues += 1;
      } else {
        failures.push(`${testCase.id}: residue ${residueIndex}`);
      }
    }
  }
}

for (const pair of corpus.pairs) {
  const observed = validateAntibodyPair(pair.vh, pair.vl);
  if (observed.ok !== pair.ok) failures.push(`${pair.id}: pair outcome`);
}

console.log(JSON.stringify({
  reference: corpus.reference,
  cases: corpus.cases.length,
  domains: { exact: exactDomains, total: domains },
  residues: { exact: exactResidues, total: residues },
  maximumAbsoluteBitScoreDifference: maximumScoreDifference,
  pairCases: corpus.pairs.length,
  failures,
}, null, 2));

if (failures.length > 0 || maximumScoreDifference > 4) process.exitCode = 1;
