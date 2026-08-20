import init, {
  type AnarcismApi,
  type ChainType,
  type NumberingOptions,
  numberFasta,
  numberSequence,
  numberSequences,
  validateAntibodyPair,
} from "../dist/index.js";

const options = {
  allowedChains: ["H", "K", "L"],
  allowedSpecies: ["human", "mouse"],
  minBitScore: 80,
  alternativeHitCount: 3,
  assignGermline: true,
} satisfies NumberingOptions;

const apiPromise: Promise<AnarcismApi> = init(new Uint8Array());
void apiPromise;
const single = numberSequence("ACDEFGHIK", options);
const chain: ChainType | undefined = single.domains[0]?.chainType;
void chain;
numberSequences([{ id: "one", sequence: "ACDEFGHIK" }], options);
numberFasta(">one\nACDEFGHIK\n", options);
validateAntibodyPair("ACDEFGHIK", "LMNPQRSTVWY", {
  ...options,
  startMax: 10,
  endMin: 100,
});
