import init, {
	type AnarcismApi,
	type ChainType,
	type NumberingOptions,
	numberFasta,
	numberSequence,
	numberSequences,
	validateAntibodyPair,
} from "../dist/index.js";
import {
	AnarcismWorkerPool,
	recommendedWorkerCount,
	type WorkerPoolMetadata,
} from "../dist/worker-pool.js";

const options = {
	allowedChains: ["H", "K", "L"],
	allowedSpecies: ["human", "mouse"],
	minBitScore: 80,
	alternativeHitCount: 3,
	assignGermline: true,
} satisfies NumberingOptions;

const apiPromise: Promise<AnarcismApi> = init(new Uint8Array());
void apiPromise;
apiPromise.then((api) => {
	const version: string = api.version;
	const chains: ChainType[] = api.chains();
	const species: string[] = api.species();
	void [version, chains, species];
});
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

const pool = new AnarcismWorkerPool();
const poolMetadata: Promise<WorkerPoolMetadata> = pool.initialize(recommendedWorkerCount());
const pooledBatch = pool.numberSequences([{ id: "one", sequence: "ACDEFGHIK" }], options);
const pooledFasta = pool.numberFasta(">one\nACDEFGHIK\n", options);
void [poolMetadata, pooledBatch, pooledFasta];
pool.terminate();
