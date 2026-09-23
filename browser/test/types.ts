import {
	Anarcism,
	type ChainType,
	type NumberingOptions,
	type PairValidationResult,
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

const anarcismPromise: Promise<Anarcism> = Anarcism.create({ source: new Uint8Array() });
anarcismPromise.then((anarcism) => {
	const version: string = anarcism.version;
	const chains: ChainType[] = anarcism.chains();
	const species: string[] = anarcism.species();
	void [version, chains, species];

	const single = anarcism.numberSequence("ACDEFGHIK", options);
	const chain: ChainType | undefined = single.domains[0]?.chainType;
	void chain;
	anarcism.numberSequences([{ id: "one", sequence: "ACDEFGHIK" }], options);
	anarcism.numberFasta(">one\nACDEFGHIK\n", options);
	const pair: PairValidationResult = anarcism.validateAntibodyPair("ACDEFGHIK", "LMNPQRSTVWY", {
		...options,
		startMax: 10,
		endMin: 100,
	});
	void pair;
});
const controller = new AbortController();
void Anarcism.create({ signal: controller.signal });
const syncAnarcism: Anarcism = Anarcism.createSync(new Uint8Array());
void syncAnarcism;
// @ts-expect-error instances come from Anarcism.create()
void new Anarcism();

const poolPromise: Promise<AnarcismWorkerPool> = AnarcismWorkerPool.create({
	workers: recommendedWorkerCount(),
	maxWorkers: 8,
});
// @ts-expect-error pools come from AnarcismWorkerPool.create()
void new AnarcismWorkerPool();
const pool = await poolPromise;
const poolVersion: string = pool.version;
const poolChains: ChainType[] = pool.chains();
const resized: Promise<WorkerPoolMetadata> = pool.resize(2, { signal: controller.signal });
void AnarcismWorkerPool.create({ signal: controller.signal });
void pool.numberSequences([], { assignGermline: true, signal: controller.signal });
void pool.validateAntibodyPair("A", "B", { startMax: 10, signal: controller.signal });
const pooledBatch = pool.numberSequences([{ id: "one", sequence: "ACDEFGHIK" }], options);
const pooledFasta = pool.numberFasta(">one\nACDEFGHIK\n", options);
const pooledPair: Promise<PairValidationResult> = pool.validateAntibodyPair(
	"ACDEFGHIK",
	"LMNPQRSTVWY",
	{ ...options, startMax: 10, endMin: 100 },
);
void [poolVersion, poolChains, resized, pooledBatch, pooledFasta, pooledPair];
pool.terminate();
