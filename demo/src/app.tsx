import { useEffect, useLayoutEffect, useRef, useState } from "preact/hooks";
import type { JSX } from "preact";
import { DEFAULT_SEQUENCE, EXAMPLES } from "./examples.ts";
import { Results } from "./results.tsx";
import {
	MAX_WORKER_COUNT,
	createAnalysisPool,
	errorCode,
	errorMessage,
	normalizeWorkerCount,
	parseFasta,
	recommendedWorkerCount,
} from "./anarcism.ts";
import type { AnalysisPool, ChainType, NumberingOptions, SequenceResult } from "./anarcism.ts";

const FIELD_LABEL = "text-xs font-semibold uppercase tracking-wide text-zinc-500";
const CONTROL =
	"mt-1.5 w-full rounded-md border border-zinc-300 bg-white px-2 py-1.5 font-mono text-xs text-zinc-900 focus:border-zinc-500 focus:outline-none";
const PRIMARY_BUTTON =
	"rounded-md bg-blue-600 px-3.5 py-2 text-sm font-semibold text-white disabled:opacity-50";
const GHOST_BUTTON =
	"rounded-md border border-zinc-300 px-2 py-1 text-xs text-blue-700 disabled:opacity-50";
const PANEL = "rounded-lg border border-zinc-200 bg-white p-4";

interface Settings {
	chains: ChainType[];
	species: string[];
	minBitScore: string;
	workerCount: number;
	assignGermline: boolean;
	alternativeHits: boolean;
}

interface Timings {
	initMs?: number;
	runMs?: number;
	perMs?: number;
}

export function App() {
	const pool = useRef<AnalysisPool | null>(null);
	const autoRan = useRef(false);
	const [version, setVersion] = useState<string | null>(null);
	const [chains, setChains] = useState<ChainType[]>([]);
	const [species, setSpecies] = useState<string[]>([]);
	const [poolSize, setPoolSize] = useState(0);
	const [ready, setReady] = useState(false);
	const [busy, setBusy] = useState<string | null>("Loading module…");
	const [error, setError] = useState<string | null>(null);
	const [results, setResults] = useState<SequenceResult[] | null>(null);
	const [timings, setTimings] = useState<Timings>({});
	const [sequence, setSequence] = useState(DEFAULT_SEQUENCE);
	const [settings, setSettings] = useState<Settings>({
		chains: [],
		species: [],
		minBitScore: "80",
		workerCount: recommendedWorkerCount(),
		assignGermline: false,
		alternativeHits: false,
	});

	// One worker keeps single-sequence latency low; FASTA batches grow the pool
	// lazily in `run`.
	useEffect(() => {
		const started = performance.now();
		let cancelled = false;

		createAnalysisPool(1).then(
			(analysisPool) => {
				if (cancelled) {
					analysisPool.terminate();
					return;
				}
				pool.current = analysisPool;
				setTimings({ initMs: performance.now() - started });
				setVersion(analysisPool.version);
				setChains(analysisPool.chains());
				setSpecies(analysisPool.species());
				setPoolSize(analysisPool.size);
				setReady(true);
				setBusy(null);
			},
			(failure: unknown) => {
				if (cancelled) return;
				setError(`Initialization failed: ${errorMessage(failure)}`);
				setBusy(null);
			},
		);

		return () => {
			cancelled = true;
			pool.current?.terminate();
			pool.current = null;
		};
	}, []);

	// Number the initial input as soon as the pool is ready. A layout effect
	// starts the run before the first ready frame paints.
	useLayoutEffect(() => {
		if (!ready || autoRan.current) return;
		autoRan.current = true;
		void run(sequence);
	}, [ready]);

	function patch(changes: Partial<Settings>) {
		setSettings((current) => ({ ...current, ...changes }));
	}

	async function growPool(target: number) {
		const analysisPool = pool.current;
		if (!analysisPool || target <= analysisPool.size) return;
		const started = performance.now();
		await analysisPool.resize(target);
		setTimings((current) => ({ ...current, initMs: performance.now() - started }));
		setPoolSize(analysisPool.size);
	}

	async function run(text: string) {
		const analysisPool = pool.current;
		const input = text.trim();
		if (!analysisPool || !ready || busy || !input) return;

		setBusy("Numbering…");
		const started = performance.now();
		try {
			const options = requestOptions(settings);
			let next: SequenceResult[];
			if (input.startsWith(">")) {
				const inputs = parseFasta(input);
				await growPool(Math.min(settings.workerCount, inputs.length));
				next = await analysisPool.numberSequences(inputs, options);
			} else {
				next = [await analysisPool.numberSequence(input, options)];
			}
			const runMs = performance.now() - started;
			setTimings((current) => ({
				...current,
				runMs,
				perMs: runMs / Math.max(next.length, 1),
			}));
			setError(null);
			setResults(next);
		} catch (failure) {
			// `code` is the stable identifier and the message never carries sequence
			// content, so both are safe to show.
			setResults(null);
			setTimings((current) => ({ ...current, runMs: undefined, perMs: undefined }));
			setError(`${errorCode(failure)}: ${errorMessage(failure)}`);
		} finally {
			setBusy(null);
		}
	}

	// Only shrinking needs an immediate resize; growth happens when a batch
	// actually arrives.
	async function changeWorkerCount(value: string) {
		const analysisPool = pool.current;
		const requested = normalizeWorkerCount(Number(value), MAX_WORKER_COUNT);
		patch({ workerCount: requested });
		if (!analysisPool || requested >= analysisPool.size) return;

		setReady(false);
		setBusy("Resizing pool…");
		const started = performance.now();
		try {
			await analysisPool.resize(requested);
			setTimings((current) => ({ ...current, initMs: performance.now() - started }));
			setError(null);
		} catch (failure) {
			patch({ workerCount: Math.max(analysisPool.size, 1) });
			setError(`Worker pool resize failed: ${errorMessage(failure)}`);
		} finally {
			setPoolSize(analysisPool.size);
			setReady(analysisPool.size > 0);
			setBusy(null);
		}
	}

	const disabled = !ready || busy !== null;

	return (
		<div class="mx-auto max-w-6xl px-5 py-6">
			<header class="mb-5 flex items-center gap-4">
				{/* The mark carries the wordmark, so it stands in for the heading text. */}
				<h1 class="shrink-0">
					<img
						src={`${import.meta.env.BASE_URL}anarcism.png`}
						alt="anarcism"
						width={500}
						height={499}
						class="h-28 w-auto"
					/>
				</h1>
				<p class="max-w-2xl text-sm text-zinc-600">
					Antibody and TCR variable-domain recognition and IMGT numbering, running in a configurable
					pool of dedicated workers. Sequences never leave this tab.
				</p>
			</header>

			<div class="grid items-start gap-5 lg:grid-cols-[21rem_1fr]">
				<section class={`${PANEL} space-y-4`}>
					<Field id="sequence" label="Sequence or FASTA">
						<textarea
							id="sequence"
							class={`${CONTROL} h-36 resize-y break-all`}
							spellcheck={false}
							placeholder="Paste one sequence, or FASTA with > headers"
							value={sequence}
							onInput={(event) => setSequence(event.currentTarget.value)}
							onKeyDown={(event) => {
								if ((event.metaKey || event.ctrlKey) && event.key === "Enter") {
									void run(event.currentTarget.value);
								}
							}}
						/>
					</Field>

					<div class="flex flex-wrap gap-1.5">
						{EXAMPLES.map((example) => (
							<button
								key={example.label}
								type="button"
								class={GHOST_BUTTON}
								disabled={disabled}
								onClick={() => {
									setSequence(example.sequence);
									void run(example.sequence);
								}}
							>
								{example.label}
							</button>
						))}
					</div>

					<Field label="Allowed chains">
						<div class="mt-1.5 flex flex-wrap gap-x-3 gap-y-1.5">
							{chains.map((chain) => (
								<Checkbox
									key={chain}
									label={chain}
									checked={settings.chains.includes(chain)}
									onChange={(checked) => {
										patch({
											chains: checked
												? [...settings.chains, chain]
												: settings.chains.filter((value) => value !== chain),
										});
									}}
								/>
							))}
						</div>
					</Field>

					<Field id="species" label="Allowed species" hint="Nothing selected means every species.">
						<select
							id="species"
							multiple
							class={`${CONTROL} h-28 font-sans`}
							onChange={(event) => {
								patch({
									species: [...event.currentTarget.selectedOptions].map((option) => option.value),
								});
							}}
						>
							{species.map((name) => (
								<option key={name} value={name}>
									{name}
								</option>
							))}
						</select>
					</Field>

					<Field id="min-bit-score" label="Minimum bit score">
						<input
							id="min-bit-score"
							type="number"
							min={0}
							step={1}
							class={CONTROL}
							value={settings.minBitScore}
							onInput={(event) => patch({ minBitScore: event.currentTarget.value })}
						/>
					</Field>

					<Field
						id="workers"
						label="Batch workers"
						hint="Extra workers start lazily for FASTA batches; a single sequence uses one worker."
					>
						<input
							id="workers"
							type="number"
							min={1}
							max={MAX_WORKER_COUNT}
							step={1}
							class={CONTROL}
							disabled={disabled}
							value={settings.workerCount}
							onChange={(event) => void changeWorkerCount(event.currentTarget.value)}
						/>
					</Field>

					<Field label="Options">
						<div class="mt-1.5 flex flex-wrap gap-x-3 gap-y-1.5">
							<Checkbox
								label="Assign germline"
								checked={settings.assignGermline}
								onChange={(checked) => patch({ assignGermline: checked })}
							/>
							<Checkbox
								label="Alternative hits"
								checked={settings.alternativeHits}
								onChange={(checked) => patch({ alternativeHits: checked })}
							/>
						</div>
					</Field>

					<button
						type="button"
						class={PRIMARY_BUTTON}
						disabled={disabled}
						onClick={() => void run(sequence)}
					>
						{busy ?? "Number"}
					</button>

					<dl class="flex flex-wrap gap-x-4 gap-y-1 border-t border-zinc-200 pt-3 text-xs text-zinc-500">
						<Stat label="pool ready" value={formatMs(timings.initMs)} />
						<Stat label="workers ready/max" value={`${poolSize} / ${settings.workerCount}`} />
						<Stat label="last run" value={formatMs(timings.runMs)} />
						<Stat label="per seq" value={formatMs(timings.perMs)} />
						<Stat label="version" value={version ?? "—"} />
					</dl>
				</section>

				<section class={PANEL}>
					{error ? (
						<p class="rounded-md border border-red-300 p-2.5 font-mono text-xs text-red-700">
							{error}
						</p>
					) : results ? (
						<Results results={results} />
					) : (
						<p class="text-sm text-zinc-500 italic">
							{busy ?? "Enter a sequence or FASTA to number."}
						</p>
					)}
				</section>
			</div>
		</div>
	);
}

function Field({
	id,
	label,
	hint,
	children,
}: {
	id?: string;
	label: string;
	hint?: string;
	children: JSX.Element;
}) {
	return (
		<div>
			{id ? (
				<label for={id} class={FIELD_LABEL}>
					{label}
				</label>
			) : (
				<p class={FIELD_LABEL}>{label}</p>
			)}
			{children}
			{hint ? <p class="mt-1 text-xs text-zinc-500">{hint}</p> : null}
		</div>
	);
}

function Checkbox({
	label,
	checked,
	onChange,
}: {
	label: string;
	checked: boolean;
	onChange: (checked: boolean) => void;
}) {
	return (
		<label class="flex items-center gap-1.5 text-sm">
			<input
				type="checkbox"
				checked={checked}
				onChange={(event) => onChange(event.currentTarget.checked)}
			/>
			{label}
		</label>
	);
}

function Stat({ label, value }: { label: string; value: string }) {
	return (
		<div class="flex gap-1.5">
			<dt>{label}</dt>
			<dd class="font-medium tabular-nums text-zinc-900">{value}</dd>
		</div>
	);
}

function requestOptions(settings: Settings): NumberingOptions {
	const options: NumberingOptions = {
		assignGermline: settings.assignGermline,
		alternativeHitCount: settings.alternativeHits ? 3 : 0,
	};
	if (settings.chains.length) options.allowedChains = settings.chains;
	if (settings.species.length) options.allowedSpecies = settings.species;
	const minBitScore = Number(settings.minBitScore);
	if (Number.isFinite(minBitScore)) options.minBitScore = minBitScore;
	return options;
}

function formatMs(value: number | undefined): string {
	return value === undefined ? "—" : `${value.toFixed(2)} ms`;
}
