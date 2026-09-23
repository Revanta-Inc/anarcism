import type { DomainResult, Region, SequenceResult } from "./anarcism.ts";

const REGIONS: readonly Region[] = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"];

const FRAMEWORK = "bg-zinc-100 text-zinc-600";
const REGION_COLOR: Record<Region, string> = {
	FR1: FRAMEWORK,
	CDR1: "bg-amber-100 text-amber-800",
	FR2: FRAMEWORK,
	CDR2: "bg-emerald-100 text-emerald-800",
	FR3: FRAMEWORK,
	CDR3: "bg-rose-100 text-rose-800",
	FR4: FRAMEWORK,
};

const BADGE =
	"rounded-full border border-zinc-200 bg-zinc-100 px-2 py-0.5 text-xs text-zinc-600 tabular-nums";

export function Results({ results }: { results: SequenceResult[] }) {
	return (
		<div>
			<div class="flex flex-wrap gap-1.5">
				{REGIONS.map((region) => (
					<span key={region} class={`rounded px-2 py-0.5 text-xs ${REGION_COLOR[region]}`}>
						{region}
					</span>
				))}
			</div>

			{results.map((result, index) => (
				<div key={`${result.id}-${index}`}>
					<div class="mt-5 flex flex-wrap items-baseline gap-2 border-b border-zinc-200 pb-2">
						<h2 class="font-mono text-sm font-medium">{result.id || "(sequence)"}</h2>
						<p class="text-xs text-zinc-500">
							{result.normalizedSequence.length} aa · {result.domains.length}{" "}
							{result.domains.length === 1 ? "domain" : "domains"}
						</p>
					</div>

					{result.warnings.map((warning) => (
						<p key={warning} class="mt-2 font-mono text-xs text-red-700">
							{warning}
						</p>
					))}

					{result.domains.length === 0 ? (
						<p class="mt-2 text-sm text-zinc-500 italic">No variable domain detected.</p>
					) : (
						result.domains.map((domain) => <Domain key={domain.domainIndex} domain={domain} />)
					)}
				</div>
			))}
		</div>
	);
}

function Domain({ domain }: { domain: DomainResult }) {
	return (
		<div class="mt-4">
			<div class="flex flex-wrap gap-1.5">
				<span class="rounded-full bg-blue-600 px-2 py-0.5 text-xs font-semibold text-white">
					{domain.chainType}
				</span>
				<span class={BADGE}>{domain.receptorType}</span>
				<span class={BADGE}>{domain.species}</span>
				<span class={BADGE}>
					{domain.start}–{domain.end}
				</span>
				<span class={BADGE}>bit {domain.bitScore.toFixed(1)}</span>
				<span class={BADGE}>E {formatE(domain.eValue)}</span>
				{domain.germline ? (
					<span class={BADGE}>
						{domain.germline.species} · {domain.germline.vGene ?? "?"}{" "}
						{formatIdentity(domain.germline.vIdentity)} · {domain.germline.jGene ?? "?"}{" "}
						{formatIdentity(domain.germline.jIdentity)}
					</span>
				) : null}
			</div>

			<div class="mt-2.5 flex flex-wrap gap-px">
				{domain.numbering.map((residue) => (
					<div
						key={`${residue.position}${residue.insertionCode}`}
						class={`w-7 rounded py-0.5 text-center font-mono leading-tight ${REGION_COLOR[residue.region]}`}
						title={`${residue.region} · IMGT ${residue.position}${residue.insertionCode} · ${residue.aminoAcid}${
							residue.sequenceIndex >= 0
								? ` · input index ${residue.sequenceIndex}`
								: " · scheme gap"
						}`}
					>
						<span class="block text-[8px] tabular-nums opacity-70">
							{residue.position}
							{residue.insertionCode}
						</span>
						<span
							class={`block text-[13px] font-semibold ${residue.aminoAcid === "-" ? "opacity-30" : ""}`}
						>
							{residue.aminoAcid}
						</span>
					</div>
				))}
			</div>

			<p class="mt-2.5 rounded-md border border-zinc-200 bg-zinc-50 p-2 font-mono text-xs break-all text-zinc-500">
				{domain.paddedImgtAlignment}
			</p>

			{domain.alternativeHits.length > 0 ? (
				<details class="mt-2">
					<summary class="cursor-pointer text-xs text-zinc-500">
						{domain.alternativeHits.length} alternative hits
					</summary>
					<table class="mt-1.5 text-xs">
						<thead class="text-zinc-500">
							<tr>
								<th class="pr-3 text-left font-semibold">profile</th>
								<th class="pr-3 text-left font-semibold">chain</th>
								<th class="pr-3 text-left font-semibold">species</th>
								<th class="pr-3 text-left font-semibold">bit</th>
								<th class="text-left font-semibold">E</th>
							</tr>
						</thead>
						<tbody class="tabular-nums">
							{domain.alternativeHits.map((hit) => (
								<tr key={hit.profile}>
									<td class="pr-3">{hit.profile}</td>
									<td class="pr-3">{hit.chainType}</td>
									<td class="pr-3">{hit.species}</td>
									<td class="pr-3">{hit.bitScore.toFixed(1)}</td>
									<td>{formatE(hit.eValue)}</td>
								</tr>
							))}
						</tbody>
					</table>
				</details>
			) : null}
		</div>
	);
}

function formatE(value: number | undefined): string {
	if (value === undefined) return "—";
	if (value === 0) return "0";
	return value < 1e-4 ? value.toExponential(1) : value.toPrecision(3);
}

function formatIdentity(value: number | undefined): string {
	return value === undefined ? "" : `${(value * 100).toFixed(1)}%`;
}
