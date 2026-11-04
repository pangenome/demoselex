# Loci-only ARG pipeline — build report

*Branch `loci-arg-pipeline`, commit of 2026-10-21. This document explains what was
built, why, and how it was validated. Companion code: `analysis/g1k_loci_build.py`,
`analysis/build_input_tables.py`, `scripts/download_g1k_vcf.sh`,
`scripts/slurm_loci_arg.sh`, and the modified `analysis/g1k_allele_analyzer.py`.*

## 1. Context and goal

The diplotype analyzer (`analysis/g1k_allele_analyzer.py`) needs the dated 1000G ARG
(`data/derived/g1k/chr<N>_<start>_<end>.dated.trees`) to turn allele labels into
statements about selection vs demography. That ARG was previously built whole-genome on
another machine and is not in the repo. We rebuild it here, but only for the loci we
actually analyze: the 265 SV loci of `svs_refined.hgsvcv3.bed` (GRCh38) touch
**167 of ~2,885 ten-Mb chunks (~6% of the genome)** — so a loci-only build saves 94%
of the compute, at the cost of reduced inference context (see §3).

## 2. What was built

### 2.1 Per-locus dated chunks (`analysis/g1k_loci_build.py`)

For each 10 Mb grid chunk that a locus touches (plus the LCT validation chunk
chr2:135–136 Mb, included in every plan):

1. **stream** the 1000G high-coverage phased panel VCF (3,202 samples, GRCh38) for a
   window = chunk ± 5 Mb flank, keeping the exact site filters of `g1k_phased.py`
   (SNVs only, panel AF with min(AF,1−AF) ≥ 0.05, duplicate positions skipped);
2. **tsinfer** the window (flank = inference context);
3. carve the **central 10 Mb**, `keep_intervals().trim().simplify()` — precisely the
   operations `g1k_phased.date()` applies to a full-chromosome ARG;
4. **tsdate** `inside_outside`, µ=1.25e-8, Ne=10,000, generations — identical parameters;
5. dump under the canonical name `chr<N>_<start>_<end>.dated.trees` with positions
   relative to the chunk start, so `chunk_for()` in the analyzer works unchanged.

Engineering details that mattered in practice:

- **Resumability**: a finalised `SampleData` for a chunk is reused; a partial one is
  re-streamed. Completed chunks are skipped by existence.
- **htslib error detection**: stderr is captured during streaming; any `[E::` line
  (e.g. `bgzf_uncompress` failures) deletes the SampleData and fails the task so the
  scheduler retries, instead of silently producing a truncated ARG. This was added
  after remote streaming from EBI produced ~30% corrupted windows under 72-way
  concurrency (some "completing" with a fraction of their sites).
- **Local VCF preference**: if `data/derived/g1k/vcf/chrN…vcf.gz` exists it is used
  instead of the HTTPS URL. All 22 chromosomes were downloaded once (27 GB,
  size-verified, segmented + resumable — `scripts/download_g1k_vcf.sh`) after which
  streaming is fast and reliable.

**Faithfulness caveat**: the main build gives tsinfer whole-chromosome context; this
build gives it 20 Mb. The LCT selftest (§4) shows the difference is small at the
statistics we read.

### 2.2 Environment

A conda env at `.venv` with a self-built htslib 1.24 (the bioconda build lacks libcurl;
the pip cyvcf2 wheel bundles OpenSSL 1.0.2k from 2017 which cannot verify modern
certificates). This is what enables remote tabix access at all; the full-chromosome
downloads then made it moot for production runs.

### 2.3 The two input tables (`analysis/build_input_tables.py`)

Builds the INPUT CONTRACT files from `docs/loci-analysis-note.md`:

**`alleles_clusters.tsv`** — one row per (locus, sample), full grid, `NA` for uncalled:

| column | source |
|---|---|
| `chrom, pos, start, end` | `svs_refined.hgsvcv3.bed` (pos = interval midpoint; start/end = locus span) |
| moderns (`HG/NA…`, all 3,202) | **read-based genotypes** `geno_3202_fa/geno_3202_master.tsv.gz` (likegt, rank-1 calls) |
| ancients (1,628 samples) | `aDNA_real_1x_geno/results/*.geno1.tsv` (likegt on aDNA) |
| `allele1, allele2` | **cluster labels** from the MDS/DBSCAN haplotype clustering (`cosi_like/mds/paths_<locus>.tsv`) |
| `path1, path2` | original path keys (`sample#copy#contig`, the `:start-end` graph interval stripped — it is not a stable identity) |
| `cluster_k, silhouette, noise_frac` | per-locus clustering QC from `mds_summary.tsv` |

Label scheme: per locus, DBSCAN clusters are numbered **by descending size (1 = the
common class)**; **noise paths get private labels after k** so nothing is dropped
(see rationale below). `allele_labels_clusters.tsv` maps label → path key → cluster id.

Why cluster labels (this is the important design decision):

- With raw path-key labels, every modern's paths are private (their own assemblies):
  each label has exactly one ARG-cohort carrier, and the analyzer's `n ≥ 3` filter
  drops everything. The analysis is vacuous.
- Cluster labels define **structural haplotype classes** with many carriers, shared
  between moderns and ancients through the paths they were genotyped onto — which is
  the operational form of the "allelic haplotype sharing is the observable" argument
  in the methods note.
- Noise→NA was rejected with data: it would delete 15.2% of allele instances
  genome-wide and **78–86% at the most interesting loci** (LPA 78%, KATNAL2 86%,
  KIR2DL1 80%) — the hyper-diverse loci where DBSCAN's noise is likely real biology.
  Private labels keep those calls; they simply sit out of class-level statistics.

**`dates.tsv`** — `sample, date_bp, date_sd, source, culture` from
`old_dups_metadata-V1_dated.tsv` (date_mean_BP, date_sd_BP, date_method,
cultural_genetic_context); moderns at date 0. 5,162 dated ancients + moderns.

Known caveats (also printed by the builder): the GYP cluster has no modern graph-truth
rows when the identity table is used (per-gene CN cannot reconstruct the merged-locus
diplotype) — moot with the read-based 3202 genotyping, which genotypes the merged
locus; `--min-qv` is available but defaults to 0; the anchor `pos` is the SV-interval
midpoint (see §2.4).

### 2.4 Analyzer upgrades (`analysis/g1k_allele_analyzer.py`)

- **Ingestion hardened**: separator sniffed (TSV or CSV), `NA`/missing alleles →
  uncalled (`-1`), and `freq` divides by **called** haplotypes.
- **Multi-anchor (`--anchors K`)**: a structural class spans many local trees across a
  megabase-scale locus, and no single midpoint tree represents it. With `--anchors K`
  the analyzer samples K trees evenly across the locus span (using the table's
  `start`/`end`, clamped into the chunk) and aggregates per allele: median `age` plus
  `age_min`/`age_max`, median `age_pct`, `n_backgrounds`, `dom_share`, `dom_freq`,
  `background_n`, and `n_anchors` (how many anchors yielded `n ≥ 3`).
  Without `--anchors` (or without start/end columns) behaviour is the previous
  single-anchor mode.

## 3. Infrastructure

Chunks are independent single-threaded jobs: `scripts/slurm_loci_arg.sh` runs them as a
SLURM array on the `workers` partition (72 concurrent, 12 GB/task), reading chunk
coordinates from `data/derived/g1k/slurm_chunks.txt`. Per-chunk wall time ≈ 15 min
stream (local VCF on MooseFS) + 1–12 min tsinfer + 14–23 min tsdate ≈ 45–55 min;
159 chunks ≈ 2.5–3 h on 6 nodes.

## 4. Validation

**LCT selftest** (windowed rebuild vs the documented full-chromosome calibration,
EUR, tau=200):

| quantity | documented | this build |
|---|---|---|
| derived freq | 0.536 | 0.5355 |
| derived age | 495 gen | 516.5 |
| derived backgrounds / dom_freq | 1 / 1.00 | **1 / 1.00** |
| ancestral age / pct / backgrounds | 6,510 / 93.3 / 3 | 6,402.6 / 92.9 / **3** |

**chr22 pilot** (4 loci, all 3,202 moderns EUR-restricted to 633, `--anchors 9`):

- The full chain ran end to end: read genotypes → cluster classes → local ARG →
  per-class freq/age/backgrounds.
- Multi-anchor exposed the recombination-scrambling it was designed for: per-class
  ages vary hugely across a locus's trees (e.g. CYP2D7 allele 1: median 8,814 gen,
  range 4,203–34,078), where the single-anchor mode had returned one identical MRCA
  for every class at a locus.
- Read-based moderns reshape the spectrum vs the 233 graph-truth samples: POTEH
  changes from an apparent 96% single class to 52%/25%/5% — the graph-panel samples
  flatten structure.

## 5. Genome-wide analysis (all 256 loci, 8 configurations)

`scripts/run_genome_wide.sh` ran the analyzer over the full table: global mixed
cohort at tau 200/500/1000, and per-superpopulation (AFR/AMR/EAS/EUR/SAS) at
tau 200, all with `--anchors 9` (~32,000 class rows).

Findings:

- **Tau x cohort interaction (measured)**: at tau 200–500 in the global cohort,
  backgrounds are fragmentary (avg dom_share 0.37→0.49 across the tau sweep;
  dominant-background size 491→1,012 haps). Background statistics only become
  informative when tau is scaled to where the cohort's structure lives.
- **Structural classes are generally not genealogical clades**: ~99% of classes
  have within-locus age percentile > 90 - a sequence-similarity class over a
  megabase locus is recombination-scrambled relative to any single local tree.
  The multi-anchor `age_min` is the candidate-discovery axis; the
  `age_min`-`age_max` spread is the scrambling measure.
- 182 classes had `age_min` <= 1,000 gen with n >= 10 - but see §6.

## 6. Anchor-profile validation: age-based candidates are regional effects

`analysis/anchor_profile.py` profiled the 57 frequency-band candidates at 21
anchors per locus. Every locus showed the same youth pattern across ALL its
classes (identical profiles), with age percentile = 100 at the young anchors too:
the youth signal belongs to the **regional SNP genealogy** (locus-wide shallow
coalescence), not to any structural class. Loci where this regional shallowness
appears (ANKRD36C region, POLR2J2, ARHGEF35) coincide with the repo's
independent SNP-level screen hits. Classes with age 1-15 gen and n=3-6 are
trivially explained by close relatives in the 3,202 panel (the 698 HGDP samples
are a related subset). Conclusion: no class-specific young structural haplotypes;
the static snapshot cannot be the selection test (see §7) - consistent with the
method's own design (selection on standing variation is invisible to age filters).

## 7. The temporal layer (`analysis/temporal_traj.py`)

Per (locus, class): equal-count time bins of the 1,627 dated called ancients
(mean bin dates back to ~12.7 kya) with the modern cohort as the present-day
endpoint; rise and biggest adjacent jump; punctuated/gradual shape; diffusion
drift null (z at Ne 5,000 and 10,000); **date-permutation null** (bin-assignment
permutation; a synthetic rise/flat/fall selftest is built in and caught a real
permutation bug during development). `--modern-pop` restricts moderns for
like-for-like comparison; `--anc-group-file/--group` stratify by country.

Result funnel (EUR-restricted moderns):

  200 RISING (global moderns)  ->  58 (EUR moderns; the rest were
  population-composition change)  ->  29 at tractable loci with real ancient
  support  ->  **13 ROBUST** (rise repeats within >= half of 7 country strata:
  Hungary, Russia, Belgium, UK, Italy, Sweden, Denmark).

The 13 robust classes: TARP.1, CYP2D7.1, KLRC1.1, OR2L2.1, CFHR1.1, ZNF146.2,
UPF2.1, ABCA13.1, ADGRG7.4, HLA-A.5, HLA-A.19, THSD7B.2, OR2J2.17 - concentrated
in immune and metabolic duplicated loci (CYP2D6/7, KLRC1/NK receptor cluster,
CFHR1 complement CNV, HLA), rising in 3-7 of 7 countries each.

## 8. Hard-to-genotype loci: a QV-free metric and the AMY1A sensitivity floor

Locus tractability is measurable without trusting likegt's graph_qv: the overlap
between the ancient and modern **class-frequency distributions** at a locus.
33/256 loci fall below 0.5 overlap (POTEH 0.04, POTEA 0.12, TAF11L11 0.17,
KATNAL2 0.28, KIR2DL1 0.32, AMY1A 0.42...). At such loci, 1x ancient calls
scatter onto noise/private labels, deflating ancient frequencies of the big
classes - trajectories there are uninterpretable (POTEH's apparent punctuated
rise is entirely this artifact, visible as the dashed 'other' line in the
figures).

AMY1A is the calibration case: Bolognini et al. (Nature 634:617-625, 2024;
doi:10.1038/s41586-024-07911-1) demonstrated a strong, real rise of
duplication-containing amylase haplotypes over 12 kyr in West Eurasians. Our
pipeline sees the direction (class 1: 13% -> 34%) but the permutation null buries
it (p = 0.98) because the ancient calls at AMY1A are too scattered. That is the
method's quantified detection floor - and their result also shows why: AMY1
dup/del mutation rates up to ~10,000x the SNP rate mean class frequencies are
shaped by mutation pressure on a state, not clonal descent, so (a) clade-age
thinking and (b) a drift-only null are both wrong at recurrent-CNV loci. Several
robust candidates (CYP2D7, KLRC1, CFHR1) are documented recurrent-CNV loci and
carry the same caveat: their rises are consistent with selection on structural
states under high mutation turnover, and a mutation-aware null (their Fig. 5d-g
equivalent: WF-diffusion/ApproxWF selection-coefficient estimation with a CNV
rate parameter) is the required next step before calling them selection in print.

## 9. The tau-free breadth readout: dom_share(tau) profiles and crossover tau*

Per the updated `docs/loci-analysis-note.md` (2026-10-22), breadth should not be
read off one tau but reported as a `dom_share(tau)` profile summarized by the
crossover depth tau* at which the dominant background comes to hold >= 50% of
carriers. `analysis/tau_star_profiles.py` computes this per class from the
existing multi-tau runs (three profile points at tau 200/500/1000, log-
interpolated crossover), with the by-construction guard: confinement is only
genealogically meaningful when tau* < clade age (at tau >= age all carriers sit
in one background by construction; exactly 1 class genome-wide trips the guard).
Output: `genome_wide_erikarg/tau_star_profiles.csv` (8,032 classes: 2,048
confined at tau<=200, 2,150 with a tau* in (200,1000], 3,833 not confined
within 1,000).

Applied to the 13 ROBUST risers, the tau* readout separates them into two
biologically distinct regimes:

- **distributed risers (8)** - TARP.1, CYP2D7.1, KLRC1.1, OR2L2.1, CFHR1.1,
  ZNF146.2, UPF2.1, ABCA13.1: flat, low profiles (dom_share 0.15-0.41 at every
  tau); carriers never concentrate into one background at any depth. Standing
  structural variation whose *frequency* rose in Europe, spread across
  ancestries - the regime the original decision table has no row for, and the
  one the mutation-aware null must cover at recurrent-CNV loci
  (CYP2D7, KLRC1, CFHR1).
- **confined risers (3)** - ADGRG7.4 and THSD7B.2 (dom_share >= 0.59 already at
  tau=200), HLA-A.19 (tau* ~ 357): old but genealogically concentrated classes
  that recently became more frequent.

Figure: `plots/tau_star.png` (`analysis/plot_tau_star.py`) - (A) genome-wide
regime map (clade age vs dom_share at tau=1000, candidates highlighted against
the 0.5 threshold); (B) per-candidate profiles with tau* crossings.

## 10. Outputs and figures

All under `data/derived/old_dups/` (gitignored, results are not in the repo):
- input tables: `alleles_clusters.tsv`, `dates.tsv`, `allele_labels_clusters.tsv`,
  `anc_country.tsv`
- static analysis: `genome_wide/{global_tau200,global_tau500,global_tau1000,
  AFR_tau200,...}/summary.csv` (our windowed ARG) and the identical
  `genome_wide_erikarg/` set on the delivered full-genome ARG, plus
  `genome_wide_erikarg/tau_star_profiles.csv` (Section 9)
- trajectories: `genome_wide/trajectories.csv` (global),
  `trajectories_EUR.csv`, `genome_wide/by_country/<Country>.csv`
- candidates: `rising_candidates_annotated.csv` (58, annotated),
  `rising_robustness.csv` (29 -> 13 ROBUST with per-country concordance)
- figures (`plots/`, `analysis/plot_results.py`, `plot_tau_star.py`,
  `plot_atlas_v2.py`): `trajectories_robust.png`, `top6_panels.png` (drift
  ribbons), `rise_overview.png`, `tau_star.png` (Section 9),
  `locus_composition.png` / `locus_composition_all.png` (per-locus class
  composition, Bolognini Fig 5C style), and the atlases: v1
  `all_loci_composition.pdf` and v2 `all_loci_composition_v2.pdf` - per panel,
  composition through time with ROBUST classes starred by tau* regime (red =
  distributed, blue = confined), per-locus confinement summary
  `confined k/N tau*`, and [HARD] flags.

## 11. Status and next steps

- ARG complete (167/167 chunks); pipeline, tests and this report on branch
  `loci-arg-pipeline` (PR to pangenome/demoselex).
- Interpretation reminders from `docs/loci-analysis-note.md` still apply: age
  first, breadth only when age > tau; no selection coefficients yet; `age_pct`
  is descriptive. Breadth is now read tau-free via the dom_share(tau) profile
  and crossover tau* (Section 9), per the updated methods note.
- Next: selection-coefficient estimation per robust candidate with a
  mutation-aware null (recurrent-CNV loci); timing of the punctuated jumps
  against the Neolithic/Steppe transitions; copy-number-state class definitions
  at hard loci (Bolognini-style architectures).

---

## Appendix A. Analysis architecture: why every cohort separation exists

Two core facts force the structure:

1. **Tree statistics are computed inside the cohort you give them.** The analyzer
   restricts the ARG to the samples present in the table and computes everything
   on that subtree: the *age* of a class is the MRCA age of its carriers in that
   cohort's tree (not an absolute number), and a *background* is a group of
   haplotypes that coalesce within tau generations *among those people*.
2. **The samples are not interchangeable.** Moderns are global; the ancients are
   mostly West Eurasian; the modern panel contains close relatives (the 698 HGDP
   samples are a related subset); Europe's population changed by layered
   migration (hunter-gatherers -> Neolithic farmers -> Steppe). Any statistic
   computed on a mixture inherits the mixture.

### A.1 What a cohort restriction actually does: the induced subtree

The restriction is explicit in the analyzer
(`analysis/g1k_allele_analyzer.py`, per locus):

```python
ts = tskit.load(str(chunk)).simplify(samples=haps_of_the_cohort)
```

`tskit.simplify()` builds the **induced subtree**: every lineage not ancestral
 to a retained haplotype is pruned, so the analysis runs on the genealogy of
the cohort present in the table - and nothing else. Three consequences follow,
and each one is why a corresponding cohort separation exists:

**(a) Node dates are inherited, not recomputed - but class ages still change
with the cohort.** The ARG is dated once, globally; simplify keeps the retained
nodes' absolute dates. The MRCA node of two haplotypes is the same dated node in
the full tree and in any restriction that keeps both. What *does* change is
which carriers define the class: a class is "all carriers of the label present
in this cohort", so a globally distributed class shrinks to its European
carriers in an EUR run, and their MRCA can be far younger than the worldwide
one. Age is therefore a property of the *class x cohort* pair, never of the
label alone. (This is the mechanism behind Section 6: age-based candidates were
regional shallow-genealogy effects - the youth was in the regional carrier set,
not in any structural haplotype.)

**(b) Backgrounds live inside the subtree - the reason the per-superpopulation
layer exists.** A background is defined by coalescences *among the cohort's
lineages* within tau. After simplify, non-cohort samples cannot connect two
cohort haplotypes. In the global tree at tau=200, two haplotypes from different
continents never coalesce (their meeting point is tens of thousands of
generations deep), so almost every haplotype is its own background and breadth
is vacuous. Restrict to one superpopulation and the pruned tree contains only
that population's lineages, whose shallow within-population coalescences now
define the backgrounds - dom_share and n_backgrounds become meaningful at
small tau. This is the measured Layer-2 confound: global dom_share 0.37 -> 0.49
and dominant-background size 491 -> 1,012 haps across the tau sweep.

**(c) The frequency denominator follows the cohort.** `freq` divides by the
cohort's *called* haplotypes, so the same label legitimately has different
frequencies in global, EUR and per-country runs - a rise seen by comparing a
restricted ancient series against a wider modern panel is partly composition,
not biology. This is why Layer 4 restricts the moderns to EUR when the ancient
series is European (200 global RISING -> 58 after matching), and why Layer 5
re-checks the survivors within single countries.

**Boundary note: the ancients are never in the tree.** The ARG contains only the
3,202 moderns; every "cohort restriction" above is over moderns. The aDNA side
enters exclusively through the temporal layer, which uses the genotype table
plus dates and no tree at all. The two halves of the analysis (static genealogy,
temporal trajectories) are therefore cleanly separable substrates, and the
ROBUST verdicts rest on the temporal half.

```
LAYER 0  genotype labels         class = structural haplotype (MDS/DBSCAN clusters
                                 over graph paths; noise keeps private labels)
LAYER 1  the ARG                  dated genealogy of the 3,202 moderns (the substrate;
                                 ages are substrate-dependent - compare percentiles,
                                 never absolute ages, across builds)
LAYER 2  static per-class stats   freq, clade age + age_pct, backgrounds(tau)
                                 RUNS: global x tau{200,500,1000} and
                                 per-superpop x tau200  = the "8 configurations"
LAYER 3  candidate validation    per-anchor youth profiles
                                 (all age-based candidates were regional effects)
LAYER 4  temporal trajectories   class frequency through ancient time
                                 RUNS: global moderns -> EUR-restricted moderns
                                 (cohort matching; 200 -> 58 RISING)
LAYER 5  within-population       per-country trajectories, 7 European countries
         replication             (58 -> 29 -> 13 ROBUST)
```

Each separation controls one specific confound, measured on this data:

| separation | confound it kills | measured effect here |
|---|---|---|
| per-superpopulation runs (Layer 2) | backgrounds need coalescences below tau; in the global tree at tau=200 they are ~all singletons, so breadth is vacuous | dom_share 0.37 -> 0.49 and dominant-background size 491 -> 1,012 across tau 200/500/1000 global; only population-restricted runs make tau=200 meaningful |
| EUR-restricted moderns vs ancients (Layer 4) | comparing *European ancients* to *global moderns* turns continent composition into fake "rises" | 200 global RISING -> 58 after cohort matching (the rest were turnover) |
| per-country stratification (Layer 5) | within-Europe migration (Steppe/Neolithic replacement) mimics selection in a pooled series | 29 -> 13 ROBUST (pooled-only classes had per-country rises of ~0.00-0.03) |
| n >= 10 / min-class 20 | close relatives in the panel share private calls | age-1-15-generation classes with n=3-6 are relatives, not biology |
| tractable-locus filter (ancient/modern class-distribution overlap >= 0.5) | 1x ancient calls scattering onto noise labels deflates ancient frequencies of big classes | 33/256 loci excluded (POTEH 0.04, AMY1A 0.42); the AMY1A case (Bolognini et al. 2024) is the quantified detection floor |

A single run-everything-together analysis would still produce numbers for all of
these questions - but each number would silently be answering a different, wrong
question. This is the same reasoning that motivates tau-free verdicts in the
updated methods note: no statistic should depend on a dial (cohort, tau) that can
be tuned until the answer looks convincing.

What tau is and is not used for (asked frequently):

- tau enters exactly one function, `backgrounds()`: walk up the local tree from a
  sample while parent age <= tau; the last node reached is that haplotype's
  background. From it: n_backgrounds, dom_share, dom_freq, background_n.
- tau is NOT used for freq, n, age, age_pct - and not anywhere in the temporal
  layer (drift and permutation nulls use real dates, Ne and generations).
- The 13 ROBUST candidates are therefore tau-independent: the tau sweep only
  affects the descriptive background statistics, and the trajectory verdicts
  never touch tau.

## Appendix B. Cross-build validation against the delivered full-genome ARG

The full-genome ARG bundle (E. Garrison, `/moosefs/erikg/demoselex-arg/`: 294 whole-
genome chunks, its own analyzer copy with bundle path resolution, cohort file with
identical sample order) was used to re-run all eight configurations
(`scripts/run_genome_wide_erikarg.sh`, output `data/derived/old_dups/genome_wide_erikarg/`).
His analyzer selftest reproduces the documented LCT calibration to the decimal
(derived 495.3 gen / pct 20.0 / 1 background; ancestral 6,510.3 / 93.3 / 3).

Genome-wide agreement, windowed rebuild vs full-genome build (all 8 configs):

- carrier counts: identical for every common class (n is table-derived)
- median-anchor age: geometric-mean ratio 1.006-1.014
- age percentile: within 10 points for ~100% of classes
- n_backgrounds: within +-1 for 98-99% of classes

A single-anchor comparison on four chr22 loci shows large position-specific age
differences (ratio 0.25-1.66 by locus); the 9-anchor median removes them. The
multi-anchor aggregation is therefore robust both across positions within a locus
and across ARG builds. All substantive conclusions hold identically on both
substrates; absolute ages remain build-dependent (compare percentiles across
builds, never ages).

## Appendix C. Genome-wide exact clade test (the pilot's Q2, answered for every class)

The delivered bundle ships `sv_phase_pilot.py`: an exact single-origin test for
unphased SV genotypes on the dated ARG (Stage A phase-free ancestry enrichment;
Stage B exact clade-consistency with induced phase; Stage C certificate over all
admissible phases). We ran it verbatim (his defaults, his ARG) on ALL 256 loci and all 8,032 structural classes - not only the candidates - via `scripts/run_sv_pilot_all.sh` + `analysis/sv_pilot_scan.py` (output: `/moosefs/laurap/old_dups/demoselex_exp/pilot_run/pilot_genome.csv`, 23,920 class x position rows - outside the repo, alongside the per-chromosome tables and logs in the same directory).

**The base rate (Q2 asked genome-wide):**

| certificate | share of evaluations |
|---|---|
| NO-CLADE | 23,803 (99.5%) |
| SINGLE-ORIGIN-CONSISTENT | 114 (0.5%) |
| NO-HET | 3 |

5 classes (0.1% of 8,032) are single-origin at every evaluated position
(TTC28 x2, MAP3K2, SPATA5, MEI4); 89 (1.1%) at >= 1 position. A megabase
structural class is essentially never an exact genealogical clade - the exact
test confirms, at ARG resolution, what the multi-anchor age percentiles said
(Section 5): **the class label is not a lineage**. Confinement (tau*) and exact
clade-consistency are different statements: 64/4,198 tau*-confined classes are
single-origin somewhere (1.5%) vs 25/3,834 distributed (0.65%) - a weak positive
association only. dom_share ~ 0.5-0.65 means majority concentration, not
cladal descent.

**The 13 ROBUST candidates under Q2:** every one is NO-CLADE (fit errors 60 -
1,830 genotypes; n_origins 4 - 5,158). None is an exact clade - including the
tau*-confined three (ADGRG7.4, THSD7B.2, HLA-A.19). What separates the
candidates instead is Stage A (phase-free ancestry enrichment): 10/13 have
carriers significantly ancestry-conriched (permutation p at the 5,000-draw
floor, p = 0.0002) - CYP2D7.1, KLRC1.1, OR2L2.1, CFHR1.1, ZNF146.2, ABCA13.1,
ADGRG7.4, THSD7B.2, HLA-A.19, OR2J2.17 - while TARP.1 (p = 0.88), UPF2.1
(p = 0.76) and HLA-A.5 (p = 1.00) rise without ancestry concentration at all:
genuinely cosmopolitan classes.

**Operational notes from the genome-wide run** (worth reporting upstream):

- ancients are dropped from the input (the ARG contains only the 3,202
  moderns), matching the static analyzer's behavior;
- 6 loci straddle two chunk files and are degraded to midpoint-only;
- TPTE2 (chr13:19.3-19.6 Mb, pericentromeric) has NO tree coverage on either
  build: the chunk named chr13_10000000_20000000 is actually trimmed to
  13,999,701 bp. The pilot crashes loudly there (`ts.at` out of bounds via the
  filename-trusting chunk lookup); the static analyzer had been silently
  skipping the locus all along (zero TPTE2 rows in every summary.csv). The
  pilot's loud failure is the better behavior.

**Interpretation.** This is the genome-scale confirmation of the lineage-scan
note's premise: the strongest available statement about a *label* is
descriptive. The candidates' rises are real as frequency changes of structural
classes (the temporal layer), but no candidate is a young clade - the
selection question moves to erik's lineage-resolved program (allele placed on
backgrounds, realized-genealogy nulls), which is the correct instrument for
what happens next. The Stage A split (ancestry-enriched vs cosmopolitan
risers) is the first axis of that next step and is already computed.

## Appendix D. Structural phase: likegt labels onto VCF/ARG haplotype space

Likegt's hap1/hap2 assignments are per-(sample, locus) arbitrary numbering
(read-based, no genome-wide phasing input) - measured: S vs SW is ~52/48 over
280,524 het (sample, locus) rows. HPRC assembly #1/#2 is equally unrelated to
VCF hap order (all 223 cohort donors mixed ~50/50 across loci). Neither
numbering can be assumed; per-locus linkage is the only valid route.

**Method** (`analysis/phase_link.py`): each assigned path is a real HPRC donor
haplotype segment. Its sequence is reconstructed from the per-locus GFA,
aligned to the GRCh38 window (minimap2 -x asm20 -c), and read at every site
where the sample's two ARG-phased haplotypes differ; the 2x2 match gives the
orientation, with concordance = 1 - (donor-sample divergence at het sites).

**Results genome-wide** (280,524 het rows, 256 loci):
- 83.4% decided (mean best-diagonal concordance 0.888; mean gap 0.776 -
  usually decisive), 10.5% AMBIG, 6.1% too-few-sites, TPTE2 no coverage
- 16,101 rows are self-assignments (223 cohort members are HPRC donors):
  assembly-vs-VCF phase agreement mean 0.959, 51% perfect at all het sites

**Validation** (`analysis/phase_link_trio_check.py`): phase-aware Mendelian
transmissions - child's paternal ARG hap identified by SNP match, the
structural cluster must ride the transmitted parent hap. 57,313 checks
(602 trios, 227 loci): **86.9% phase-consistent**. Mismatches: 6.4%
phase-swaps, 6.7% de-novo/call errors. Hard-list loci median 0.649 vs
0.938 elsewhere - independent confirmation of the QV-free hardness metric by
Mendelian logic. Candidate loci: THSD7B 100%, UPF2 100%, OR2J2 98.9%,
CFHR1 97.3%, ADGRG7 95.9%, HLA-A 96.2%, KLRC1 96.0%, ABCA13 95.4%;
weaker: ZNF146 64.1% (hard), TARP 74.7%, CYP2D7 79.6%, OR2L2 84.8%.

**Combined VCFs** (`analysis/build_combined_vcf.py`, runner
`scripts/run_combined_vcf_all.sh`, output
`data/derived/old_dups/combined_vcf/<LOCUS>.combined.vcf.gz`, tabix-indexed):
per locus, all phased 1000G panel records in [lo-5kb, hi+5kb] unchanged, plus
one biallelic symbolic record per structural class at the locus midpoint
(REF = reference base, ALT = <CLS:k>), GT phased in VCF hap space from the
link ('|' where decided, '/' where not - never guessed; '.' for uncalled).
256 loci, 9,865 class records, 807 MB. Any standard tool can now treat the
structural classes as ordinary phased variants (LD, haplotype matching,
imputation checks) - the direct input to the lineage program's within-/
between-background decomposition.

**Outputs of this workstream** (all gitignored results, outside the repo tree
at `/moosefs/laurap/old_dups/demoselex_exp/`):
- `phase_link_run/<LOCUS>.csv` + merged `phase_link_genome.csv`: per
  (sample, locus) orientation and phased cluster labels (`arg_cluster_hap0/1`)
- `phase_link_run/trio_check.csv`: the 57,313 Mendelian checks with per-locus
  agreement rates (use as the per-locus phase-reliability weight)
- `data/derived/old_dups/combined_vcf/`: the 256 combined VCFs

Cross-references: this is the WP2 donor-side ingredient of
`docs/next-steps-arg-sv-phasing.md` (the 223-donor assembly-to-VCF-hap map is
`phase_link_genome.csv` itself), and the input substrate for Q1 of
`docs/lineage-resolved-selection.md` (within- vs between-background rise),
which can now be computed for every class genome-wide.

## Appendix E. Class construction in detail: from likegt paths to cluster labels

How the integer class labels in `alleles_clusters.tsv` / `allele_labels_clusters.tsv`
were produced: the exact inputs, the command, and every design decision.

### E.1 The five input files

| input | what it is | upstream |
|---|---|---|
| `svs_refined.hgsvcv3.bed` | the 256 locus spans (GRCh38): chrom, start, end, gene - defines locus identity, midpoint `pos`, and the `start`/`end` columns | refined SV/dup call set |
| `geno_3202_fa/geno_3202_master.tsv.gz` | modern genotypes, one row per (sample, locus): sample, locus, hap1, hap2, cosine, graph_qv, reads_realigned, nodes_used; hap1/hap2 are **assigned donor paths** like `HG00621#2#CM087873.1:38485010-38771168` | `likegt geno` (SLURM array): each of the 3,202 moderns' reads realigned to the per-locus path catalog, best pair by coverage cosine |
| `aDNA_real_1x_geno/results/*.geno1.tsv` | ancients: one file per sample (1,628), same schema - likegt on the 1x reads (pseudo-haploid, effectively one haplotype draw; the table still reports the best hap1/hap2 pair) | the aDNA genotyping run |
| `old_dups_metadata-V1_dated.tsv` | date_mean_BP, date_sd_BP, date_method, cultural_genetic_context per ancient (feeds `dates.tsv`, not the labels) | AADR-derived dating |
| `cosi_like/mds/paths_<locus>.tsv` (+ `mds_summary.tsv`) | `idx, path, cluster` for every catalog path - the MDS/DBSCAN clustering output; summary QC: k, silhouette, noise_frac per locus | odgi cosine distances -> exact MDS -> DBSCAN on the distance matrix (full methods: `cosi_like/mds/MDS_CLUSTERING_REPORT.md`) |

### E.2 The command

```bash
.venv/bin/python analysis/build_input_tables.py \
    --bed          /moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed \
    --modern-geno  /moosefs/raid5/lpignata/old_dups/geno_3202_fa/geno_3202_master.tsv.gz \
    --geno-dir     /moosefs/laurap/old_dups/aDNA_real_1x_geno/results \
    --metadata     /moosefs/laurap/old_dups/old_dups_metadata-V1_dated.tsv \
    --mds          /moosefs/laurap/old_dups/cosi_like/mds \
    --min-qv       0
```

`--min-qv 0` is deliberate: likegt's `graph_qv` is NOT trusted as a filter
(early design decision). Call difficulty is instead measured data-internally
by the ancient/modern class-distribution overlap (Section 8, the [HARD] flag).

### E.3 Steps (code references: `analysis/build_input_tables.py`)

1. **path_key normalization** - every likegt path ID loses its `:start-end`
   suffix (`HG00621#2#CM087873.1:38485010-38771168` ->
   `HG00621#2#CM087873.1`): the interval is the *graph* interval, not a
   stable identity.
2. **Load genotypes**: 3,202 moderns from the master; 1,628 ancients from
   the per-sample files (3-retry read for moosefs transient errors; unreadable
   -> all-NA rows). Loci = union, required to exist in the BED.
3. **Load cluster maps**: per locus, path_key -> cluster_id (raw MDS/DBSCAN
   id); QC from mds_summary.tsv.
4. **Label map per locus** (`cluster_labels`): over the union of *assigned*
   path keys (moderns UNION ancients) - count distinct assigned paths per
   MDS cluster; renumber clusters by descending count -> **class labels 1..k**
   (ties by cluster id); then every noise path (cluster -1) and any path
   missing from the map gets a **private label k+1, k+2, ...** in sorted
   order. Never merged, never NA.
5. **Write**:
   - `allele_labels_clusters.tsv`: locus, path_key, **cluster_id (raw MDS
     id, -1 for noise)**, **allele (renumbered class)** - the label map
   - `alleles_clusters.tsv`: the full grid, 256 loci x (3,202 + 1,628) =
     1,236,481 rows: chrom, pos(midpoint), sample, allele1, allele2, path1,
     path2, cluster_k, silhouette, noise_frac, start, end; NA where uncalled
   - `dates.tsv`: dates from the metadata; moderns = 0 BP

### E.4 The three orderings that are easy to conflate

Per locus there are three different "sizes", and none implies another:

1. **catalog paths** per MDS cluster (what the MDS plot shows: ~466 paths
   per locus distributed over k clusters + noise)
2. **assigned paths** per cluster (what the class renumbering actually
   ranks): a cluster with 465 catalog paths can lose the "class 1" slot to a
   smaller cluster with more assigned paths. This flips at 7 loci: GOLGA6D,
   LPA, KATNAL2, HLA-DQA1, CHST5, C2orf78, ANKRD36C
3. **carrier haplotypes** per class (what selection analyses see): a single
   popular donor path can carry hundreds of haplotypes (ZNF146 class 2:
   2,111 haps via 5 assigned paths of a 6-path cluster)

The reconciliation table `data/derived/old_dups/cluster_class_map.tsv`
(`analysis/cluster_class_map.py`, decoder) joins all three levels per locus:
MDS cluster -> catalog paths -> assigned paths -> class label -> carrier haps,
plus one row per private (noise-path) class.

### E.5 Design decisions and their reasons

- **Noise paths -> private labels, never NA**: dropping would delete every
  sample whose best match is an outlier path (~15% of calls genome-wide,
  78-86% at the hardest loci); merging would fabricate a heterogeneous
  pseudo-class. Private labels keep all calls. Cost: the class count inflates
  (8,032 "classes" genome-wide vs ~2-10 real clusters per locus); most
  classes are single-donor-path families with 1-40 carriers.
- **Class 1 is not guaranteed to be the biggest cluster** (E.4.2) - and class
  number is not comparable across loci; always join through
  `allele_labels_clusters.tsv` / the decoder.
- **Real clusters vs private classes stratify naturally** (cluster rows vs
  private rows in the decoder) - recommended for any downstream analysis,
  since single-path classes can look artificially lineage-like.

### E.6 Historical fossils (why the names are what they are)

The script was first run WITHOUT `--mds` (labels 1..k by sorted path key;
moderns from `gene_cn/all_copies_identity.tsv`, only 233 assembly donors).
Two pivots produced the current scheme: (a) moderns switched to the
read-based `geno_3202` (all 3,202 in the ARG cohort), because path-key labels
turned out to be almost all private - 1 modern carrier each, analysis
vacuous; (b) labels switched to MDS-cluster identity, which defines shared
structural families. The `_clusters` suffix in the output filenames and the
`suffix = '_clusters' if cl_of else ''` branch in the code are fossils of
that history.
