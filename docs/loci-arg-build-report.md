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
