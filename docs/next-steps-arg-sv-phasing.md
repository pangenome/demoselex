# Next work: SV phasing onto the dated ARG, and the steps that depend on it

Status: proposed. Owner: chat agent. Supersedes nothing; extends `docs/MISSION.md` §17 outputs
without amending the objective.

This document says what to build next, in what order, why each step exists, and what would
count as it having failed. The ordering is not arbitrary: four of the steps are *gates*,
and quoting a number before its gate is passed is the failure mode this plan exists to
prevent.

---

## 1. Where we are

Validated and on `main` (commit series ending `f95872a`):

| Tool | What it settles | Evidence |
|---|---|---|
| `analysis/g1k_allele_analyzer.py` (merged with PR #1) | class frequency/age/backgrounds per locus; TSV ingestion; `--anchors K` | LCT derived age 516.5 vs 495 documented; ancestral 3 backgrounds |
| `analysis/sv_phase_pilot.py` | **exact** single-origin verdict over every admissible phase; phase-free enrichment | 715-carrier clade with 615 unknown phase bits → 1 origin at true age 1039 gen, p=0.001; scattered control → NO-CLADE, p=0.93 |
| `analysis/sv_anchor_match.py` | anchor/stitch an assembly haplotype onto the panel; recombination-aware segmentation | recombinant anchor recovered as 2 tracts at 0/97 and 0/113 mismatches, breakpoint within one site spacing |
| `analysis/g1k_loci_build.py` (PR #1) | loci-only ARG rebuild, 167 chunks | 167/167 chunks; chr20 length fix |
| `analysis/temporal_traj.py` (PR #1) | ancient class-frequency trajectories with drift + date-permutation nulls | 13 robust rising classes, 3–7 of 7 countries |

Measured facts from the merged analysis phase that constrain the plan:

- **τ × cohort interaction:** global cohort, `dom_share` 0.37→0.49 and dominant background
  491→1,012 haplotypes as τ goes 200→500. Background statistics are only informative when
  τ is scaled to where the cohort's structure lives.
- **~99% of DBSCAN classes are not genealogical clades** (`age_pct` > 90 for nearly all):
  a sequence-similarity class over a megabase locus is recombination-scrambled relative to
  any single local tree.
- Local tree depth in this ARG reaches **54,327 generations** (chr22:15 Mb), i.e. the
  cohort's structure is far deeper than any single background cutoff.

## 2. What is now established as an open problem

1. **Windowing makes ages mutually inconsistent.** Everything dated is a 10 Mb window
   (294 chunks in `data/derived/g1k`, 312 in `argphased`, 313 in `arggenome` — three
   near-duplicate families). tsdate dates each window independently, so the same node gets
   two different ages in two adjacent windows. Every cross-locus age claim inherits this.
2. **A per-chromosome dating attempt already failed for a fixable reason.** `argmerge/chr2.trees`
   is a whole chromosome (1,392 samples, 59,889 sites, 12.3M edges, 244 Mb) but dating it
   died in tsdate 0.2.7: `AssertionError: Use fewer rescaling intervals`
   (`argmerge/chr2_date.log`). Not memory — this machine has 1 TB / 256 cores.
3. **Unphased carriers inflate scatter.** A sample-level carrier set puts the allele on the
   wrong haplotype half the time, which biases *against* the statistic we use to reject
   artifacts, and makes the age too old.
4. **τ is a knob, not a property.** Tuning it until `dom_share` looks convincing optimises
   the reported statistic. The usable window is bounded below by cohort structure and above
   by the carrier clade age (breadth → 1 by construction; already recorded in `a587215`).
5. **Anchor matching cannot phase non-donors.** Measured: 28/6404 haplotypes matched at
   ≥0.99 with a donor anchor, but only 148/200 (74%) of heterozygous non-donor samples got
   the correct haplotype. That 74% is the boundary of the method, not a bug.

---

## 3. Work packages

Each has an aim, a method, and a validation criterion that must pass before the package is
called done. Effort is person-days for someone who knows the code.

### WP0 — Pin the ARG (gate G0)

**Aim.** One canonical ARG per analysis, with its provenance recorded. Three windowed
families in one namespace will be mixed by accident, and PR #1's 167-chunk build is a fourth.

**Method.** Add `data/derived/README.md` naming the canonical directory, the chunk grid, and
the build provenance (tsinfer/tsdate versions, sample set, flanking). Add the build hash to
every result CSV. Free QC: date the same locus in `g1k` and `argphased` and report the
build-to-build age difference — we tested determinism (bit-identical re-runs) but never
cross-build agreement, and that number is the error bar on every age we have quoted.

**Validation.** A documented canonical path; cross-build age deltas reported for ≥5 shared
loci; every result file carries the ARG identity. **Effort: 0.5 day.**

### WP1 — Per-chromosome dating (gate G1)

**Aim.** Make ages comparable across the genome. A chromosome is a *closed coalescent unit* —
no recombination between chromosomes — so dating `chr<N>` as one tree sequence is not an
approximation, it is the exact restriction of a genome-wide ARG. The 10 Mb grid is the only
arbitrary cut in the current pipeline.

**Method.** Fix the rescaling failure first (fewer rescaling intervals, or a current tsdate
with `inside_outside`); date `argmerge/chr2.trees` end to end, then build and date the
remaining 21 chromosomes as background jobs. Re-derive the Ne prior genome-wide and re-check
the LCT anchor, because a change in dating scale shifts **every** age we have quoted,
including the AMY 8–17 gen and LCT 495 gen calibrations.

**Validation.** LCT derived age within ~5% of 495; ancestral allele still resolves to 3
backgrounds; per-chromosome vs windowed ages compared at ≥500 shared nodes with the
difference distribution reported. If the difference is small, windows stay usable for set-based
scans and ages get a documented correction; if it is large, per-chromosome becomes canonical
for anything age-based. **Effort: 2–4 days, mostly unattended compute.**

### WP2 — Modern phasing, route A: anchor and mosaic (part of gate G2)

**Aim.** Attach each SV call to the specific panel haplotypes that carry it, using HPRC v2
named haplotype subsequences.

**Method.** Align the subsequence to the reference to get a local VCF → compare **base
strings** against the ARG (never allele indexes: the ARG is ancestral/derived, the VCF is
REF/ALT) → windowed candidate selection → binary segmentation with a minimum tract length →
per-segment haplotype set and mismatch count. Exclude the SV interior (±2 kb) where the
alternate path is structurally divergent. Widen flanks to hundreds of kb: only ~30% of sites
are informative (498 of 708 carried the ancestral base in the chr22 test), and that sets the
resolution.

**Validation.** **Hard gate:** for every anchor haplotype that is itself a panel member,
require rank-1 recovery with 0 mismatches. Hold out a subset of HPRC haplotypes, anchor from
flanks only, and require the held-out haplotype to be recovered at rank 1 above a stated
flank width. Report the flank width at which rank-1 recovery reaches 95%. **Effort: 3 days.**

### WP3 — Modern phasing, route B: exact clade consistency at scale (part of gate G2)

**Aim.** For all 265 loci, produce the τ-free verdicts: `single_origin`, `n_origins`,
`clade_fit_errors`, and carrier-clade age with its phase-induced range.

**Method.** Run `sv_phase_pilot.py` over the merged input tables (moderns, all 3,202). No
sampling is needed — the consistency test is polynomial — so cost is dominated by I/O.

**Validation.** Selftest must pass; verdicts must be stable across the 3 τ values and across
the 9 anchors (a verdict that flips with the anchor is a mosaic, not a clade, and must be
labelled as such). **Effort: 2 days.**

### WP4 — Reconcile A against B, and adjudicate the class list (gate G2 closes here)

**Aim.** Two independent routes to the same quantity. Agreement is the paper's confidence
statement; disagreement localises a bad call.

**Method.** Build the 2×2: anchor-derived background versus clade-consistency verdict, per
SV. Locus classes: (i) both agree single origin → report with a date; (ii) anchor says one
background but carriers are not a clade → mosaic or merged call, segment it; (iii) clade but
anchor disagrees → suspect the anchor (novel variants, sparse flanks); (iv) neither → artifact
candidate, check segdups and mappability.

Then apply the verdicts to the existing candidate lists: the 39 screen loci
(`docs/sv-screen-passes.md`) and the 13 robust rising classes (PR #1). The question is narrow
and answerable: **is each candidate a clade, a mosaic, or a mixture?** A class that is a
mixture cannot carry a convincing rise window, because the label is not a lineage.

**Validation.** Disagreement rate reported; every one of 13 + 39 adjudicated with a stated
class; the counts must reproduce from the committed scripts. **Effort: 3 days.**

### WP5 — Replace τ with a profile and a crossover time

**Aim.** Stop reporting breadth at a hand-tuned cutoff.

**Method.** For each locus compute `dom_share(τ)`, `n_backgrounds(τ)` over the sweep already
run (200/500/1000, extended to 50…20000), and derive **τ\***, the crossover at which carriers
become confined (e.g. `dom_share` crosses 0.5). Report `dom_share(τ)` curves for the headline
set. Keep the local empirical null for significance; do not anchor τ to a central quantile of
the mixed-cohort TMRCA distribution (the median is dominated by cross-continental splits and
lands near tens of thousands of generations, at which point breadth is 1 by construction).
If one cohort-level scale is needed for ranking, use a stated low quantile of the
**within-superpopulation** pairwise-TMRCA distribution and report which.

**Validation.** τ\* reproducible; conclusions rank-consistent across τ ∈ {200, 500, 1000};
the τ\* difference between mixed and superpopulation cohorts reported as a result.
**Effort: 2 days.**

### WP6 — Ancient placement, and the ascertainment model (gate G3)

**Aim.** Ancient SV presence through time, per lineage. **Not phasing** — a 1× pseudo-haploid
call is one haplotype draw, so there is no phase to resolve; the ancient is *placed* on the
modern genealogy.

**Method.** (a) SV alleles per ancient from read evidence, with the genotype class recorded:
homozygous/hemizygous deletion and CNV class are reliable at 1×, **heterozygous deletion is
the residual uncertainty**. (b) Placement from variation **independent of the SV** (flanking
markers only — placing with the SV is circular). (c) Per-sample placement confidence, and the
analysis re-run with ambiguous placements dropped. (d) **Detection-probability model as a
function of coverage and date**, with negative-control loci (common neutral deletions matched
for span and repeat content) pushed through the identical pipeline.

**Validation.** The rise window must exceed what detection-rate drift alone produces on the
controls. An ancient carrying an SV yields a hard age bound; report the fraction of loci where
that bound is consistent with the WP1 tsdate age (this doubles as clock calibration).
**Effort: 4 days — the largest single item, and the one with the highest false-positive risk.**

### WP7 — Lineage-level trajectories

**Aim.** Turn PR #1's class-label trajectories into lineage trajectories, which is what
MISSION §17 asks for ("frequency by time bin within ARG-derived background").

**Method.** Replace the class-frequency-in-ancients layer with node/segment-level lineage
frequencies from WP6 placement, binned by date. Test the discriminating question directly:
does the SV rise *within* a background (in-situ selection) or does the background itself
replace others (turnover)? This is decision point 5 on a genealogical axis instead of a label
axis.

**Validation.** In-situ versus turnover verdicts on simulated trajectories with known
behaviour; the 13 robust classes must either survive on the lineage axis or be reclassified,
and the reclassification reasons stated. **Effort: 3 days.**

### WP8 — Homoplasy as a measurement

**Aim.** At the fast-mutating repeat/CNV loci, homoplasy is not a nuisance — it is the signal.

**Method.** Count independent origins per allele (`n_origins`) and compare with a
mutation-rate null (stepwise for copy-number loci, infinite-sites for biallelic deletions).
Decompose identity into identity-by-descent (same clade) versus identity-by-state (convergence
on distant branches) and derive a per-locus μ.

**Validation.** Recover a known stepwise μ from simulated copy-number loci on the real local
trees; report the realised precision. **Effort: 3 days.**

---

## 4. Sequencing and gates

| Gate | Claim it guards | Passes when | Blocks |
|---|---|---|---|
| **G0** ARG pinned | "we can reproduce this number" | WP0 | everything |
| **G1** dating consistent | any cross-locus age comparison | WP1 | WP4, WP6, WP8 |
| **G2** phasing valid | any single-origin claim | WP2 + WP3 + WP4 | WP4 headline, WP7 |
| **G3** ascertainment modelled | any ancient "rise" claim | WP6 | WP7 |

Recommended order, exploiting parallelism (this box has 256 cores, so WP1's compute runs
unattended alongside analysis):

```
WP0  (0.5d)
  ├─ WP1 per-chromosome dating ──────── unattended, 2–4d ──┐
  ├─ WP2 anchor/mosaic phasing (3d)                        │
  │    └─ WP3 exact clade verdicts (2d)                    │
  │         └─ WP4 reconciliation + class adjudication (3d)─┤ G2
  ├─ WP5 tau* profiles (2d)  ← can start immediately, needs no gate
  └─ WP6 ancient placement + ascertainment (4d) ────────────┤ G3
       └─ WP7 lineage trajectories (3d)                     │
WP8 homoplasy/mutation rate (3d) ← after G1
```

## 5. Phasing: what we will and will not do

**Will.** Anchor assembly haplotypes to the panel by base-string matching with recombination
segmentation; use the exact clade-consistency test as the referee, since it is polynomial and
needs no sampling; place (not phase) pseudo-haploid ancients; treat the donor/non-donor
asymmetry explicitly — donor anchors are observations, non-donor anchors are background
identifications with a measured 74% per-sample phase ceiling.

**Will not.** Use Beagle / SHAPEIT / MVNcall. A deletion leaves carriers hemizygous (one
copy); an inversion is an LD island with suppressed recombination; and the statistic we want —
confinement to a clade — is exactly what a Li–Stephens panel phaser optimises toward, so it
would assume the answer. Homoplasy makes it worse: such a phaser has no model for the same
allele arising on two branches, so it resolves genuine recurrence into a single-copy phase and
erases the signal.

## 6. Anti-goals

- No novel SV discovery in ancients (1× gives no breakpoint support).
- No age comparison across independently dated windows without a calibration curve.
- No τ tuned until `dom_share` looks convincing.
- No rise window reported without the detection-probability control.
- No headline claim about a class before WP4 has said whether it is a clade.

## 7. Decision-register entries to add (`docs/DECISIONS.md`)

1. Canonical ARG directory and build provenance (WP0).
2. Per-chromosome dating as the age standard, or windows-with-correction, with the measured
   difference that decided it (WP1).
3. Single-origin verdicts are τ-free; breadth is reported as a τ profile with τ\* (WP5).
4. Ancients are placed, never phased; heterozygous deletion calls are the reported uncertainty
   (WP6).
5. Class-level candidate lists are provisional until adjudicated clade/mosaic/mixture (WP4).

## 8. Reproduction

```bash
cd analysis
../.venv/bin/python sv_phase_pilot.py --selftest          # 5 s, prints PASS
../.venv/bin/python sv_anchor_match.py --recomb-selftest --flank 200000
../.venv/bin/python sv_phase_pilot.py --table loci.tsv --out pilot.csv
G1K_ARG_DIR=/path/to/chunks ../.venv/bin/python g1k_allele_analyzer.py --table t.tsv
```

Collaborator bundle (analyzer + both SV tools + ARG chunks): `/moosefs/erikg/demoselex-arg`.
