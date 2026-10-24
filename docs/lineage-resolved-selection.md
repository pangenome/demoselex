# Lineage-resolved selection: what becomes askable

Status: proposal / framing note. Not a task list (that is
`docs/next-steps-arg-sv-phasing.md`) and not a mission amendment. This note states the
objective and the space of questions it opens.

## 1. The objective, stated plainly

We are not trying to reproduce a candidate list. The existing candidate funnel
(200 RISING -> 58 -> 29 -> **13 ROBUST**) is a good screen built on the only unit
available to it — a haplotype *class label* — and its top rows are the ones where that
unit is most likely to be wrong.

The objective is to change the unit of analysis.

> Ask whether an allele rose **because it was selected**, by following it as an inherited
> lineage through time, rather than asking whether a label's frequency changed faster than
> a demographic model permits.

Both halves of that sentence are load-bearing. "Inherited lineage" replaces the label.
"Through time" means ancients are evidence about descent, not just a frequency series.
"Because it was selected" means the claim must survive a null built from the genealogy that
the allele actually lived in.

## 2. Why the current unit cannot answer the question

A frequency over time is anonymous. `rise > 0` with a drift z and a date-permutation p is
compatible with at least four different worlds:

1. the allele rose inside a persistent background — the selection-like case;
2. the background itself was replaced by another — ancestry turnover, which is structure,
   not selection;
3. the class definition drifted, so its frequency changed while nothing biological did;
4. calling quality changed over time, so detectability changed.

The current funnel removes (2) by restricting moderns to one region and removes (4) only
implicitly. It cannot remove (3) at all, because a class has no identity beyond its
definition. And the measured evidence says (3) is the common case: ~99% of classes have a
within-locus age percentile above 90, i.e. a megabase similarity class is not a genealogical
clade. A rise measured on a non-lineage is a rise of a label.

The strongest available statement about a *label* is therefore descriptive. The strongest
available statement about a *lineage* is mechanistic. That is the whole argument.

## 3. The three substitutions

**Unit: allele -> allele on a lineage.**
An allele is placed on the specific haplotype backgrounds that carry it, at the leaves of
the dated ARG, with a carrier clade where one exists. Phase is not invented: for moderns,
heterozygote assignment is exact given a clade (the consistency test is polynomial); where
an HPRC v2 assembly haplotype is available the placement is an observation; for ancients,
a pseudo-haploid call *is* one haplotype draw, so there is nothing to phase — the sample is
placed.

**Null: demographic model -> realized genealogy.**
Today the null is `sqrt(p(1-p)t/2Ne)` at fixed Ne, genome-wide. After placement the null is
the genealogy the allele actually sat in: its clade's age, its clade's branching structure,
its background's size, and size-matched clades from the same local tree. Drift stops being
an assumption and becomes an observation.

**Time: bin index -> dated interval.**
An allele on a clade has an age interval (carrier-clade MRCA below, that node's parent
above) and, with ancients, an independent bound from presence/absence. Two unrelated data
sources constrain one number.

## 4. What the data contribute, exactly

| data | contribution | what it cannot do |
|---|---|---|
| 1000G/HGDP phased haplotypes (the ARG leaves) | the genealogy; dates; the background partition | nothing about the past directly |
| HPRC v2 haplotype subsequences | anchor a call to a named haplotype: an observation, not an inference, for donors | phase non-donors (measured ceiling ~74% per-sample) |
| SV/CNV calls at the 265 loci | a variant class that SNP scans are blind to, in the SNP-poorest regions | rare SVs (ascertainment: common only) |
| SV calls in 1,627 dated ancients | presence/absence through time; hard lower bounds on age | novel SV discovery; het deletion calls are the weak class at 1x |
| loci-only + full-genome ARG | two independent windowed builds agreeing to ~1% on age | window-boundary date consistency (unresolved) |

## 5. The questions this makes askable

Each entry: the question, the statistic, the null, and the artifact it kills.

**Q1. Did the allele rise inside its background, or did the background replace others?**
Decompose allele frequency change per time bin into within-lineage and between-lineage parts,
using placed ancients. Null: the same decomposition under date permutation. Kills turnover.

**Q2. Is the carrier set even a lineage?**
Exact clade-consistency test over every admissible phase. Returns single-origin, the clade,
the induced phase, and `clade_fit_errors` (the minimum number of genotypes that must be wrong
to make the carriers cladal). Kills label drift, merged calls, mismapping. **This is the
cheapest question and it should be asked of every candidate first.**

**Q3. How old is it, and does the ancient record agree?**
Age interval from the edge bracket plus a binomial absence bound per time bin
(zero carriers in N samples => frequency below ~3/N at that date). Disagreement is
informative: an old clade with a late frequency bound is standing variation, not a sweep.
Also calibrates the dating itself, independently of the LCT anchor.

**Q4. How many origins, and on which backgrounds?**
`n_origins` and the per-origin clades. Multiple origins of the same state on distinct
backgrounds is the strongest evidence class available, because the replication is across
lineages and no population label is involved.

**Q5. What mode of selection?**
Hard sweep (one young clade, star-like branching) versus soft sweep (carriers across several
older backgrounds, each young-ish) versus standing variation (old background, young allele on
it) versus recurrent (multiple origins). Statistic: clade branching times against size-matched
neutral clades in the same tree; and age(allele) versus age(background).

**Q6. Which variant in the interval is causal?**
Order markers by age along the locus — the placed SV plus flanking SNPs — and locate the age
minimum. The causal allele should sit at or near it. Statistic: local age profile; null: the
distribution of age minima in matched neutral intervals.

**Q7. Is the signal replicated on independent lineages?**
The mission's own question, asked on an axis the structure correction does not touch. Report
per-lineage rise and its date; replication = same-direction rise on disjoint lineages.
Compare directly against the label-based and country-based versions of the same statistic.

**Q8. Is it hypermutability rather than selection?**
Per-locus mutation rate from independent origins against a stepwise (CNV) or infinite-sites
(deletion) null. Mandatory at repeat-rich loci: a high-frequency copy-number state is the
expected outcome of a hot locus, not of selection.

**Q9. Is the call set itself good?**
Clade-consistency per call, aggregated per locus, profiles mappability deserts and reference
bias. Retro-fits every other number downward or upward; should run first and last.

## 6. The output that does not exist today

One row per (allele, locus, lineage), with time on it:

```
locus  chrom  start  end  allele  n_carriers  freq_modern
single_origin  n_origins  n_feasible_clades  clade_fit_errors
age_lo  age_hi  age_median  phase_induced
background_n  dom_share  dom_share_tau200  dom_share_tau500  tau_star
first_ancient_bp  ancient_n_bin  ancient_freq_bound  detection_adj_p
mu_locus  origins_per_unit_time
within_rise  between_rise  mode  calls  provenance_hash
```

This is the mission's §17 "frequency by time bin within ARG-derived background" made
concrete, with the confound columns that let a reader check each claim rather than trust it.

## 7. What must be true for any of it to hold

- **G2 first.** Placement validity is the gate for everything in §5. No single-origin claim
  before the anchor and the exact test agree, and before the HPRC self-recovery check passes.
- **G3 before any ancient rise.** Detection probability is a function of sample age; without
  matched negative controls, every rise is confounded with improving aDNA quality.
- **G1 for cross-locus age comparisons.** The two windowed builds agree to ~1%, which is
  encouraging but is not the same as window-boundary consistency; per-chromosome dating
  settles it.
- **One canonical ARG**, hash recorded in every output.

## 8. What we will not claim

- Not that any of the 13 risers is selected. They are re-derived or retired, not defended.
- Not that structural variants are generally under selection. The expected outcome is that
  most SVs are old and drift-maintained.
- Not selection coefficients as point estimates. Bounds, and only under stated assumptions.
- Not absence of an allele in ancients as proof of absence: only a frequency bound, with the
  bin's N and sampled region attached.
- Not cross-window age comparison without a calibration curve.

## 9. Why this is a different question from the ones already answered

Akbari et al. ask whether an allele-frequency change survives a structure correction built
from population labels. The existing funnel asks whether a class-frequency change exceeds a
demographic drift model, replicated across countries. Both are label-based, and both must
therefore treat ancestry as nuisance to be subtracted.

This asks whether the allele rose *inside a lineage*, which turns ancestry from a nuisance
into the measurement. That is not a stronger version of the same test; it is the same
question asked in the coordinate system where the confound does not exist.

## 10. What success looks like

Success is not more candidates. It is:

1. a candidate list where every row has a lineage, a date interval, a mode, an origin count,
   a mutation-rate null, and two independent age constraints;
2. the 13 re-derived into whatever survives — plausibly a few lineages and a documented
   explanation for the rest;
3. the eight "distributed" classes resolved as either segmented lineages or explicitly
   retired as labels without identity;
4. a clean negative for the class as a whole: structural variants are largely old and
   drift-maintained, with a small number of exceptions that are mechanistically specified.

A defensible negative plus a handful of mechanisms is the outcome that survives review. That
is the bar, and it is a higher bar than a longer list.
