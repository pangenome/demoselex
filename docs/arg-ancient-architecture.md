# The architecture: what we infer, what we place, what each can establish

Status: methods statement. Companion to `docs/lineage-resolved-selection.md` (the objective)
and `docs/next-steps-arg-sv-phasing.md` (the work plan). This note fixes the design and, more
importantly, the *claim hierarchy* — how much weight each kind of evidence is allowed to
carry.

## 1. The architecture in one sentence

Infer the genealogy from **modern phased haplotypes only**; **place** ancient samples onto it;
use the ancients as evidence about **time**, never as members of the tree.

## 2. Why the moderns must be the substrate

Forced by data quality, not preference. tsinfer and tsdate need dense, phased, low-error
genotypes across thousands of samples. The ancient set is ~1x, pseudo-haploid, post-mortem
damaged, and irregularly distributed in time and space. Including them in inference would
either be impossible or would degrade the tree everywhere to accommodate the noisiest samples.
The modern-only ARG is the only version buildable at this scale.

## 3. Why the ancients stay out — and why that is the asset

If ancients were in the ARG, age estimates and temporal evidence would share their errors, and
the test would be circular: validating the inference with itself. Keeping them out makes the
ancient record a genuine **out-of-sample** check. This is the reason "cross-validation" is the
right word, and it is worth protecting deliberately — e.g. never tune the ARG or its dating on
the same loci that will be used as the ancient test.

## 4. It is a calibration loop, not a one-way check

The relationship runs both ways and this is the part usually left implicit:

- **Ancients validate the ARG's dates.** A clade dated to 500 generations should be present in
  ~500-generation-old samples and absent in samples several times older.
- **Ancients calibrate the absolute clock.** ARG times are relative, scaled by a mutation-rate
  prior and demographic assumptions. Radiocarbon-dated carriers are the only absolute anchors
  available. LCT is one textbook anchor; 1,627 dated ancients are a distributed one.

So the sequence is infer -> constrain -> revise, not infer -> check. This is stronger than
validation, and it is what makes a published age defensible rather than merely internally
consistent.

## 5. What each side can and cannot establish

| claim | establishment | source |
|---|---|---|
| a clade's **age** | direct test: presence/absence versus sample date | ancients |
| an allele's **frequency over time** | binomial bound per time bin (0 carriers in N => p below ~3/N) | ancients |
| the **absolute clock** | calibration against dated carriers | ancients |
| that a **background is real** | weak: needs confident placement | ancients + placement |
| the ARG's **topology** | **not testable this way** — ancients are not in the tree | neither |

The fourth and fifth rows are where overclaiming happens. "We cross-validated the ARG" is
false as stated. The true statements are: the **dates** are cross-validated, the
**frequencies** are bounded, the topology is inferred and unverified by this route.

## 6. The claim hierarchy (mandatory)

Every ancient-derived statement carries the strength of the evidence behind it. In order:

1. **Time.** Presence of an allele in a dated sample proves it existed by that date; absence in
   N samples bounds its frequency at that date. Robust to placement.
2. **Calibration.** The same presence/absence structure tests and scales the ARG's clock.
   Robust to placement.
3. **Ancestry identity.** "This ancient carried the allele on background B" depends entirely
   on placement, and placement is a nearest-relative match. Weakest of the three; must carry a
   placement confidence, and conclusions must be re-run with low-confidence placements dropped.

The ordering is not stylistic. Rows 1 and 2 survive placement error; row 3 does not. Reporting
them at the same strength is the most likely way this project produces a wrong paper.

## 7. The limitation to state up front

**Placement finds the nearest modern relative, not the actual ancestor.** Ancient gene pools
contained lineages that left no modern descendants. Those lineages have no branch in a
modern-only ARG, so an ancient can be placed on background B merely because B is the closest
*surviving* relative while its true lineage is extinct.

Consequences, all of which follow the hierarchy in section 6:

- presence/absence bounds are unaffected — they do not depend on where the ancient sits;
- clock calibration is unaffected, for the same reason;
- **"same background" claims are conditional** and should be phrased as "nearest modern
  background", never as descent.

There is a second, smaller limitation in the same direction: ancients can only be placed at
variation present in the modern panel, so variation private to ancient samples is invisible.

## 8. Alternatives considered

**Serial-sample ARG (ancients inside the tree).** Theoretically cleaner: one tree sequence whose
samples come from different times. It breaks the coalescent assumptions that tsinfer and tsdate
rely on, and is tractable today only for a modest number of high-coverage genomes (~20-30x), not
1,627 samples at ~1x. This is a research programme, not a step in this plan.

**No ARG; IBD / haplotype matching only.** Simpler, and it loses exactly what makes the
questions answerable: the dates.

**Current design** (moderns inferred, ancients placed) is the only one that both builds at this
scale and retains absolute time. It is the standard architecture for exactly this reason.

## 9. What this forces on every reported result

1. Each claim names its evidence class (section 6) and its data source.
2. Ancient ancestry statements carry a placement confidence, with a sensitivity analysis that
   drops low-confidence placements.
3. Any "genome-wide" statement carries the ARG's actual coverage, not its nominal grid (four
   chunks are truncated mid-chromosome: chr13/chr14 chr14-20 Mb, chr15 chr13-20 Mb, chr22
   chr19.5-20 Mb).
4. Presence and absence are never reported as symmetric evidence: absence is a frequency bound
   with the bin's N and sampled region attached.
5. The extinct-lineage caveat appears wherever a background assignment drives a conclusion.

## 10. Summary

Modern ARG for inference, ancients for time. The value is not that the ancients check the tree
— they cannot — but that they turn relative times into absolute ones and bound frequencies over
time, independently of the inference. Where the plan asks an ancient to say *which* background
it was on, it is asking the one question this architecture answers least well, and the plan
should keep asking it only with the confidence attached.
