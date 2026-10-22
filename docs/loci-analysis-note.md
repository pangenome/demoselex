# What we do with the loci — operational note

The 256 segmental-duplication / CNV loci, and how each one is processed into a statement
about selection versus demography. This is the "how we run it" note; the reasoning is in
`docs/methods-note.md`, the code map in `docs/pipeline-entrypoint.md`.

## 1. What a locus is, in data terms

Not a copy-number dosage. At each locus every individual carries **two alleles**, each
labelled with an integer (e.g. 1..10, or 1..20). The diplotype (a, b) *is* the copy-number /
haplotype configuration. One row per (locus, individual):

```
chrom,pos,sample,allele1,allele2
2,135851076,HG00096,1,2
2,135851076,HG00097,2,2
...
```

The allele labels come from our own genotyping, so they carry meaning independent of the
ARG (this matters: the ARG's internal allele orientation is arbitrary — see §7).

## 2. What we compute at each locus

Using the dated ARG at that locus (GRCh38), for **each distinct allele label**:

| quantity | meaning |
|---|---|
| `freq` | allele count / called haplotypes |
| `age` | tsdate time of the MRCA of that allele's carriers |
| `age_pct` | percentile of `age` among clades with the same descendant count in the same local tree (empirical, within-locus null) |
| `n_backgrounds` | how many ARG-derived local backgrounds hold ≥5% of that allele's carriers |
| `dom_share` | fraction of the allele's carriers that sit in its dominant background |
| `dom_freq` | the allele's frequency inside that dominant background |
| `background_n` | number of haplotypes in the dominant background |

A **local background** is the most recent clade, in the local tree at that position,
containing a haplotype whose root is younger than `tau` (default 200 generations). It is
derived from the genealogy alone — no population labels are used anywhere.

## 3. How to run it

```
.venv/bin/python analysis/g1k_allele_analyzer.py \
    --table diplo.csv --tau 200 --out results.csv
```

- `--table`: the diplotype CSV above; may contain many loci (`chrom`,`pos` group them).
- `--tau`: the background depth in generations. 200 is the calibrated default; breadth is
  the tau-stable quantity, so a sensitivity run at 100 and 500 is cheap and worth reporting.
- The analyzer restricts the ARG to the samples present in the table, so the genotyped and
  ARG cohorts must overlap. Samples absent from the ARG are skipped (and counted).
- `--selftest` recodes a known SNP as a two-allele diplotype; use it to confirm the install
  and to see the expected contrast.

## 4. Reading the output

Per locus, the interpretable pattern is in the **breadth** (how many backgrounds carry the
allele) crossed with **age** and **within-background frequency**:

| breadth | clade age | within-background freq | reading |
|---|---|---|---|
| 1 background | young | high (fixed in it) | candidate recent selection on that haplotype |
| 1 background | old | high | old haplotype, ancestrally confined (drift, or ancient selection) |
| many backgrounds | old | — | structure / carriage: the allele's distribution tracks ancestry turnover |
| many backgrounds | young | — | recent expansion across ancestries; check artefacts before interpreting |

Worked anchor (self-test, LCT rs4988235 recoded in EUR, tau = 200):

| allele | freq | age | age_pct | n_backgrounds | dom_freq |
|---|---|---|---|---|---|
| derived (ALT) | 0.536 | 495 gen | 20.0 | **1** | **1.00** |
| ancestral (REF) | 0.465 | 6,510 gen | 93.3 | 3 | 1.00 |

So at a locus like LCT, the derived allele is young, confined to one background and fixed
within it, while the ancestral allele is old and distributed — exactly the contrast the
method is designed to make visible.

## 5. QC per locus, before believing any row

1. **Enough called genotypes**: require a minimum number of called haplotypes; report the
   rate. Sparse loci give unstable backgrounds.
2. **Background count is meaningful**: if the locus has very few large backgrounds (say <3
   with ≥20 haplotypes), "breadth" is not informative at that locus — report it as such
   rather than as "confined".
3. **tau sensitivity**: report breadth at 100 / 200 / 500 generations; a reading that flips
   with tau is not a claim.
4. **Mapping/paralogy**: the allele labels are genotyped, so this is less of an issue than
   for SNP-based analysis, but still record each locus's mappability / duplication structure.
5. **Cohort overlap**: state how many genotyped samples were present in the ARG.

## 6. What we present

Per locus: allele label, frequency, clade age (with percentile), breadth, dominant-background
share and within-background frequency, and the number of backgrounds. Plus, in aggregate:

- the distribution of breadth across the 256 loci (how many CNV alleles are haplotype-confined
  versus ancestry-distributed),
- the joint distribution of age percentile and breadth (the two axes together separate
  regimes),
- the control rows from §4 as calibration,
- the skyline Ne(t) as the demographic frame for any absolute age quoted.

## 7. What we do not claim

- Not selection coefficients, and not dates beyond the resolution of the ARG (ages are
  substrate-dependent; compare percentiles, not absolute ages, across loci).
- Not "derived" versus "ancestral" from the ARG's own allele states: the build passed
  placeholder alleles to `tsinfer`, so its internal orientation is arbitrary. The diplotype
  labels sidestep this, but any statement about which allele is derived must come from the
  genotyping/annotation, not from the tree.
- Not a p-value from `age_pct`: it is an empirical percentile within the local tree, not a
  test statistic with a sampling distribution. Breadth and within-background frequency are
  descriptive; treat them as such.

---

## BED-level scan (available now, before the diplotype data)

`analysis/g1k_bed_scan.py loci.bed --pop EUR --breadth` summarises any BED whose 4th column
holds a gene/locus name. Per locus it reports: number of scored variants, how many fall in
the young-clade bottom 5% (k5) versus expectation, the Poisson p, the most extreme variant's
percentile and age, whether the locus is in the segdup compartment, and (with `--breadth`)
the local-background breadth of that extreme variant.

Demonstration on six gene-named loci (GRCh38, +/-20 kb, EUR):

| locus | n variants | k5/expected | enrichment | Poisson p | best variant pct | segdup | breadth |
|---|---|---|---|---|---|---|---|
| LCT | 162 | 17 / 8.1 | **2.10x** | **0.004** | 1.4 | no | 1 |
| AMY1A | 1 | 1 / 0.1 | 20x (n=1) | 0.049 | 10.9 | yes | 1 |
| HLA-B | 1,862 | 29 / 93.1 | 0.31x | 1.00 | 5.5 | no | 1 |
| CYP2D6 | 190 | 5 / 9.5 | 0.53x | 0.96 | 3.7 | yes | 1 |
| GSTM1 | 117 | 3 / 5.9 | 0.51x | 0.93 | 2.6 | yes | 1 |
| SLC24A5 | 16 | 0 / 0.8 | 0.00x | 1.00 | 73.2 | no | 3 |

Two things this shows, and they are the argument for structural genotyping:

1. The SNP-based local-ancestry metric fires where a SNP sweep exists (LCT, 2.1x, p=0.004)
   and correctly stays flat at the old-allele locus (SLC24A5).
2. It is **depleted** at exactly the loci whose biology is structural rather than SNP-level
   (CYP2D6 0.53x, GSTM1 0.51x). Those loci need the structural haplotypes - i.e. COSI-GT
   output - not SNP clade ages. An earlier version of this note added "and near-blind at AMY
   because the duplication region carries almost no common SNPs (1 variant in AMY1A +/-20 kb)".
   **That measurement was made around the gene, not the locus, and it is wrong at locus scale**
   - see the 2026-10-22 correction below.

Note the breadth column refers to the locus's *most extreme* variant, which at a locus like
SLC24A5 is not the canonical SNP; when the diplotype data arrive, breadth should be computed
per allele label with `g1k_allele_analyzer.py`, which is the authoritative path.

---

## CORRECTION (2026-10-20): breadth is only valid when tau < the carrier clade's age

A haplotype's background is its most recent clade *under* tau. If a variant's carriers share
an ancestor younger than tau, all their lineages merge below tau and therefore terminate at
the same node: breadth = 1 BY CONSTRUCTION. Demonstrated: six COSI-GT SV loci with carrier
clade ages 7-9 generations give breadth = 1 at tau = 200, 2,000 and 10,000 alike.

Consequence for interpretation:

- **Young alleles (clade age < tau)**: breadth is vacuous. The readout is the clade age and
  its percentile. Do not report breadth as evidence of "confinement".
- **Old alleles (clade age > tau)**: breadth is informative and is the carriage/selection
  discriminator - SLC24A5 (age ~8,270 gen, tau 200) spans 3-6 backgrounds; LCT (age 495) is
  confined to 1.

So the two statistics cover different regimes, and the decision table in section 4 must be
read with the allele's age in hand: age first, breadth only where age > tau.

Measured cohort interaction (2026-10-22): at tau 200-500 in the global cohort, backgrounds are
fragmentary - average dom_share 0.37 -> 0.49, dominant background 491 -> 1,012 haplotypes -
because a mixed cohort holds structure far deeper than any single cutoff. Breadth therefore
depends on a tunable constant, and tuning tau until dom_share looks convincing is optimising
the statistic being reported. Preferred practice: keep the verdict tau-free (clade age and
n_origins do not use tau at all), and report breadth as a PROFILE dom_share(tau) over the
sweep, summarised by the crossover time tau* at which carriers become confined. Do not anchor
tau to the median pairwise TMRCA of a mixed cohort: that median is dominated by
cross-continental splits and lands near tens of thousands of generations (local tree depth at
chr22:15 Mb is 54,327 generations), where breadth is 1 by construction.

---

## CORRECTION (2026-10-22): the LOCUS is the inference frame, not the gene

The statement retired above ("near-blind at AMY, 1 variant in +/-20 kb") was measured around
the gene. Measured on the BED intervals instead:

| locus | BED interval (GRCh38) | span | variants (MAC >= 5) | spacing |
|---|---|---|---|---|
| AMY1A | chr1:103,304,997-103,901,127 | 596 kb | **667** | 1 / 893 bp |
| GSTM1 | chr1:109,561,616-109,802,331 | 241 kb | **509** | 1 / 472 bp |

AMY1A is not SNP-poor. It is one of the strongest screen passes in the 265-locus scan
(1.69x, 629 variants scored, q = 3.0e-3). What is real is a HOLE in the middle of the locus,
which is the duplicated gene cluster itself:

```
20 kb bins  103.305-103.625 Mb   4-124 variants per bin    <- clean left flank, 249 kb
            103.625-103.745 Mb   3,0,1,0,0,0-9             <- the hole (~125 kb)
            103.745-103.885 Mb   6-48 variants per bin     <- clean right flank, 126 kb
```

Segdup corroboration: the contiguous paralogous block is chr1:103,554,220-103,774,845
(221 kb). (38 records overlap the span summing to 89.5%, which over-counts - use the
contiguous block, not the summed overlap.)

**Rules that follow, and that all locus inference must obey:**

1. **Frame = the BED locus interval.** Not the SV call, not the gene. Flank widths, masks and
   anchors are all derived from the frame per locus - never a fixed +/-20 kb or +/-200 kb
   constant, which lands inside the hole at exactly the loci we care about.
2. **Mask the paralogous core, then infer from the clean flanks.** Build the mask from data
   first (bin the frame; flag bins below a fraction of the locus median density) because that
   is usable at all 265 loci with no external input; use segdup as corroboration, not the
   primary criterion.
3. **The core is bracketed, which is stronger than flanked.** At AMY1A the CNV sits between a
   249 kb clean block and a 126 kb clean block. Anchor from both sides and ask whether the two
   sides agree on the ancestry background: agreement means the CNV arose once on one lineage;
   disagreement means the CNV is itself an ancestry boundary (mis-assembly, or a genuinely
   recombinant haplotype). This is a sharper test than a single flank permits.
4. **Publish a per-locus resolution budget**: frame span, masked-core span, usable informative
   variants (excluding sites where the allele is ancestral), spacing, and each flank's width.
   That number is the resolution available for anchoring and segmentation, and it decides per
   locus whether the marker set is sufficient.
5. **Trust is separate from count.** Variants that exist inside the core can still be
   paralogous mismapping. The exact clade-consistency test is a free mappability diagnostic:
   a mismapped marker shows inflated n_origins against the local tree. Report compatibility per
   bin.

Consequence for the tools: `sv_anchor_match.py` currently takes fixed flank windows around the
call, which is the wrong parameterisation for these loci; flanks should come from the frame
minus the mask. Where a locus's resolution budget is genuinely poor, the fallback is to put
the phased SV alleles into the ARG build itself (as biallelic sites in the tsinfer input),
restricted to SVs that pass the clade-consistency test so that only infinite-sites-compatible
markers are injected.

---

## COSI-GT SV locus scan (265 loci, HGSVC v3 refined set)

Source: `~/cosigt_paper/resources/benchmarking/svs/svs_refined.annotation.bed` (265 loci;
identical coordinates to `svs_refined.hgsvcv3.bed`, plus DEL/INS/COMPLEX classification and
top-allele type/length/sample counts; GRCh38). Script:
`analysis/g1k_bed_scan.py ... --pop EUR`; outputs
`data/derived/g1k/sv265_scan_annotated.csv`.

Result over 264 scorable loci:

- **As a class, SV loci are DEPLETED** in SNP-level young-clade signal: median enrichment
  0.78 (DEL), 0.78 (COMPLEX), 0.86 (INS) against the null expectation of 1.0. Segdup and
  non-segdup loci alike (0.83 vs 0.81). Structural loci are not, in general, sites of recent
  SNP sweeps.
- **But 39 loci pass FDR q < 0.05** (32 at q < 0.01), i.e. a strong tail of SV-adjacent
  regions with young-clade excess. Top by q: GNAI1 (4.12x, q=3e-26), NNT (3.64x, q=2e-19),
  MLXIPL (4.19x, q=1e-16), KATNAL2 (4.03x, q=2e-16), POLR2J2 (3.02x), HILPDA (2.59x),
  BMPR2 (3.06x), ADK (2.57x, best variant at the 0.02 percentile), CTAGE6, ARHGEF35, SCFD2,
  ACTA2, SESTD1, GEN1, **LPA (2.04x, q=2.6e-4)**, FBXO22, SPATA31D4, GOLGA6L4, RPL24, PRR20A.
- LPA is notable: KIV-2 copy number is one of the best-documented CNVs under selection, and
  its SV region shows young-clade excess in an analysis that never used copy number.

Caveats: the 39 hits are SNP-level signals *near* SV loci (not the SV alleles themselves);
some top loci are paralogous families (POLR2J2, CTAGE6, ARHGEF35, SPATA31D4, GOLGA6L4) where
paralogy contributes; and the window statistic's low-rho fragility applies. To test the SV
alleles themselves, the diplotype analyzer is the route.

---

## INPUT CONTRACT (agreed 2026-10-20)

Two files; nothing else required (groupings, population labels and CN values are not needed -
the alleles are just markers on the genome).

**1. Allele table** - one row per sample per locus:

```
chrom  pos  sample  allele1  allele2
1      103304997  HG00096  1  3
```

Integer labels consistent within a site; unordered diplotype (structural haplotypes are
unphased); missing = NA or -1. This is exactly what `g1k_allele_analyzer.py --table` reads.

**Unphased is not a blocker** (2026-10-22). The only ambiguity is one bit per heterozygote,
and the question "is this allele consistent with a single origin" is exactly decidable: a clade
C is consistent iff it contains both haplotypes of every homozygous carrier, no haplotype of a
homozygous non-carrier, and exactly one haplotype of every heterozygote. `sv_phase_pilot.py`
resolves this in polynomial time with no sampling - 715 carriers with 615 unknown phase bits
were recovered as one origin at the true age - and the induced phase comes out as a by-product.
So the allele table remains the input; phasing is not required from the caller. External
phasers (Beagle/SHAPEIT/MVNcall) are inappropriate here: a deletion leaves carriers
hemizygous, an inversion is an LD island, and such phasers have no model for the same allele
arising on two branches, which is exactly the recurrence we want to measure.

**2. Date table** - one row per sample:

```
sample  date_bp  date_sd  source  culture
```

Dates are NOT missing for AADR material: `data/raw/aadr-v66.p1/v66.p1_1240K.aadr.PUB.anno`
carries date_bp, date_sd and date_method (that is where the 1,000 dated ancients came from).
For samples outside AADR, use culture-based date INTERVALS with wide uncertainty - never
point dates - so precision is not manufactured. This is the one bounded literature/annotation
job standing between us and the temporal analysis.

Optional third file (already in hand from COSI-GT): `chrom start end gene classification`.

## ULTIMATE OUTPUT (agreed 2026-10-20)

One trajectory row per (locus, allele):

| field | meaning |
|---|---|
| locus, allele | which structural haplotype |
| freq_by_bin_by_background | frequency in each time slice, within each ARG-derived ancestry background |
| rise_window | when the frequency rose, and how narrow (punctuated vs gradual) |
| specificity | rise inside one background (in-situ selection) vs across backgrounds (turnover) |
| drift_null | expected frequency change over the interval from the local Ne(t) skyline, plus a date-permutation null |
| verdict | young + confined + within-background rise = candidate; old + broad = structure |

Rationale (why this and not a SNP scan): a single marker cannot tag these alleles - homoplasy
across duplicated copies destroys any single-SNP signal, which is why GWAS and SNP-based
selection scans fail at these loci. Allelic haplotype sharing, and its rate of increase
relative to the generations available, survives that and is the observable to use.

STATUS: allele analysis works today (`g1k_allele_analyzer.py`). The temporal layer
(`--dates`, time-binned within-background frequencies, drift nulls) is NOT yet written; it is
a ~150-line addition, testable by faking dates on moderns before real ancient data arrives.
