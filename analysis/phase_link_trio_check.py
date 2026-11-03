#!/usr/bin/env python3
"""Phase-aware Mendelian validation of the structural phase link.

For every (trio, locus) where the child and at least one parent have decided
orientations:
1. identify the child's PATERNAL/MATERNAL ARG haplotypes at the locus by
   SNP matching (child hap vs father's and mother's phased ARG haps at sites
   where the parents differ);
2. do the same for the transmitted parent haplotype (parent hap closest to
   the child's inherited hap);
3. CHECK: the structural cluster riding the child's inherited hap must equal
   the parent's cluster on the transmitted hap (same cluster label), or, when
   the parent's own phase is decided, the labels must match on the matched
   haps. Any mismatch is a phase or call error.

Inputs: phase_link_genome.csv (orientations), geno master, ARG chunks,
        samples_population.txt (FamilyID, SampleID, FatherID, MotherID).
Output: per-check rows + summary: agreement rate of phase-aware transmission,
        stratified by self-assigned vs catalog rows and by locus.
"""
import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import tskit

PHASE_DIR = Path('/moosefs/laurap/old_dups/demoselex_exp/phase_link_run')
ARG_DIR = Path('/moosefs/erikg/demoselex-arg/arg')
POPFILE = Path('/moosefs/erikg/demoselex-arg/samples_population.txt')
MASTER = Path('/moosefs/raid5/lpignata/old_dups/geno_3202_fa/geno_3202_master.tsv.gz')
CHUNK_RE = re.compile(r'chr(\d+)_(\d+)_(\d+)\.dated\.trees')


def load_arg_sites(chrom, lo, hi):
    sites = []
    for f in sorted(ARG_DIR.glob(f'chr{chrom}_*.dated.trees')):
        m = CHUNK_RE.match(f.name)
        if m and int(m.group(2)) <= hi and lo < int(m.group(3)):
            ts = tskit.load(str(f))
            for v in ts.variants():
                p = int(m.group(2)) + int(v.position)
                if lo <= p < hi:
                    sites.append((p, v.alleles, v.genotypes))
    sites.sort(key=lambda x: x[0])
    if not sites:
        return None
    return (np.array([p for p, a, g in sites]),
            {p: a for p, a, g in sites},
            np.array([g for p, a, g in sites]).T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='phase_link_trio_check.csv')
    a = ap.parse_args()

    pop = {}
    for line in open(POPFILE):
        f = line.split()
        if f[0] == 'FamilyID' or len(f) < 7:
            continue
        pop[f[1]] = (f[2], f[3])   # father, mother
    in_cohort = set(pop)

    # trios: child with both parents in cohort
    trios = {c: (fa, mo) for c, (fa, mo) in pop.items()
             if fa != '0' and mo != '0' and fa in in_cohort and mo in in_cohort}
    print(f'trios with both parents in cohort: {len(trios)}', file=sys.stderr)

    # phase-link rows keyed (sample, locus)
    pl = {}
    for f in PHASE_DIR.glob('*.csv'):
        if f.name.startswith('phase_link'):
            continue
        for r in csv.DictReader(open(f)):
            if r['orient'] in ('S', 'SW'):
                pl[(r['sample'], r['locus'])] = r
    print(f'decided phase-link rows: {len(pl)}', file=sys.stderr)

    # locus spans
    bed = {}
    for line in open('/moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed'):
        f = line.split()
        if len(f) >= 4:
            c = f[0].split('#')[-1].replace('chr', '')
            bed[f[3]] = (int(c), int(f[1]), int(f[2]))

    # group needed (trio, locus) by locus
    need = defaultdict(list)
    for child, (fa, mo) in trios.items():
        for locus in bed:
            if (child, locus) in pl and ((fa, locus) in pl or (mo, locus) in pl):
                need[locus].append((child, fa, mo))
    print(f'loci with usable trios: {len(need)}', file=sys.stderr)

    order = [l.split()[1] for l in open(POPFILE)][1:]
    idx = {s: i for i, s in enumerate(order)}

    out_rows = []
    for locus, trios_here in sorted(need.items()):
        chrom, lo, hi = bed[locus]
        arg = load_arg_sites(chrom, lo, hi)
        if arg is None:
            continue
        P, A, G = arg
        for child, fa, mo in trios_here:
            kc, kf, km = idx[child], idx[fa], idx[mo]
            chA, chB = G[2*kc], G[2*kc+1]
            # parent-hap closeness: for each child hap, count matches vs each parent's haps
            # use sites where father's two haps differ OR mother's differ
            inf = [j for j in range(len(P))
                   if G[2*kf][j] != G[2*kf+1][j] or G[2*km][j] != G[2*km+1][j]]
            if len(inf) < 15:
                continue
            def close(hap, kh):
                # returns (n_matching_sites) of hap vs parent kh's two haps
                a = sum(1 for j in inf if hap[j] == G[2*kh][j])
                b = sum(1 for j in inf if hap[j] == G[2*kh+1][j])
                return a, b
            chA_f, chA_m = close(chA, kf), close(chA, km)
            chB_f, chB_m = close(chB, kf), close(chB, km)
            # child hap A paternal if it is closer to father than to mother AND
            # closer to father than child hap B is
            A_pat = max(chA_f) > max(chA_m) and max(chA_f) >= max(chB_f)
            B_pat = max(chB_f) > max(chB_m) and max(chB_f) > max(chA_f)
            if A_pat == B_pat:      # undetermined parent-of-origin
                continue
            pat_hap = 0 if A_pat else 1          # child ARG hap index from father
            mat_hap = 1 - pat_hap
            # transmitted father hap: father's hap closest to child's paternal hap
            fA = sum(1 for j in inf if chA[j] == G[2*kf][j]) if pat_hap == 0 else sum(1 for j in inf if chB[j] == G[2*kf][j])
            fB = sum(1 for j in inf if chA[j] == G[2*kf+1][j]) if pat_hap == 0 else sum(1 for j in inf if chB[j] == G[2*kf+1][j])
            f_tx = 0 if fA >= fB else 1
            mA = sum(1 for j in inf if (chA if mat_hap == 0 else chB)[j] == G[2*km][j])
            mB = sum(1 for j in inf if (chA if mat_hap == 0 else chB)[j] == G[2*km+1][j])
            m_tx = 0 if mA >= mB else 1

            rc = pl[(child, locus)]
            ch_cls = [rc['arg_cluster_hap0'], rc['arg_cluster_hap1']]
            for parent, p_id, tx_hap in ((fa, kf, f_tx), (mo, km, m_tx)):
                rp = pl.get((parent, locus))
                if rp is None:
                    continue
                p_cls = [rp['arg_cluster_hap0'], rp['arg_cluster_hap1']]
                # the transmitted child hap for this parent
                child_hap = pat_hap if parent == fa else mat_hap
                child_cluster = ch_cls[child_hap]
                parent_cluster = p_cls[tx_hap]
                if child_cluster is None or parent_cluster is None:
                    continue
                # de novo call error: child's cluster absent from BOTH parent haps
                absent = child_cluster not in p_cls
                out_rows.append(dict(
                    locus=locus, child=child, parent=parent, which='F' if parent == fa else 'M',
                    child_hap=child_hap, tx_hap=tx_hap,
                    child_cluster=child_cluster, parent_cluster=parent_cluster,
                    transmitted_match=int(child_cluster == parent_cluster),
                    absent_from_parent=int(absent), n_sites=len(inf)))
    with open(a.out, 'w', newline='') as f:
        if out_rows:
            w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
            w.writeheader()
            w.writerows(out_rows)
    n = len(out_rows)
    ok = sum(r['transmitted_match'] for r in out_rows)
    absent = sum(r['absent_from_parent'] for r in out_rows)
    print(f'checks: {n} | phase-consistent transmissions: {ok} ({100*ok/n if n else 0:.1f}%) '
          f'| child cluster absent from parent entirely (de novo/call error): {absent} '
          f'({100*absent/n if n else 0:.2f}%)')


if __name__ == '__main__':
    main()
