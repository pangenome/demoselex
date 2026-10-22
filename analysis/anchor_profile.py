#!/usr/bin/env python3
"""Per-anchor youth profiles for candidate structural haplotype classes.

A class is "young at an anchor" if the MRCA of its carriers at that anchor's
local tree is recent. A real young class (one birth event, linkage to the
flanking genealogy near that event, degraded by recombination with distance)
shows a CONTIGUOUS block of young anchors; a class that is young at one
scattered anchor is likely noise. The profile also localizes where along the
locus the youth signal lives.

Candidates are selected from an analyzer summary CSV (default: the
global_tau1000 run) by --max-age-min, --min-n and a frequency band, or given
explicitly as chrom,pos,allele lines via --pick.

Usage:
  anchor_profile.py [--summary global_tau1000/summary.csv] [--anchors 21]
                   [--young-age 2000] [--min-n 10] [--freq-lo 0.005] [--freq-hi 0.30]
                   [--max-age-min 1000] [--pick candidates.txt] [--out profiles.tsv]
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from g1k_allele_analyzer import (POPFILE, mrca_time, leaf_counts, backgrounds)  # noqa: E402


def samples_in_order():
    ids = []
    with open(POPFILE) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) >= 7:
                ids.append(p[1])
    return ids


def chunk_for(chrom, pos):
    OUT = ROOT / 'data/derived/g1k'
    for f in sorted(OUT.glob(f'chr{chrom}_*.dated.trees')):
        s = int(f.name.split('_')[1]); e = int(f.name.split('_')[2].split('.')[0])
        if s <= pos < e:
            return f, s, e
    raise FileNotFoundError(f'no dated chunk containing chr{chrom}:{pos}')


def anchors_for(g, start, cend, n_anchors, chrom, pos):
    if all(c in g.columns for c in ('start', 'end')):
        lo = max(int(g['start'].iloc[0]), start + 1)
        hi = min(int(g['end'].iloc[0]), cend - 1)
    else:
        lo = hi = int(pos)
    if hi <= lo:
        return [int(pos) - start]
    return [int(round(lo + i * (hi - lo) / (n_anchors - 1))) - start
            for i in range(n_anchors)]


def profile_locus(ts, start, rels, labels_by_hap, tau, young_age):
    """rows: per (label, anchor) statistics."""
    t = ts.tables.nodes.time
    samples = ts.samples()
    rows = []
    for i, rel in enumerate(rels):
        tree = ts.at(rel)
        par = tree.parent_array
        cnt = leaf_counts(tree, t, ts.num_nodes, samples)
        nodes = np.where(par != -1)[0]
        internal = np.unique(par[nodes])
        null_age = {}
        for v in internal:
            n = int(cnt[v])
            if n < 2:
                continue
            lb = int(np.floor(np.log2(n)))
            null_age.setdefault(lb, []).append(float(t[v]))
        null_age = {k: np.sort(np.array(v)) for k, v in null_age.items()}
        bg = backgrounds(par, t, samples, tau)
        ubg, bgcnt = np.unique(bg, return_counts=True)
        bg_size = {int(u): int(c) for u, c in zip(ubg, bgcnt)}
        for lab in labels_by_hap:
            carriers = np.where(labels_by_hap[lab] == lab)[0]
            n = len(carriers)
            if n < 3:
                rows.append(dict(allele=lab, anchor=i, pos_abs=start + rel,
                                 age=np.nan, age_pct=np.nan, n=n,
                                 n_backgrounds=np.nan, dom_share=np.nan,
                                 young=0))
                continue
            age = mrca_time(par, t, [samples[k] for k in carriers])
            lb = int(np.floor(np.log2(n)))
            pool = []
            for b in (lb - 1, lb, lb + 1):
                pool.extend(null_age.get(b, []))
            pool = np.array(pool)
            pct = 100 * float((pool < age).mean()) if len(pool) else np.nan
            bcar = bg[carriers]
            cnt_b = {}
            for b in bcar:
                cnt_b[int(b)] = cnt_b.get(int(b), 0) + 1
            dom = max(cnt_b, key=cnt_b.get)
            rows.append(dict(allele=lab, anchor=i, pos_abs=start + rel,
                             age=round(float(age), 1), age_pct=round(pct, 1), n=n,
                             n_backgrounds=len([b for b, c in cnt_b.items()
                                                if c >= max(1, 0.05 * n)]),
                             dom_share=round(cnt_b[dom] / n, 3),
                             young=1 if float(age) <= young_age else 0))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--table', default=str(ROOT / 'data/derived/old_dups/alleles_clusters.tsv'))
    ap.add_argument('--summary', default=str(ROOT / 'data/derived/old_dups/genome_wide/global_tau1000/summary.csv'))
    ap.add_argument('--pick', default=None, help='file of chrom,pos,allele lines to profile')
    ap.add_argument('--anchors', type=int, default=21)
    ap.add_argument('--tau', type=float, default=1000.0)
    ap.add_argument('--young-age', type=float, default=2000.0)
    ap.add_argument('--min-n', type=int, default=10)
    ap.add_argument('--freq-lo', type=float, default=0.005)
    ap.add_argument('--freq-hi', type=float, default=0.30)
    ap.add_argument('--max-age-min', type=float, default=1000.0)
    ap.add_argument('--out', default=str(ROOT / 'data/derived/old_dups/genome_wide/anchor_profiles.tsv'))
    a = ap.parse_args()

    # candidate selection
    if a.pick:
        cands = []
        with open(a.pick) as f:
            for line in f:
                p = line.rstrip('\n').split(',')
                if len(p) == 3:
                    cands.append((p[0], p[1], p[2]))
    else:
        cands = []
        for r in csv.DictReader(open(a.summary)):
            if (int(r['n']) >= a.min_n and float(r['age_min']) <= a.max_age_min
                    and a.freq_lo <= float(r['freq']) <= a.freq_hi):
                cands.append((r['chrom'], r['pos'], r['allele']))
    print(f'{len(cands)} candidate classes to profile')

    order = samples_in_order()
    idx = {s: i for i, s in enumerate(order)}
    df = pd.read_csv(a.table, sep=None, engine='python')
    for col in ('allele1', 'allele2'):
        df[col] = pd.to_numeric(df[col], errors='coerce')

    by_locus = {}
    for c, p, lab in cands:
        by_locus.setdefault((c, p), []).append(lab)

    out_rows = []
    for (chrom, pos), labels in sorted(by_locus.items(), key=lambda kv: (int(kv[0][0]), int(kv[0][1]))):
        g = df[(df['chrom'] == int(chrom)) & (df['pos'] == int(pos))]
        present = [s for s in order if s in set(g['sample'])]
        haps = [h for s in present for h in (2 * idx[s], 2 * idx[s] + 1)]
        f, start, cend = chunk_for(int(chrom), int(pos))
        ts = tskit.load(str(f)).simplify(samples=haps)
        sample_pos = {s: k for k, s in enumerate(present)}
        labels_by_hap = {}
        for lab in labels:
            arr = np.full(len(haps), -1, dtype=np.int64)
            for r in g.itertuples():
                k = sample_pos.get(r.sample)
                if k is None:
                    continue
                a1, a2 = r.allele1, r.allele2
                if (not pd.isna(a1)) and int(a1) == int(lab):
                    arr[2 * k] = int(lab)
                if (not pd.isna(a2)) and int(a2) == int(lab):
                    arr[2 * k + 1] = int(lab)
            labels_by_hap[int(lab)] = arr
        rels = anchors_for(g, start, cend, a.anchors, chrom, pos)
        rows = profile_locus(ts, start, rels, labels_by_hap, a.tau, a.young_age)
        for r in rows:
            r.update(chrom=int(chrom), pos=int(pos))
        out_rows.extend(rows)
        ny = {int(lab): sum(1 for r in rows if r['allele'] == int(lab) and r['young'])
              for lab in labels}
        print(f'chr{chrom}:{pos}: {len(labels)} labels, {len(rels)} anchors, '
              f'young-anchor counts {ny}', flush=True)

    out = pd.DataFrame(out_rows)
    out.to_csv(a.out, sep='\t', index=False)
    print(f'wrote {a.out}: {len(out)} rows')

    # summary: contiguity per candidate
    print('\ncandidate summary (young-age <= %.0f gen):' % a.young_age)
    print(f'{"locus":>22} {"allele":>6} {"n_young":>7} {"max_run":>8} {"span_frac":>9}  verdict')
    for (chrom, pos), labels in sorted(by_locus.items(), key=lambda kv: (int(kv[0][0]), int(kv[0][1]))):
        for lab in labels:
            pr = out[(out['chrom'] == int(chrom)) & (out['pos'] == int(pos))
                     & (out['allele'] == int(lab))].sort_values('anchor')
            y = list(pr['young'])
            n_young = sum(y)
            max_run = best = 0
            for v in y:
                best = best + 1 if v else 0
                max_run = max(max_run, best)
            span_frac = n_young / len(y) if len(y) else 0
            verdict = ('CONTIGUOUS' if max_run >= 3 and n_young >= 3
                       else 'scattered' if n_young >= 1 else 'never-young')
            print(f'chr{chrom}:{pos:>12} {lab:>6} {n_young:>7} {max_run:>8} {span_frac:>9.2f}  {verdict}')


if __name__ == '__main__':
    main()
