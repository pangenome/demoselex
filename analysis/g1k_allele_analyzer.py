#!/usr/bin/env python3
"""Per-locus analyzer for multiallelic diplotype data (segdup / CNV loci).

Input is NOT copy-number dosage. Each individual carries TWO alleles, each labelled with
an integer (e.g. 1..10, or 1..20); the diplotype (a,b) is the copy-number/haplotype
configuration at that locus. A file has one row per (locus, individual):

    chrom, pos, sample, allele1, allele2

For each locus and each distinct allele label, the analyzer reports, using the dated ARG
at that locus (no site injection needed - the tree is already there):

  freq              allele count / called haplotypes
  age               tsdate time of the MRCA of the allele's carriers
  age_pct           percentile of that age among clades of the same descendant-count bin
                    in the same local tree(s)  (empirical, within-locus null)
  n_backgrounds     number of ARG-derived local backgrounds (>=5% of carriers) holding it
  dom_share         fraction of the allele's carriers in its dominant background
  dom_freq          allele frequency inside that dominant background
  background_n      size of the dominant background (haplotypes)

Local background of a haplotype = most recent clade containing it whose root is younger
than tau (walk up until a node older than tau). Label-free: no population labels used.

Usage:
  g1k_allele_analyzer.py --table diplo.csv [--tau 200] [--out result.csv]
  g1k_allele_analyzer.py --selftest        # LCT SNP recoded as a 2-allele diplotype
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/g1k'
POPFILE = OUT / 'samples_population.txt'
TAU = 200.0


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
    for f in sorted(OUT.glob(f'chr{chrom}_*.dated.trees')):
        s = int(f.name.split('_')[1]); e = int(f.name.split('_')[2].split('.')[0])
        if s <= pos < e:
            return f, s
    raise FileNotFoundError(f'no dated chunk containing chr{chrom}:{pos}')


def leaf_counts(tree, t, n_nodes, samples):
    cnt = np.zeros(n_nodes, dtype=np.int32)
    cnt[samples] = 1
    par = tree.parent_array
    for v in np.argsort(t, kind='stable'):
        p = par[v]
        if p != -1:
            cnt[p] += cnt[v]
    return cnt


def mrca_time(par, t, nodes):
    counts = {}
    for u in nodes:
        v = int(u)
        while v != -1:
            counts[v] = counts.get(v, 0) + 1
            v = int(par[v])
    n = len(nodes)
    best = None
    for v, c in counts.items():
        if c == n and (best is None or t[v] < t[best]):
            best = v
    return (float(t[best]) if best is not None else np.nan)


def backgrounds(par, t, samples, tau):
    out = np.zeros(len(samples), dtype=np.int64)
    for k, s in enumerate(samples):
        v = int(s)
        while par[v] != -1 and (t[int(par[v])] - t[v]) < 1e-9 or (par[v] != -1 and t[int(par[v])] <= tau):
            if par[v] == -1:
                break
            v = int(par[v])
        out[k] = v
    return out


def analyse_locus(ts, start, rels, allele_by_hap, tau):
    """allele_by_hap: array over haplotype index (in ts sample order) -> allele label or -1.
    rels: anchor positions relative to the chunk start; per-label quantities are
    computed at each anchor and aggregated (median + range), because a structural
    locus spans many local trees and no single tree represents it."""
    t = ts.tables.nodes.time
    samples = ts.samples()
    labels = sorted({int(a) for a in allele_by_hap if a >= 0})
    called_haps = int((allele_by_hap >= 0).sum())
    per_anchor = []
    for rel in rels:
        tree = ts.at(rel)
        par = tree.parent_array
        cnt = leaf_counts(tree, t, ts.num_nodes, samples)
        # within-locus null: ages of all clades by leaf-count bin
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
        # local backgrounds
        bg = backgrounds(par, t, samples, tau)
        ubg, bgcnt = np.unique(bg, return_counts=True)
        bg_size = {int(u): int(c) for u, c in zip(ubg, bgcnt)}
        rows_at = {}
        for lab in labels:
            carriers = np.where(allele_by_hap == lab)[0]
            n = len(carriers)
            if n < 3:
                continue
            age = mrca_time(par, t, [samples[i] for i in carriers])
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
            broad = [b for b, c in cnt_b.items() if c >= max(1, 0.05 * n)]
            dom = max(cnt_b, key=cnt_b.get)
            rows_at[lab] = dict(
                age=round(float(age), 1) if np.isfinite(age) else np.nan,
                age_pct=round(pct, 1),
                n_backgrounds=len(broad), n_backgrounds_all=len(cnt_b),
                dom_share=round(cnt_b[dom] / n, 3),
                dom_freq=round(cnt_b[dom] / max(bg_size.get(dom, 0), 1), 3),
                background_n=bg_size.get(dom, 0))
        per_anchor.append(rows_at)
    rows = []
    for lab in labels:
        n = int((allele_by_hap == lab).sum())
        if n < 3:
            continue
        sts = [d[lab] for d in per_anchor if lab in d]
        if not sts:
            continue
        med = lambda k: float(np.median([s[k] for s in sts if not pd.isna(s[k])]))
        ages = [s['age'] for s in sts if not pd.isna(s.get('age', np.nan))]
        rows.append(dict(allele=lab, n=n, freq=round(n / max(called_haps, 1), 4),
                          age=round(med('age'), 1),
                          age_min=round(min(ages), 1) if ages else np.nan,
                          age_max=round(max(ages), 1) if ages else np.nan,
                          age_pct=round(med('age_pct'), 1),
                          n_backgrounds=round(med('n_backgrounds')),
                          n_backgrounds_all=round(med('n_backgrounds_all')),
                          dom_share=round(med('dom_share'), 3),
                          dom_freq=round(med('dom_freq'), 3),
                          background_n=round(med('background_n')),
                          n_anchors=len(sts)))
    return rows


def run(table, tau, n_anchors=1):
    df = pd.read_csv(table, sep=None, engine='python')   # sniff comma or tab
    need = {'chrom', 'pos', 'sample', 'allele1', 'allele2'}
    if not need.issubset(df.columns):
        raise SystemExit(f'table needs columns {sorted(need)}; got {list(df.columns)}')
    for col in ('allele1', 'allele2'):
        df[col] = pd.to_numeric(df[col], errors='coerce')  # NA -> NaN
    order = samples_in_order()
    idx = {s: i for i, s in enumerate(order)}
    out_rows = []
    for (chrom, pos), g in df.groupby(['chrom', 'pos']):
        # restrict the ARG to the cohort present in the table, keeping VCF order
        present = [s for s in order if s in set(g['sample'])]
        if not present:
            print(f'chr{chrom}:{pos}: no samples match the ARG cohort'); continue
        haps = [h for s in present for h in (2 * idx[s], 2 * idx[s] + 1)]
        f, start = chunk_for(int(chrom), int(pos))
        end = int(f.name.split('_')[2].split('.')[0])
        ts = tskit.load(str(f)).simplify(samples=haps)
        # map: haplotype position in ts -> allele label (NA/-1 = missing)
        sample_pos = {s: k for k, s in enumerate(present)}
        allele_by_hap = np.full(len(haps), -1, dtype=np.int64)
        for r in g.itertuples():
            k = sample_pos.get(r.sample)
            if k is None:
                continue
            a1, a2 = r.allele1, r.allele2
            allele_by_hap[2 * k] = -1 if pd.isna(a1) else int(a1)
            allele_by_hap[2 * k + 1] = -1 if pd.isna(a2) else int(a2)
        # anchors: spread across the locus span (clamped into this chunk)
        if n_anchors > 1 and all(c in g.columns for c in ('start', 'end')):
            lo = max(int(g['start'].iloc[0]), start + 1)
            hi = min(int(g['end'].iloc[0]), end - 1)
            if hi > lo:
                rels = [int(round(lo + i * (hi - lo) / (n_anchors - 1))) - start
                        for i in range(n_anchors)]
            else:
                rels = [int(pos) - start]
        else:
            rels = [int(pos) - start]
        rows = analyse_locus(ts, start, rels, allele_by_hap, tau)
        for r in rows:
            r.update(chrom=int(chrom), pos=int(pos), n_samples=len(present))
            out_rows.append(r)
        print(f'chr{chrom}:{pos}: {len(present)} samples, alleles {sorted(set(allele_by_hap[allele_by_hap>=0]))}', flush=True)
    res = pd.DataFrame(out_rows)
    return res


def selftest(tau):
    """Recode LCT rs4988235 (ALT allele) as a 2-allele diplotype and check the known result."""
    order = samples_in_order()
    eur = []
    with open(POPFILE) as f:
        next(f)
        for i, line in enumerate(f):
            p = line.split()
            if len(p) >= 7 and p[6] == 'EUR':
                eur.append(p[1])
    chrom, pos = 2, 135851076
    f, start = chunk_for(chrom, pos)
    haps = [h for s in eur for h in (2 * order.index(s), 2 * order.index(s) + 1)]
    ts = tskit.load(str(f)).simplify(samples=haps)
    j = int(np.argmin(np.abs(np.array([s.position for s in ts.sites()]) - (pos - start))))
    G = ts.genotype_matrix()[j]
    rows = []
    for k, s in enumerate(eur):
        a1 = int(G[2 * k]); a2 = int(G[2 * k + 1])
        # encode allele state 2 (=1 in the ARG's ALT orientation) as label 2
        rows.append(dict(chrom=chrom, pos=pos, sample=s,
                         allele1=a1 + 1, allele2=a2 + 1))
    tmp = OUT / 'selftest_diplo.csv'
    pd.DataFrame(rows).to_csv(tmp, index=False)
    res = run(str(tmp), tau)
    print('\nSELFTEST (LCT rs4988235 recoded as a 2-allele diplotype in EUR):')
    print(res.to_string(index=False))
    print('\nexpected: label 2 (ALT, freq ~0.54) young clade, confined to one background, '
          'high within-background frequency; label 1 (REF) the complement.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--table')
    ap.add_argument('--tau', type=float, default=TAU)
    ap.add_argument('--anchors', type=int, default=1,
                    help='sample this many trees across each locus span and aggregate')
    ap.add_argument('--out', default=None)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        selftest(a.tau); return
    if not a.table:
        raise SystemExit('need --table or --selftest')
    res = run(a.table, a.tau, a.anchors)
    out = a.out or (OUT / 'allele_locus_results.csv')
    res.to_csv(out, index=False)
    print(f'\nwrote {out}')
    print(res.to_string(index=False))


if __name__ == '__main__':
    main()
