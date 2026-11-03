#!/usr/bin/env python3
"""Figure: tau-free breadth readout.

Panel A: genome-wide confinement regimes - all classes (age vs dom_share at
tau=1000) with the ROBUST risers highlighted, split by regime (distributed vs
confined) against the 0.5 confinement threshold.
Panel B: dom_share(tau) profiles of the 13 ROBUST risers with gold stars at
the interpolated crossover tau* (confinement threshold 0.5).

Reads data/derived/old_dups/genome_wide_erikarg/tau_star_profiles.csv and
data/derived/old_dups/genome_wide/rising_robustness.csv.
Usage: plot_tau_star.py [--out data/derived/old_dups/plots/tau_star.png]
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
GWE = ROOT / 'data/derived/old_dups/genome_wide_erikarg'
ROB = ROOT / 'data/derived/old_dups/genome_wide/rising_robustness.csv'
BED = Path('/moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed')
CMAP = {'not-confined<=1000': '#d73027', 'CONFINED<200': '#4575b4',
        'confined-in-profile': '#4575b4'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(ROOT / 'data/derived/old_dups/plots/tau_star.png'))
    a = ap.parse_args()
    rows = list(csv.DictReader(open(GWE / 'tau_star_profiles.csv')))
    for r in rows:
        for c in ('dom_share_200', 'dom_share_500', 'dom_share_1000'):
            r[c] = float(r[c])
        r['tau_star'] = float(r['tau_star']) if r['tau_star'] else None
    bed = {}
    for line in open(BED):
        f = line.split()
        if len(f) >= 4:
            c = f[0].split('#')[-1].replace('chr', '')
            bed[(int(c), (int(f[1]) + int(f[2])) // 2)] = f[3]
    rob = {}
    for r in csv.DictReader(open(ROB)):
        if r['robustness'] == 'ROBUST':
            rob[(r['chrom'], r['pos'], r['allele'])] = float(r['eur_rise'])
    rob_rows = sorted((r for r in rows if (r['chrom'], r['pos'], r['allele']) in rob),
                      key=lambda r: -rob[(r['chrom'], r['pos'], r['allele'])])

    fig = plt.figure(figsize=(16, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.35], wspace=0.25)
    ax = fig.add_subplot(gs[0])
    other = [r for r in rows if (r['chrom'], r['pos'], r['allele']) not in rob]
    ax.scatter([float(r['age']) for r in other], [r['dom_share_1000'] for r in other],
               s=3, alpha=0.15, color='0.6', edgecolor='none')
    dist = [r for r in rob_rows if r['status'] == 'not-confined<=1000']
    conf = [r for r in rob_rows if r['status'] != 'not-confined<=1000']
    ax.scatter([float(r['age']) for r in dist], [r['dom_share_1000'] for r in dist],
               s=70, color='#d73027', edgecolor='black', lw=0.5, zorder=5,
               label=f'ROBUST rising, distributed ({len(dist)})')
    ax.scatter([float(r['age']) for r in conf], [r['dom_share_1000'] for r in conf],
               s=70, color='#4575b4', marker='s', edgecolor='black', lw=0.5, zorder=5,
               label=f'ROBUST rising, confined ({len(conf)})')
    for r in rob_rows:
        g = bed.get((int(r['chrom']), int(r['pos'])), '?')
        ax.annotate(f"{g}\u00b7{r['allele']}", (float(r['age']), r['dom_share_1000']),
                    fontsize=7.5, xytext=(4, 4), textcoords='offset points')
    ax.axhline(0.5, color='black', ls='--', lw=0.8, alpha=0.6)
    ax.set_xscale('log')
    ax.set_xlabel('clade age (generations, log)')
    ax.set_ylabel('dom_share at tau=1000')
    ax.set_ylim(-0.03, 1.02)
    ax.legend(loc='upper left', fontsize=9)
    ax.set_title('A. Genome-wide confinement regimes (erik-ARG)', fontsize=11)

    axb = fig.add_subplot(gs[1])
    seen = set()
    for r in rob_rows:
        g = bed.get((int(r['chrom']), int(r['pos'])), '?')
        taus = [200, 500, 1000]
        ys = [r['dom_share_200'], r['dom_share_500'], r['dom_share_1000']]
        axb.plot(taus, ys, 'o-', color=CMAP[r['status']], lw=1.6, ms=4, alpha=0.85,
                 label=('_' if g in seen else f"{g}\u00b7{r['allele']}"))
        seen.add(g)
        if r['tau_star']:
            axb.plot([r['tau_star']], [0.5], marker='*', ms=13, color='gold',
                     markeredgecolor='black', zorder=6)
    axb.axhline(0.5, color='black', ls='--', lw=0.8, alpha=0.6)
    axb.set_xscale('log')
    axb.set_xticks([200, 500, 1000])
    axb.set_xticklabels(['200', '500', '1000'])
    axb.set_xlabel('tau (background depth, generations)')
    axb.set_ylabel('dom_share(tau)')
    axb.set_ylim(0, 1.0)
    axb.legend(fontsize=8, ncol=2, loc='upper left', framealpha=0.7)
    axb.set_title('B. dom_share(tau) profiles of the 13 ROBUST risers\n'
                  '(star = crossover tau*; red = distributed, blue = confined)', fontsize=11)
    fig.suptitle('Tau-free breadth readout: dom_share(tau) profiles and crossover tau* '
                 'per structural class', fontsize=13)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print('wrote', out)


if __name__ == '__main__':
    main()
