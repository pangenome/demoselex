#!/usr/bin/env python3
"""Figures for the temporal-layer results.

  trajectories_robust.png   13 ROBUST candidates: pooled-EUR trajectory (bold,
                            coloured) + per-country thin lines; frequency vs
                            date BP. The class id and gene annotate each panel.
  rise_overview.png         all EUR-RISING candidates: oldest-bin vs modern
                            frequency, coloured by robustness, annotated genes
                            for the ROBUST set.
  top6_panels.png           the six largest robust rises, enlarged, with the
                            modern EUR endpoint and drift-null ribbon.

Usage: plot_results.py [--out-dir data/derived/old_dups/plots]
"""
import argparse
import csv
import glob
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
GW = ROOT / 'data/derived/old_dups/genome_wide'
COUNTRIES = ['Hungary', 'Russia', 'Belgium', 'United_Kingdom', 'Italy', 'Sweden', 'Denmark']
CCOL = {'Hungary': '#b2182b', 'Russia': '#2166ac', 'Belgium': '#1b7837',
        'United_Kingdom': '#762a83', 'Italy': '#e08214', 'Sweden': '#4393c3',
        'Denmark': '#8c510a'}


def bed_names():
    bed = {}
    for line in open(Path('/moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed')):
        f = line.split()
        if len(f) >= 4:
            c = f[0].split('#')[-1].replace('chr', '')
            bed[(int(c), (int(f[1]) + int(f[2])) // 2)] = f[3]
    return bed


def traj_points(r):
    pts = [(int(s.split(':')[0]), float(s.split(':')[1]))
           for s in r['freqs_by_bin'].split(';')]
    pts.sort(key=lambda t: -t[0])          # oldest first
    return pts


def load(key):
    """pooled EUR row + per-country rows for (chrom,pos,allele)."""
    pooled = None
    for r in csv.DictReader(open(GW / 'trajectories_EUR.csv')):
        if (r['chrom'], r['pos'], r['allele']) == key:
            pooled = r
    per = {}
    for c in COUNTRIES:
        try:
            for r in csv.DictReader(open(GW / f'by_country/{c}.csv')):
                if (r['chrom'], r['pos'], r['allele']) == key:
                    per[c] = r
        except FileNotFoundError:
            pass
    return pooled, per


def robust_list():
    return [r for r in csv.DictReader(open(GW / 'rising_robustness.csv'))
            if r['robustness'] == 'ROBUST']


def panel(ax, pooled, per, gene, allele, big=False):
    for c, r in per.items():
        pts = traj_points(r)
        if len(pts) >= 2:
            xs = [-p[0] / 1000 for p in pts] + [0]
            ys = [p[1] for p in pts] + [float(r['freq_modern'])]
            ax.plot(xs, ys, '-', color=CCOL.get(c, 'grey'), lw=0.8, alpha=0.7)
    if pooled:
        pts = traj_points(pooled)
        xs = [-p[0] / 1000 for p in pts] + [0]
        ys = [p[1] for p in pts] + [float(pooled['freq_modern'])]
        ax.plot(xs, ys, 'o-', color='black', lw=2.2, ms=3.5, zorder=5)
        # drift ribbon around the oldest-bin frequency
        f0 = float(pooled['freq_first'])
        if f0 > 0:
            import math
            dt = pts[0][0] / 28.0
            for ne in (5000, 10000):
                sd = math.sqrt(f0 * (1 - f0) * dt / (2 * ne))
                ax.fill_between([-xs[0], 0], f0 - 2 * sd, f0 + 2 * sd,
                                color='grey', alpha=0.10, zorder=1)
    ax.set_title(f'{gene} class {allele}', fontsize=10)
    ax.set_ylim(-0.04, 1.0)
    if big:
        ax.set_xlabel('kiloyears BP')
        ax.set_ylabel('class frequency')
    ax.tick_params(labelsize=8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out-dir', default=str(ROOT / 'data/derived/old_dups/plots'))
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    bed = bed_names()

    rows = robust_list()
    rows.sort(key=lambda r: -float(r['eur_rise']))
    n = len(rows)
    ncol = 4
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 2.6 * nrow))
    for i, r in enumerate(rows):
        ax = axes[i // ncol, i % ncol] if nrow > 1 else axes[i % ncol]
        pooled, per = load((r['chrom'], r['pos'], r['allele']))
        panel(ax, pooled, per, r['gene'], r['allele'])
    for j in range(n, nrow * ncol):
        fig.delaxes(axes[j // ncol, j % ncol] if nrow > 1 else axes[j % ncol])
    handles = [plt.Line2D([], [], color='black', lw=2.2, marker='o', ms=3.5,
               label='pooled EUR'),
               plt.Line2D([], [], color='grey', alpha=0.3, lw=6, label='drift null (2sd)')]
    handles += [plt.Line2D([], [], color=CCOL[c], lw=1.0, alpha=0.8,
               label=c.replace('_', ' ')) for c in COUNTRIES]
    fig.legend(handles=handles, loc='lower center', ncol=9, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle('ROBUST rising structural classes: pooled EUR (bold) vs per-country series',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(out / 'trajectories_robust.png', dpi=150)
    print('wrote', out / 'trajectories_robust.png')

    top6 = rows[:6]
    fig, axes = plt.subplots(2, 3, figsize=(13, 6))
    for i, r in enumerate(top6):
        ax = axes[i // 3, i % 3]
        pooled, per = load((r['chrom'], r['pos'], r['allele']))
        panel(ax, pooled, per, r['gene'], r['allele'], big=True)
        ax.set_title(f"{r['gene']} class {r['allele']}  "
                     f"({float(r['eur_rise']):+.2f}, {r['n_countries_rising']})", fontsize=11)
    handles = [plt.Line2D([], [], color='black', lw=2.2, marker='o', ms=3.5,
               label='pooled EUR'),
               plt.Line2D([], [], color='grey', alpha=0.3, lw=6, label='drift null (2sd)')]
    handles += [plt.Line2D([], [], color=CCOL[c], lw=1.0, alpha=0.8,
               label=c.replace('_', ' ')) for c in COUNTRIES]
    fig.legend(handles=handles, loc='lower center', ncol=9, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle('Six largest robust rises (grey ribbon: 2-sd drift null, Ne 5k-10k)',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(out / 'top6_panels.png', dpi=150)
    print('wrote', out / 'top6_panels.png')

    # overview scatter: all EUR candidates
    fig, ax = plt.subplots(figsize=(8, 7))
    colmap = {'ROBUST': '#d73027', 'single-country': '#fdae61', 'pooled-only': '#7fbf7b'}
    allr = list(csv.DictReader(open(GW / 'rising_robustness.csv')))
    for r in allr:
        pooled, _ = load((r['chrom'], r['pos'], r['allele']))
        if not pooled:
            continue
        x = float(pooled['freq_first'])
        y = float(pooled['freq_modern'])
        col = colmap[r['robustness']]
        ax.scatter(x, y, s=40, color=col, alpha=0.8,
                   edgecolor='black', linewidth=0.4, zorder=3)
        if r['robustness'] == 'ROBUST':
            ax.annotate(f"{r['gene']}·{r['allele']}", (x, y), fontsize=7,
                        xytext=(3, 3), textcoords='offset points')
    ax.plot([0, 1], [0, 1], '--', color='grey', lw=0.8)
    ax.set_xlabel('frequency in oldest ancient bin')
    ax.set_ylabel('frequency in modern EUR')
    ax.set_xlim(-0.02, 0.6); ax.set_ylim(-0.02, 1.0)
    for lab, col in colmap.items():
        ax.scatter([], [], color=col, label=lab)
    ax.legend(loc='lower right', fontsize=9)
    ax.set_title('EUR rising candidates: ancient vs modern frequency')
    fig.tight_layout()
    fig.savefig(out / 'rise_overview.png', dpi=150)
    print('wrote', out / 'rise_overview.png')


if __name__ == '__main__':
    main()
