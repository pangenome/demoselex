#!/usr/bin/env python3
"""22-page atlas: class composition through time + tau* confinement per locus.

Per locus panel (256 duplicated/CNV loci, 12 per page):
- top-7 classes by modern frequency, through time (ancient bins -> modern),
  dashed line = remainder ("other classes")
- ROBUST risers drawn as stars, red = distributed regime (never confined),
  blue = confined regime (CONFINED<200 or confined-in-profile from the
  erik-ARG tau_star_profiles.csv)
- title: gene (span kb), * = locus holds a robust riser, [HARD] = ancient/modern
  class-distribution overlap < 0.5, and the locus's confinement summary
  "confined k/N tau*" from the tau* profiles

Reads: genome_wide/trajectories_EUR.csv, rising_robustness.csv,
genome_wide_erikarg/tau_star_profiles.csv, alleles_clusters.tsv,
/moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed
Usage: plot_atlas_v2.py [--out data/derived/old_dups/plots/all_loci_composition_v2.pdf]
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
GW = ROOT / 'data/derived/old_dups/genome_wide'
GWE = ROOT / 'data/derived/old_dups/genome_wide_erikarg'
TABLE = ROOT / 'data/derived/old_dups/alleles_clusters.tsv'
BED = Path('/moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed')
NCOL, NROW = 4, 3
PER = NCOL * NROW


def load():
    bed, span = {}, {}
    for line in open(BED):
        f = line.split()
        if len(f) >= 4:
            c = f[0].split('#')[-1].replace('chr', '')
            bed[(int(c), (int(f[1]) + int(f[2])) // 2)] = f[3]
            span[f[3]] = (int(f[2]) - int(f[1])) / 1000
    rows = defaultdict(list)
    for r in csv.DictReader(open(GW / 'trajectories_EUR.csv')):
        g = bed.get((int(r['chrom']), int(r['pos'])))
        if g:
            rows[g].append(r)
    rob = {}
    for r in csv.DictReader(open(GW / 'rising_robustness.csv')):
        if r['robustness'] == 'ROBUST':
            rob[(r['chrom'], r['pos'], r['allele'])] = r
    tstar = {}
    loc_conf = defaultdict(lambda: [0, 0])
    for r in csv.DictReader(open(GWE / 'tau_star_profiles.csv')):
        tstar[(r['chrom'], r['pos'], r['allele'])] = r
        loc_conf[(r['chrom'], r['pos'])][1] += 1
        if r['status'] in ('CONFINED<200', 'confined-in-profile'):
            loc_conf[(r['chrom'], r['pos'])][0] += 1
    anc, mod = defaultdict(lambda: defaultdict(int)), defaultdict(lambda: defaultdict(int))
    for r in csv.DictReader(open(TABLE), delimiter='\t'):
        key = (r['chrom'], r['pos'])
        tgt = mod[key] if r['sample'][:2] in ('HG', 'NA') else anc[key]
        for a in (r['allele1'], r['allele2']):
            if a != 'NA':
                tgt[a] += 1

    def overlap(key):
        A, M = anc[key], mod[key]
        na, nm = sum(A.values()), sum(M.values())
        if not na or not nm:
            return 1.0
        return sum(min(c / na, M.get(k, 0) / nm) for k, c in A.items())

    return bed, span, rows, rob, tstar, loc_conf, overlap


def draw_panel(ax, g, D, cmap):
    bed, span, rows, rob, tstar, loc_conf, overlap = D
    key0 = [k for k, v in bed.items() if v == g][0]
    cls = sorted(rows[g], key=lambda r: -float(r['freq_modern'] or 0))[:7]
    other = None
    for i, r in enumerate(cls):
        pts = sorted([(int(s.split(':')[0]), float(s.split(':')[1]))
                      for s in r['freqs_by_bin'].split(';')], key=lambda t: -t[0])
        xs = [-x / 1000 for x, _ in pts] + [0]
        ys = [y for _, y in pts] + [float(r['freq_modern'])]
        rk = (r['chrom'], r['pos'], r['allele'])
        if rk in rob:
            regime = tstar.get(rk, {}).get('status', 'not-confined<=1000')
            col = '#4575b4' if regime != 'not-confined<=1000' else '#d73027'
            ax.plot(xs, ys, '*', color=col, lw=1.6, ms=11, markeredgecolor='black',
                    markeredgewidth=0.4, label=f"{r['allele']}\u2605")
        else:
            ax.plot(xs, ys, 'o', color=cmap(i % 10), lw=0.8, ms=2.5, alpha=0.55,
                    label=r['allele'])
        if other is None:
            other = {x: 1.0 - y for x, y in zip(xs, ys)}
        else:
            for x, y in zip(xs, ys):
                other[x] = other.get(x, 1.0) - y
    if other:
        xs = sorted(other)
        ax.plot(xs, [max(other[x], 0) for x in xs], '--', color='0.55', lw=1.0)
    ov = overlap(key0)
    nc, nt = loc_conf.get((str(key0[0]), str(key0[1])), [0, 0])
    hard = ' [HARD]' if ov < 0.5 else ''
    robm = ' *' if any((r['chrom'], r['pos'], r['allele']) in rob
                       for r in rows[g] if r['pos'] == str(key0[1])) else ''
    ax.set_title(f'{g} ({span[g]:.0f} kb){robm}{hard} | confined {nc}/{nt} tau*',
                 fontsize=7.5, color='#b30000' if ov < 0.5 else 'black')
    ax.set_ylim(-0.04, 1.0)
    ax.set_xlim(-14, 1)
    ax.tick_params(labelsize=6)
    if cls:
        ax.legend(fontsize=5, ncol=4, loc='upper left', framealpha=0.4,
                  title='class (\u2605 robust)', title_fontsize=5, handlelength=1.0,
                  columnspacing=0.8, handletextpad=0.3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(ROOT / 'data/derived/old_dups/plots/all_loci_composition_v2.pdf'))
    a = ap.parse_args()
    D = load()
    bed, span, rows, rob = D[0], D[1], D[2], D[3]
    rob_loci = {r['gene'] for r in rob.values()}
    order = sorted(rows, key=lambda g: [k for k, v in bed.items() if v == g][0])
    cmap = plt.cm.tab10
    npages = (len(order) + PER - 1) // PER
    with PdfPages(a.out) as pdf:
        for p in range(npages):
            page = order[p * PER:(p + 1) * PER]
            fig, axes = plt.subplots(NROW, NCOL, figsize=(15, 10))
            for ax, g in zip(axes.flat, page):
                draw_panel(ax, g, D, cmap)
            for j in range(len(page), PER):
                fig.delaxes(axes.flat[j])
            for ax in axes.flat[::NCOL]:
                if ax.get_visible():
                    ax.set_ylabel('class frequency', fontsize=8)
            for ax in axes[-1]:
                if ax.get_visible():
                    ax.set_xlabel('kya BP', fontsize=8)
            fig.suptitle(
                f'Class composition through time + tau* confinement (erik-ARG statics) '
                f'- page {p + 1}/{npages}; \u2605 robust riser '
                f'(red=distributed, blue=confined); dashed = other classes; '
                f'[HARD] = low call-distribution overlap; * = locus holds a robust riser',
                fontsize=9)
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            pdf.savefig(fig)
            plt.close(fig)
    print(f'wrote {a.out}: {npages} pages, {len(order)} loci')


if __name__ == '__main__':
    main()
