#!/usr/bin/env python3
"""Self-contained tests for the loci-ARG pipeline helpers.

No large data required: everything runs on synthetic fixtures in /tmp or in
memory. Run with:

    .venv/bin/python analysis/test_loci_pipeline.py

Data-dependent validation (the LCT anchor against a built chunk) lives in the
analyzer itself: `g1k_allele_analyzer.py --selftest` (needs the dated chunks).
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
FAILS = []


def check(name, cond, detail=''):
    status = 'ok' if cond else 'FAIL'
    print(f'  [{status}] {name}' + (f'  ({detail})' if detail and not cond else ''))
    if not cond:
        FAILS.append(name)


# ---------------------------------------------------------------- g1k_loci_build
def test_loci_build():
    print('g1k_loci_build:')
    from g1k_loci_build import norm_chrom, chunks_for_bed
    check('norm_chrom: GRCh38-style contig', norm_chrom('GRCh38#0#chr7') == 7)
    check('norm_chrom: chr prefix', norm_chrom('chrX'[:0] + 'chr12') == 12)
    check('norm_chrom: bare number', norm_chrom('3') == 3)
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
        f.write('#comment\ntrack name=x\n')
        f.write('GRCh38#0#chr1\t15000000\t15900000\tPOTEH\n')
        f.write('GRCh38#0#chr1\t29000000\t31000000\tX\n')   # spans two chunks
        f.write('GRCh38#0#chr22\t15326560\t15928214\tP\n')
        path = f.name
    rows = chunks_for_bed(path, with_selftest=False)
    starts = {(c, s) for c, s, e in rows}
    check('grid: locus inside one 10 Mb chunk', (1, 10000000) in starts, str(sorted(starts)))
    check('grid: boundary-spanning locus gets both chunks',
          (1, 20000000) in starts and (1, 30000000) in starts, str(sorted(starts)))
    check('grid: chunk end clipped at contig length',
          any(c == 22 and s == 10000000 and e == 20000000 for c, s, e in rows))
    check('grid: no off-grid chunk starts',
          all(s % 10_000_000 == 0 for _, s, _ in rows))
    from g1k_loci_build import CLEN
    check('CLEN: chr20 is the corrected 64,444,167', CLEN[20] == 64_444_167)
    check('CLEN: all 22 autosomes present', sorted(CLEN) == list(range(1, 23)))


# ------------------------------------------------------------- build_input_tables
def test_input_tables():
    print('build_input_tables:')
    from build_input_tables import path_key, cluster_labels, load_bed
    check('path_key: interval suffix stripped',
          path_key('HG01934#1#CM088841.1:48425716-48827166') == 'HG01934#1#CM088841.1')
    check('path_key: NA passthrough', path_key('NA') is None)
    check('path_key: no colon = as-is', path_key('CHM13#0#chr7') == 'CHM13#0#chr7')

    # label scheme: clusters by descending size 1..k, then private labels
    cl = {'a': 0, 'b': 0, 'c': 0, 'd': 1, 'e': 1, 'f': 2, 'g': -1}
    relab, k, qc = cluster_labels('X', set(cl), {'X': cl},
                                  {'X': dict(k='3', silhouette='0.5', noise_frac='0.1')})
    check('labels: cluster 0 (size 3) -> 1', relab['a'] == relab['b'] == relab['c'] == 1)
    check('labels: cluster 1 (size 2) -> 2', relab['d'] == relab['e'] == 2)
    check('labels: cluster 2 (size 1) -> 3', relab['f'] == 3)
    check('labels: noise path gets private label > k', relab['g'] == 4 and k == 3)
    check('labels: qc tuple carried through', qc == ('3', '0.5', '0.1'))

    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
        f.write('GRCh38#0#chr4\t100\t300\tGYP\n')
        bpath = f.name
    bed = load_bed(bpath)
    check('load_bed: midpoint + span', bed['GYP'] == (4, 200, 100, 300), str(bed))


# ---------------------------------------------------------------- allele analyzer
def test_analyzer_core():
    print('g1k_allele_analyzer (tree helpers on a synthetic tree):')
    import tskit
    from g1k_allele_analyzer import mrca_time, backgrounds, leaf_counts

    def synth_ts(n=8):
        t = tskit.TableCollection(sequence_length=1000)
        for _ in range(n):
            t.nodes.add_row(flags=1, time=0.0)
        active = list(range(n))
        time, node = 1.0, n
        while len(active) > 1:
            nxt = []
            for i in range(0, len(active) - 1, 2):
                a, b = active[i], active[i + 1]
                t.nodes.add_row(flags=0, time=time)
                t.edges.add_row(parent=node, child=a, left=0, right=1000)
                t.edges.add_row(parent=node, child=b, left=0, right=1000)
                nxt.append(node)
                node += 1
            if len(active) % 2:
                nxt.append(active[-1])
            active = nxt
            time += 1.0
        t.sort()
        return t.tree_sequence()

    ts = synth_ts(8)
    tree = ts.first()
    par = tree.parent_array
    t = ts.tables.nodes.time
    samples = ts.samples()
    # MRCA of two sister tips = their parent
    s0, s1 = int(samples[0]), int(samples[1])
    if par[s0] == par[s1]:
        check('mrca_time: sister tips -> parent node time',
              mrca_time(par, t, [s0, s1]) == float(t[par[s0]]))
    # MRCA of all samples = root
    root = tree.root
    check('mrca_time: all samples -> root time',
          mrca_time(par, t, list(samples)) == float(t[root]))
    # leaf counts
    cnt = leaf_counts(tree, t, ts.num_nodes, samples)
    check('leaf_counts: root holds all samples', cnt[root] == ts.num_samples)
    check('leaf_counts: each sample holds 1', all(cnt[s] == 1 for s in samples))
    # backgrounds: tau=0 -> every sample is its own background
    bg = backgrounds(par, t, samples, 0.0)
    check('backgrounds: tau=0 -> singleton backgrounds', len(set(bg)) == len(samples))
    # backgrounds: tau=inf -> one background (the root clade) per tree
    bg = backgrounds(par, t, samples, 1e9)
    check('backgrounds: tau=inf -> single background', len(set(bg)) == 1)
    # multi-anchor wrapper: analyse_locus on a label present in >3 haps
    from g1k_allele_analyzer import analyse_locus
    hap = np.full(ts.num_samples, -1, dtype=np.int64)
    hap[:6] = 5                       # a 6-carrier class
    rows = analyse_locus(ts, 0, [0, 500], hap, 200.0)
    r = [x for x in rows if x['allele'] == 5]
    check('analyse_locus: multi-anchor row produced with n and freq',
          len(r) == 1 and r[0]['n'] == 6 and abs(r[0]['freq'] - 1.0) < 1e-9)
    check('analyse_locus: aggregates both anchors',
          r and r[0]['n_anchors'] == 2)


def main():
    test_loci_build()
    test_input_tables()
    test_analyzer_core()
    print()
    if FAILS:
        print(f'{len(FAILS)} FAILED: {FAILS}')
        sys.exit(1)
    print('all tests passed')


if __name__ == '__main__':
    main()
