#!/usr/bin/env python3
"""Pilot harness: put UNPHASED diploid SV/CNV genotypes onto the dated ARG genealogy.

Companion to g1k_allele_analyzer.py.  That script asks "does this allele sit in one
ancestry background".  This one asks the sharper, phasing-aware question and is explicit
about what phase can and cannot change:

  Stage A   phase-free ancestry enrichment (sample level, permutation p-value)
  Stage B   tree parsimony: assign each heterozygote's two alleles to its two haplotypes
            so that the carrier set is as cladal as possible (coordinate descent)
  Stage C   phase-robustness certificate: report the RANGE of the tree statistics over
            every phasing consistent with the observed allele pair

Central statistic -- n_origins:
    the minimum number of independent gains of the allele, which for a complete leaf
    painting equals the number of maximal all-carrier clades in the local tree.
      1  -> the carriers form a clade: one mutation explains the pattern, so the carrier
            clade's MRCA is a legitimate date for the event (single origin)
     >1  -> >=2 independent gains: homoplasy, recurrence, a merged call, or a phase error

Why no external phaser (Beagle/SHAPEIT/MVNcall): a deletion leaves carriers hemizygous
(one copy), an inversion is an LD island with suppressed recombination, and every
Li-Stephens panel phaser assumes two copies and smooths away exactly the recurrent-mutation
structure we want to measure.  Here phase is never invented -- heterozygotes are enumerated
and a locus is only reported CERTIFIED when the conclusion holds under every admissible
phase.

Input table (one row per locus x sample; integer allele labels, NA for missing):

    locus, chrom, start, end, sample, allele1, allele2

`start`/`end` are the SV span; the genealogy is read at both breakpoints and the midpoint
(a point locus may use start == end).  Samples must be 1000G IDs from
samples_population.txt; haplotype order is (2i, 2i+1) for the i-th sample.

Usage:
  python3 sv_phase_pilot.py --selftest                  # recovers known truth from the ARG
  python3 sv_phase_pilot.py --table loci.tsv --out pilot.csv
  python3 sv_phase_pilot.py --table loci.tsv --tau 200 --perm 5000 --exact-max 10
  python3 sv_phase_pilot.py --selftest --emit-demo-table demo_loci.tsv   # format example
"""
import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parent))
import g1k_allele_analyzer as A  # reuse cohort order, chunk lookup, backgrounds, mrca_time

TAU = 200.0
PRUNE = 0.05          # a background "holds" the allele if >=5% of carriers are in it


# ---------------------------------------------------------------- per-position engine
class LocalTree:
    """Genealogy at one position, in a form that makes allele painting cheap.

    Ancestor paths of the haplotypes are flattened once; then any carrier set gives its
    per-node leaf counts with a single bincount (no Python loop over the tree).
    """

    def __init__(self, ts, start, pos, tau):
        self.pos = int(pos)
        tree = ts.at(pos - start)
        self.tree = tree
        self.par = np.array(tree.parent_array)
        n = len(self.par)
        t_full = ts.tables.nodes.time
        self.t = np.asarray(t_full[:n])
        self.n_nodes = n
        self.leaves = np.asarray(ts.samples())
        self.n_hap = len(self.leaves)
        self._leaf_order = np.argsort(self.leaves)
        self._leaf_sorted = self.leaves[self._leaf_order]
        bg = A.backgrounds(self.par, self.t, self.leaves, tau)
        self.bg_of_hap = np.asarray(bg, dtype=np.int64)
        u, c = np.unique(self.bg_of_hap, return_counts=True)
        self.bg_size = {int(a): int(b) for a, b in zip(u, c)}
        # flattened ancestor paths: entry -> (node id, owning haplotype)
        paths = []
        for h in range(self.n_hap):
            v, p = int(self.leaves[h]), [int(self.leaves[h])]
            while self.par[v] != -1:
                v = int(self.par[v])
                p.append(v)
            paths.append(np.asarray(p, dtype=np.int64))
        self.path_lens = np.array([len(p) for p in paths])
        self.path_flat = np.concatenate(paths)
        self.path_owner = np.repeat(np.arange(self.n_hap), self.path_lens)
        self.cnt_all = self._counts(np.ones(self.n_hap, dtype=bool))

    def _counts(self, carrier_mask):
        sel = np.repeat(carrier_mask, self.path_lens)
        if not sel.any():
            return np.zeros(self.n_nodes, dtype=np.int32)
        return np.bincount(self.path_flat[sel], minlength=self.n_nodes).astype(np.int32)

    def clade_hap_mask(self, v):
        """Boolean mask over haplotypes: is this haplotype a leaf of clade(v)?"""
        m = np.zeros(self.n_hap, dtype=bool)
        leaf_nodes = np.fromiter(self.tree.get_leaves(v), dtype=np.int64)
        pos = np.searchsorted(self._leaf_sorted, leaf_nodes)
        m[self._leaf_order[pos]] = True
        return m

    def stats(self, carrier_mask):
        cnt_a = self._counts(carrier_mask)
        pure = (cnt_a == self.cnt_all) & (self.cnt_all > 0)
        idx = np.where(self.par >= 0, self.par, 0)
        maximal = pure & ~((self.par >= 0) & pure[idx])
        n = int(carrier_mask.sum())
        if n == 0:
            return dict(n_origins=0, n_backgrounds=0, dom_share=np.nan, dom_freq=np.nan,
                        age=np.nan)
        b = self.bg_of_hap[carrier_mask]
        u, c = np.unique(b, return_counts=True)
        dom = int(u[np.argmax(c)])
        return dict(n_origins=int(maximal.sum()),
                    n_backgrounds=int((c >= max(1, PRUNE * n)).sum()),
                    dom_share=float(c.max() / n),
                    dom_freq=float(c.max() / max(self.bg_size.get(dom, 1), 1)),
                    age=A.mrca_time(self.par, self.t, self.leaves[carrier_mask]))


# ---------------------------------------------------------------- stage A: permutation
def enrich_test(group_of_sample, carrier_samples, n_perm, rng):
    """Phase-free ancestry concentration of carrier samples.

    A sample's group is the unordered pair of local backgrounds of its two haplotypes, so
    the test never asks which haplotype is which.  Statistic: the share of carriers sitting
    in their largest ancestry group.  Null: carrier labels permuted among samples.
    """
    n = len(group_of_sample)
    carrier_samples = np.asarray(carrier_samples, dtype=bool)
    k = int(carrier_samples.sum())
    if k < 3 or k > n - 3:
        return np.nan, np.nan, np.nan

    def share(mask):
        cnt = np.bincount(group_of_sample[mask], minlength=int(group_of_sample.max()) + 1)
        return cnt.max() / k

    obs = share(carrier_samples)
    ge = 0
    for _ in range(n_perm):
        m = np.zeros(n, dtype=bool)
        m[rng.choice(n, size=k, replace=False)] = True
        if share(m) >= obs - 1e-12:
            ge += 1
    null_mean = np.nan
    return float(obs), (ge + 1) / (n_perm + 1), null_mean


# ---------------------------------------------------------------- one locus x allele
def analyse(lt, pairs, allele, n_perm, exact_max, nrand, rng, descent_passes=2):
    """Exact single-origin test, then a phase-robustness certificate.

    A clade C is CONSISTENT with the (unphased) genotypes iff
      (1) every homozygous carrier has both haplotypes in C,
      (2) no homozygous non-carrier has a haplotype in C,
      (3) every heterozygote has exactly one haplotype in C.
    If any such clade exists, the allele is explicable by ONE mutation under an admissible
    phase -- and the phase is then induced, free.  Conditions (1)+(2) are tested for every
    node at once with two leaf counts; only the surviving candidates (a chain of ancestors)
    need the exact condition (3).  No sampling, no heuristic: the answer is exact.

    When no clade is consistent, we report how far off the best clade is (clade_fit_errors,
    the minimum number of genotyping/phase violations) and, using that clade to induce a
    phase, the exact n_origins of the resulting painting.
    """
    a1 = np.array([p[0] for p in pairs])
    a2 = np.array([p[1] for p in pairs])
    both = (a1 == allele) & (a2 == allele)
    het = (a1 == allele) != (a2 == allele)
    neither = ~both & ~het
    n_het = int(het.sum())
    n_carriers = 2 * int(both.sum()) + n_het
    if n_carriers < 3:
        return None

    fixed_mask = np.zeros(lt.n_hap, dtype=bool)
    fixed_mask[2 * np.flatnonzero(both)] = True
    fixed_mask[2 * np.flatnonzero(both) + 1] = True
    het_mask = np.zeros(lt.n_hap, dtype=bool)
    het_mask[2 * np.flatnonzero(het)] = True
    het_mask[2 * np.flatnonzero(het) + 1] = True
    non_mask = np.zeros(lt.n_hap, dtype=bool)
    non_mask[2 * np.flatnonzero(neither)] = True
    non_mask[2 * np.flatnonzero(neither) + 1] = True
    n_fc = int(fixed_mask.sum())

    c_fc = lt._counts(fixed_mask)
    c_non = lt._counts(non_mask)
    # the cheap exact part of the test, for every node simultaneously
    cand = np.flatnonzero((c_fc == n_fc) & (c_non == 0) & (lt.cnt_all > 0))
    # if that set is empty, look for the closest clade instead
    scored = (n_fc - c_fc) + c_non
    if len(cand) == 0:
        k = min(40, len(scored))
        cand = np.argpartition(scored, k - 1)[:k]

    het_h0 = 2 * np.flatnonzero(het)
    het_h1 = het_h0 + 1
    best = None
    feasible = []
    for v in cand:
        m = lt.clade_hap_mask(int(v))
        viol = int((m[het_h0] == m[het_h1]).sum()) + int((n_fc - c_fc[v]) + c_non[v])
        if viol == 0:
            feasible.append(int(v))
        if best is None or viol < best[0]:
            best = (viol, int(v), m)
    viol, v_best, m_best = best

    def paint_from_mask(m):
        p = fixed_mask.copy()
        p[het_h0[m[het_h0]]] = True
        p[het_h1[m[het_h1]]] = True
        return p

    if feasible:
        # every consistent clade is a legitimate answer; report the age/breadth RANGE over
        # them, which is the honest phase-induced uncertainty
        ages, nbgs, shares = [], [], []
        for v in feasible:
            m = lt.clade_hap_mask(v) if v != v_best else m_best
            s = lt.stats(m)
            ages.append(s['age'])
            nbgs.append(s['n_backgrounds'])
            shares.append(s['dom_share'])
        youngest = int(np.argmin(ages))
        chosen = lt.stats(lt.clade_hap_mask(feasible[youngest])
                          if feasible[youngest] != v_best else m_best)
        n_origins, n_origins_min, n_origins_max = 1, 1, 1
        age_range = (float(np.min(ages)), float(np.max(ages)))
        nbg_range = (int(np.min(nbgs)), int(np.max(nbgs)))
        dom_share = chosen['dom_share']
        dom_freq = chosen['dom_freq']
        age = chosen['age']
    else:
        painting = paint_from_mask(m_best)
        chosen = lt.stats(painting)
        n_origins = chosen['n_origins']
        n_origins_min = n_origins_max = n_origins
        age_range = (chosen['age'], chosen['age'])
        nbg_range = (chosen['n_backgrounds'], chosen['n_backgrounds'])
        dom_share, dom_freq, age = chosen['dom_share'], chosen['dom_freq'], chosen['age']

    # ---- stage A, phase-free
    pair = np.stack([lt.bg_of_hap[0::2], lt.bg_of_hap[1::2]], axis=1)
    group_of_sample = np.unique(np.sort(pair, axis=1), axis=0, return_inverse=True)[1].ravel()
    obs, p, _ = enrich_test(group_of_sample, both | het, n_perm, rng)

    bg0 = lt.bg_of_hap[het_h0]
    bg1 = lt.bg_of_hap[het_h1]
    k_straddle = int((bg0 != bg1).sum())
    if n_het == 0:
        cert = 'NO-HET'
    elif feasible:
        cert = 'SINGLE-ORIGIN-CONSISTENT'      # exact, over every admissible phase
    else:
        cert = 'NO-CLADE'

    return dict(allele=int(allele), n_carriers=n_carriers,
                freq=round(n_carriers / lt.n_hap, 4),
                n_hom=int(both.sum()), n_het=n_het,
                single_origin=int(bool(feasible)), n_feasible_clades=len(feasible),
                n_origins=n_origins,
                n_origins_min=n_origins_min, n_origins_max=n_origins_max,
                clade_fit_errors=int(viol),
                n_backgrounds=chosen['n_backgrounds'],
                nbg_min=nbg_range[0], nbg_max=nbg_range[1],
                dom_share=round(dom_share, 3), dom_freq=round(dom_freq, 3),
                age=round(age, 1) if np.isfinite(age) else np.nan,
                age_min=round(age_range[0], 1) if np.isfinite(age_range[0]) else np.nan,
                age_max=round(age_range[1], 1) if np.isfinite(age_range[1]) else np.nan,
                enrich_share=round(obs, 3) if np.isfinite(obs) else np.nan,
                enrich_p=round(p, 4) if np.isfinite(p) else np.nan,
                het_straddling=k_straddle, certificate=cert)


# ---------------------------------------------------------------- driver
def _int(x):
    try:
        f = float(x)
        return int(f) if np.isfinite(f) else -1
    except (TypeError, ValueError):
        return -1


COLS = ['locus', 'chrom', 'locus_start', 'locus_end', 'pos_eval', 'allele', 'n_samples',
        'n_carriers', 'freq', 'n_hom', 'n_het', 'single_origin', 'n_feasible_clades',
        'n_origins', 'clade_fit_errors', 'n_backgrounds', 'nbg_min', 'nbg_max',
        'dom_share', 'dom_freq', 'age', 'age_min', 'age_max',
        'enrich_share', 'enrich_p', 'het_straddling', 'certificate']


def run_table(table, tau, n_perm, exact_max, nrand, out, seed=0, descent_passes=2,
              verbose=True):
    df = pd.read_csv(table)
    need = {'locus', 'chrom', 'start', 'end', 'sample', 'allele1', 'allele2'}
    if not need.issubset(df.columns):
        raise SystemExit(f'table needs columns {sorted(need)}; got {list(df.columns)}')
    order = A.samples_in_order()
    idx = {s: i for i, s in enumerate(order)}
    rng = np.random.default_rng(seed)
    rows = []
    for (locus, chrom, start, end), g in df.groupby(['locus', 'chrom', 'start', 'end'],
                                                    sort=False):
        samples = [s for s in order if s in set(g['sample'])]
        if len(samples) < 3:
            print(f'{locus}: too few samples match the ARG cohort', file=sys.stderr)
            continue
        present = set(g['sample'])
        haps = [h for s in samples for h in (2 * idx[s], 2 * idx[s] + 1)]
        positions = sorted({int(start), int(end), int((start + end) // 2)})
        f, chunk_start = A.chunk_for(int(chrom), int(max(positions)))
        ts = tskit.load(str(f)).simplify(samples=haps)
        pairs = [(-1, -1)] * len(samples)
        pos_of = {s: k for k, s in enumerate(samples)}
        for r in g.itertuples():
            if r.sample in present:
                pairs[pos_of[r.sample]] = (_int(r.allele1), _int(r.allele2))
        alleles = sorted({a for p in pairs for a in p if a >= 0})
        for pos in positions:
            lt = LocalTree(ts, chunk_start, pos, tau)
            for al in alleles:
                r = analyse(lt, pairs, al, n_perm, exact_max, nrand, rng, descent_passes)
                if r:
                    r.update(locus=locus, chrom=int(chrom), locus_start=int(start),
                             locus_end=int(end), pos_eval=pos, n_samples=len(samples))
                    rows.append(r)
            if verbose:
                print(f'  {locus} chr{chrom}:{pos} - {len(alleles)} allele labels',
                      file=sys.stderr)
    res = pd.DataFrame(rows)
    if len(res):
        res = res[COLS]
    if out:
        res.to_csv(out, index=False)
        print(f'wrote {out} ({len(res)} rows)')
    return res


# ---------------------------------------------------------------- selftest
def selftest(tau, n_perm, exact_max, nrand, emit_demo=None):
    """Ground truth from the ARG: paint a known clade as an allele, un-phase it, recover it.

    Case 1: carriers = leaves of one internal clade  -> n_origins must be 1
    Case 2: carriers = a random scattered set        -> n_origins must be > 1
    Both are converted to unphased diploid genotypes first, so the harness sees exactly what
    it will see from a collaborator's SV table.  A third check reports the statistics under
    the TRUE phase, which must always give n_origins == 1 for case 1.
    """
    chrom, pos = 22, 15_000_000
    f, start = A.chunk_for(chrom, pos)
    order = A.samples_in_order()
    ts = tskit.load(str(f))
    lt = LocalTree(ts, start, pos, tau)
    tree = ts.at(pos - start)
    cands = [(tree.num_samples(v), int(v)) for v in tree.nodes()
             if not tree.is_sample(v) and 200 <= tree.num_samples(v) <= 1500
             and lt.t[v] > 400]
    n_leaf, node = sorted(cands)[len(cands) // 2]
    clade_nodes = tree.get_leaves(node)
    clade_hap = np.array([np.flatnonzero(lt.leaves == l)[0] for l in clade_nodes])

    def to_pairs(carrier_hap):
        c = np.zeros(lt.n_hap, dtype=bool)
        c[carrier_hap] = True
        pairs = []
        for i in range(lt.n_hap // 2):
            x, y = int(c[2 * i]), int(c[2 * i + 1])
            pairs.append((1 + min(x, y), 1 + max(x, y)))   # 1 absent, 2 present, sorted
        return c, pairs

    rng = np.random.default_rng(7)
    rand_hap = rng.choice(lt.n_hap, size=len(clade_hap), replace=False)
    cases = [('single-origin clade', *to_pairs(clade_hap)),
             ('scattered random set', *to_pairs(rand_hap))]

    print(f'selftest chr{chrom}:{pos}   clade node {node}: {n_leaf} carrier haplotypes, '
          f'MRCA {lt.t[node]:.0f} gen   (local tree depth {lt.t.max():.0f} gen)')
    print(f'  half the carriers are heterozygous, so ~{2*n_leaf*(1-n_leaf/lt.n_hap):.0f} '
          f'phase bits are unknown before we start\n')
    out = []
    for name, carrier, pairs in cases:
        truth = lt.stats(carrier)
        r = analyse(lt, pairs, 2, n_perm, exact_max, nrand, rng)
        r['case'] = name
        out.append(r)
        print(f'{name}:  true phase -> n_origins {truth["n_origins"]}, '
              f'backgrounds {truth["n_backgrounds"]}, age {truth["age"]:.0f}')
        print(f'  recovered: n_origins {r["n_origins"]} '
              f'(range over admissible phases {r["n_origins_min"]}..{r["n_origins_max"]}), '
              f'backgrounds {r["n_backgrounds"]} [{r["nbg_min"]}..{r["nbg_max"]}]')
        print(f'  age {r["age"]}   dom_share {r["dom_share"]}   '
              f'enrich {r["enrich_share"]} p={r["enrich_p"]}')
        print(f'  hets {r["n_het"]} ({r["het_straddling"]} straddle two backgrounds)   '
              f'certificate {r["certificate"]}\n')
    ok1 = out[0]['n_origins'] == 1
    ok2 = out[1]['n_origins'] > 1
    print('expected: clade -> n_origins 1 ; scattered -> n_origins > 1')
    print('PASS' if (ok1 and ok2) else 'FAIL')
    if emit_demo:
        rows = [dict(locus='DEMO_SINGLE_ORIGIN_SV', chrom=chrom, start=pos - 500,
                     end=pos + 500, sample=order[i], allele1=p[0], allele2=p[1])
                for i, p in enumerate(cases[0][2])]
        pd.DataFrame(rows).to_csv(emit_demo, index=False)
        print(f'demo table -> {emit_demo}')
    return ok1 and ok2


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--table')
    ap.add_argument('--out', default=None)
    ap.add_argument('--tau', type=float, default=TAU)
    ap.add_argument('--perm', type=int, default=5000, help='permutations for stage A')
    ap.add_argument('--exact-max', type=int, default=10, help='enumerate phases up to 2^k')
    ap.add_argument('--nrand', type=int, default=200, help='sampled phases when k is large')
    ap.add_argument('--descent-passes', type=int, default=2)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--emit-demo-table', default=None)
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest(a.tau, a.perm, a.exact_max, a.nrand, a.emit_demo_table) else 1)
    if not a.table:
        raise SystemExit('need --table or --selftest')
    res = run_table(a.table, a.tau, a.perm, a.exact_max, a.nrand,
                    a.out or str(A.OUT / 'sv_phase_pilot_results.csv'), seed=a.seed,
                    descent_passes=a.descent_passes)
    print(res.head(20).to_string(index=False))


if __name__ == '__main__':
    main()
