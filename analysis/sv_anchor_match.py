#!/usr/bin/env python3
"""Anchor matching: let an assembled SV haplotype tell us which ARG haplotype carries it.

The simple version of phasing, and deliberately so.  An SV call ships with TWO region
haplotypes (the reference path and the alternate path, usually from an HPRC assembly).
Take the variants the alternate haplotype implies in the FLANKS of the call, score every
1000G haplotype in the ARG by how many of those alleles it shares, and the closest ones
are the haplotypes that sit on the same ancestry background as the SV.  That is "stitching"
the call onto the panel by nearest local match -- no HMM, no 2^k phase space, no
Li-Stephens panel assumptions.

What it does and does not give you:

  * For the DONOR (the HPRC sample the alt path came from) the assembly IS one of its
    haplotypes, so the match is a real observation: it sets phase exactly.
  * For any other sample the match identifies the ANCESTRY BACKGROUND of the alt
    haplotype, not that sample's phase.  It says "an SV that looks like this lives on
    haplotypes resembling these", which is exactly what the enrichment test in
    sv_phase_pilot.py needs -- and the clade test there is the check on whether that
    background assignment is consistent with a single origin.

Why flanks and not the SV interior: the alternate path is structurally divergent inside
the call, so variants there do not align with the panel.  The flanking variants are shared
ancestry and are the informative part.

Usage:
  python3 sv_anchor_match.py --selftest
  # library use:
  #   scores, informative = anchor_scores(ts, start, end, flank, donor_hap)
"""
import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parent))
import g1k_allele_analyzer as A


def region_genotypes(ts, lo, hi):
    """Genotype matrix for the sites in [lo, hi) (chunk-relative), alleles as small ints."""
    pos = np.asarray(ts.sites_position)
    i0, i1 = int(np.searchsorted(pos, lo)), int(np.searchsorted(pos, hi))
    if i1 <= i0:
        return np.zeros((0, ts.num_samples), dtype=np.int8), np.zeros(0, dtype=np.int64)
    G = np.array([v.genotypes for v in itertools.islice(ts.variants(), i0, i1)],
                 dtype=np.int8)
    return G, pos[i0:i1]


def anchor_scores(ts, start, pos, flank, donor_hap, exclude_inner=True, inner=2000):
    """Score every haplotype by sharing with the donor's alleles in the flanking windows.

    flank       window width on each side
    donor_hap   positional haplotype index whose implied variant pattern is the anchor
    exclude_inner  drop sites within +/-inner of the call (structural divergence)
    Returns (scores, n_informative, top_haps) with scores in [0,1].
    """
    rel = pos - start
    Gs = []
    for lo, hi in ((rel - flank, rel - inner if exclude_inner else rel + flank),
                   (rel + inner if exclude_inner else rel - flank, rel + flank)):
        G, _ = region_genotypes(ts, lo, hi)
        if len(G):
            Gs.append(G)
    if not Gs:
        return np.zeros(ts.num_samples), 0, np.zeros(0, dtype=int)
    G = np.vstack(Gs)
    alt = G[:, donor_hap]
    inf = alt != 0                        # sites where the anchor differs from the reference
    if inf.sum() == 0:
        return np.zeros(ts.num_samples), 0, np.zeros(0, dtype=int)
    match = (G[inf] == alt[inf][:, None]).sum(axis=0)
    scores = match / inf.sum()
    order = np.argsort(-scores)
    return scores, int(inf.sum()), order


def selftest(flank=20000, inner=2000, tau=200.0):
    chrom, pos = 22, 15_000_000
    f, start = A.chunk_for(chrom, pos)
    ts = tskit.load(str(f))
    tree = ts.at(pos - start)

    # a real clade to act as the ancestry background that carries the SV
    cands = [(tree.num_samples(v), int(v)) for v in tree.nodes()
             if not tree.is_sample(v) and 200 <= tree.num_samples(v) <= 1500]
    node = sorted(cands)[len(cands) // 2][1]
    clade_hap = np.array([np.flatnonzero(ts.samples() == l)[0] for l in tree.get_leaves(node)])
    in_clade = np.zeros(ts.num_samples, dtype=bool)
    in_clade[clade_hap] = True

    print(f'anchor selftest chr{chrom}:{pos}   clade node {node}: {len(clade_hap)} haplotypes, '
          f'MRCA {ts.tables.nodes.time[node]:.0f} gen')

    donor = int(clade_hap[0])             # pretend this hap's assembly is the "alt path"
    scores, n_inf, order = anchor_scores(ts, start, pos, flank, donor, inner=inner)
    print(f'  {n_inf} informative flanking variant alleles within +/-{flank} bp '
          f'(excluding +/-{inner} bp)\n')
    close = scores >= 0.99
    print(f'  donor haplotype {donor}: score {scores[donor]:.3f} (rank '
          f'{int(np.flatnonzero(order == donor)[0]) + 1}/{len(scores)})')
    print(f'  haplotypes matching the anchor at >=0.99: {int(close.sum())}, of which '
          f'{int((close & in_clade).sum())} are in the true carrier clade')
    print(f'  mean score inside the clade {scores[in_clade].mean():.3f} vs outside '
          f'{scores[~in_clade].mean():.3f}')

    # does the anchor pick the right haplotype for heterozygotes?
    donor_out = int(np.flatnonzero(~in_clade)[0])
    s_out, _, _ = anchor_scores(ts, start, pos, flank, donor_out, inner=inner)
    ok = tot = 0
    het_samples = [i for i in range(ts.num_samples // 2)
                   if in_clade[2 * i] != in_clade[2 * i + 1]]
    for i in het_samples[:200]:
        s = scores if in_clade[2 * i] else s_out        # anchor from the carrying side
        pick = int(s[2 * i + 1] > s[2 * i])
        truth = int(in_clade[2 * i + 1])
        ok += pick == truth
        tot += 1
    print(f'\n  phase by anchor on {tot} heterozygous samples: {ok}/{tot} correct '
          f'({100 * ok / max(tot, 1):.1f}%)')

    # negative control: an anchor from a distant branch should not match the clade
    distant = int(order[-1])
    s_d, n_inf_d, _ = anchor_scores(ts, start, pos, flank, distant, inner=inner)
    print(f'  control, anchor from the least-related haplotype {distant}: '
          f'mean score inside the clade {s_d[in_clade].mean():.3f}')
    print(f'\n  read: a single 20 kb flanking anchor pins the carrier background to '
          f'{int(close.sum())} haplotypes out of {len(scores)}')
    return ok == tot


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--flank', type=int, default=20000)
    ap.add_argument('--inner', type=int, default=2000)
    ap.add_argument('--recomb-selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest(a.flank, a.inner) else 1)
    if a.recomb_selftest:
        sys.exit(0 if recomb_selftest(a.flank) else 1)
    raise SystemExit('nothing to do: use --selftest, or import anchor_scores()')




# ---------------------------------------------------------------- recombination-aware matching
def load_region(ts, lo, hi):
    """Sites in [lo, hi) (chunk-relative): positions, allele strings, genotype matrix."""
    pos = np.asarray(ts.sites_position)
    i0, i1 = int(np.searchsorted(pos, lo)), int(np.searchsorted(pos, hi))
    if i1 <= i0:
        return np.zeros(0, np.int64), [], np.zeros((0, ts.num_samples), np.int8)
    G, alleles = [], []
    for v in itertools.islice(ts.variants(), i0, i1):
        G.append(v.genotypes)
        alleles.append(v.alleles)
    return pos[i0:i1], alleles, np.array(G, dtype=np.int8)


def match_costs(G, alleles, alt_alleles):
    """0/1 mismatch matrix (kept sites x haplotypes) against an alt haplotype's BASE strings.

    Compares bases, not allele indexes: the ARG stores ancestral/derived orientation while a
    locally built VCF stores REF/ALT, so index comparison would silently scramble the match.
    Sites where the alt carries the ancestral/REF base are uninformative and dropped; sites
    where the alt carries a base absent from the panel are counted as novel.
    """
    keep, alt_idx, n_ref, n_novel = [], [], 0, 0
    for i, (al, a) in enumerate(zip(alleles, alt_alleles)):
        if a is None or i >= len(G):
            n_novel += 1
            continue
        if a not in al:
            n_novel += 1
            continue
        k = al.index(a)
        if k == 0:
            n_ref += 1
            continue
        keep.append(i)
        alt_idx.append(k)
    if not keep:
        return np.zeros((0, G.shape[1])), np.zeros(0, np.int64), dict(
            n_informative=0, n_ref=n_ref, n_novel=n_novel)
    keep = np.array(keep)
    alt_idx = np.array(alt_idx)
    M = (G[keep] != alt_idx[:, None]).astype(np.int8)
    return M, keep, dict(n_informative=len(keep), n_ref=n_ref, n_novel=n_novel)


def segment_mosaic(M, lam=None, top_k=20, n_windows=10, min_sites=8):
    """Split the alt haplotype into contiguous stretches, each copying one panel haplotype.

    Recombinant haplotypes are mosaics, so a single best match is wrong by construction.
    Two details make it work:

      * CANDIDATES are the union of the best matches inside each of `n_windows` windows.
        A panel haplotype that matches only the left half has a terrible whole-region score,
        so a global candidate list would throw away the very parent we are looking for.
      * BINARY SEGMENTATION on the mismatch matrix then finds the split saving the most
        mismatches, recursing while the saving exceeds `lam`; segments with fewer than
        `min_sites` informative sites are merged into their better-matching neighbour
        (a 3-site tract is not evidence of a breakpoint).

    Returns [(first_kept_site, last_kept_site, hap, mismatches, n_sites)].
    """
    if M.size == 0:
        return []
    lam = 1 if lam is None else lam
    n = len(M)
    edges = np.linspace(0, n, n_windows + 1).astype(int)
    cand = set()
    for a, b in zip(edges[:-1], edges[1:]):
        if b > a:
            cand.update(np.argsort(M[a:b].sum(axis=0))[:top_k].tolist())
    cand = np.array(sorted(cand))
    Mk = M[:, cand]
    out = []

    def rec(i, j):
        cs = np.cumsum(Mk[i:j], axis=0)
        whole = int(cs[-1].min())
        if j - i < 2:
            out.append([i, j - 1, int(cand[int(cs[-1].argmin())]), whole, j - i])
            return
        best_gain, best_s = 0.0, None
        for s_ in range(1, j - i):
            cost = int(cs[s_ - 1].min()) + int((cs[-1] - cs[s_ - 1]).min())
            if whole - cost > best_gain:
                best_gain, best_s = whole - cost, i + s_
        if best_s is None or best_gain <= lam:
            out.append([i, j - 1, int(cand[int(cs[-1].argmin())]), whole, j - i])
        else:
            rec(i, best_s)
            rec(best_s, j)

    rec(0, n)
    out.sort()
    # merge under-determined islands into the better-matching neighbour
    changed = True
    while changed and len(out) > 1:
        changed = False
        for k, seg in enumerate(out):
            if seg[4] >= min_sites:
                continue
            others = [o for o in out if o is not seg]
            nb = min(others, key=lambda o: min(abs(o[0] - seg[1]), abs(seg[0] - o[1])))
            a, b = min(seg[0], nb[0]), max(seg[1], nb[1])
            cs = np.cumsum(Mk[a:b + 1], axis=0)
            merged = [a, b, int(cand[int(cs[-1].argmin())]), int(cs[-1].min()), b - a + 1]
            out = [o for o in out if o is not seg and o is not nb] + [merged]
            out.sort()
            changed = True
            break
    return [tuple(o) for o in out]


def recomb_selftest(flank=50000):
    """Build a deliberate recombinant anchor from two real haplotypes, recover the mosaic."""
    chrom, pos = 22, 15_000_000
    f, start = A.chunk_for(chrom, pos)
    ts = tskit.load(str(f))
    rel = pos - start
    positions, alleles, G = load_region(ts, rel - flank, rel + flank)
    print(f'recombination selftest chr{chrom}:{pos}   {len(positions)} sites in +/-{flank} bp')

    tree = ts.at(rel)
    cands = sorted([(tree.num_samples(v), int(v)) for v in tree.nodes()
                    if not tree.is_sample(v) and 200 <= tree.num_samples(v) <= 1500])
    node = cands[len(cands) // 2][1]
    A_haps = [int(np.flatnonzero(ts.samples() == l)[0]) for l in tree.get_leaves(node)]
    hapA = A_haps[0]
    # hapB must be a genuinely different background, else the mosaic is undetectable
    insite = [i for i, al in enumerate(alleles) if len(al) > 1]
    if insite:
        sub = G[insite]
        dist = (sub != sub[:, [hapA]]).sum(axis=0)
        hapB = int(np.argmax(np.where(np.arange(ts.num_samples) == hapA, -1, dist)))
    else:
        hapB = next(h for h in range(ts.num_samples) if h not in set(A_haps))

    mid = len(positions) // 2
    brk = int(positions[mid]) + start
    alt = [alleles[i][int(G[i, hapA])] if i < mid else alleles[i][int(G[i, hapB])]
           for i in range(len(positions))]
    M, keep, info = match_costs(G, alleles, alt)
    print(f'  synthetic anchor = haplotype {hapA} left of {brk:,}, haplotype {hapB} right')
    print(f'  {info["n_informative"]} informative sites ({info["n_ref"]} ancestral/REF, '
          f'{info["n_novel"]} novel)')
    segs = segment_mosaic(M)
    for (i0, i1, h, mm, n) in segs:
        p0, p1 = positions[keep[i0]] + start, positions[keep[i1]] + start
        tag = '<- true A' if h == hapA else ('<- true B' if h == hapB else '')
        print(f'  segment {p0:>12,}-{p1:>12,}  {n:4d} sites  copies hap {h:5d}  '
              f'mismatch {mm}/{n}  {tag}')
    got = {h for (_, _, h, _, _) in segs}
    boundaries = [positions[keep[i0]] + start for (i0, _, _, _, _) in segs[1:]]
    ok = {hapA, hapB}.issubset(got)
    print(f'  breakpoints detected at: {boundaries if boundaries else "none"} '
          f'(true {brk:,})')
    print('PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    main()
