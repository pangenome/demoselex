#!/usr/bin/env python3
"""Temporal trajectories: structural-class frequencies through time in the ancients.

The static analyzer gives today's snapshot; this gives the selection test -
whether a class ROSE through time, faster than drift allows, which also catches
selection on standing (old) variation that age-based filtering cannot see.

Per (locus, allele class):
  freq_by_bin   frequency in each time bin of the dated, called ancients
                (equal-count quantile bins) + the modern cohort as the bin at 0
  rise          freq_modern - freq_oldest_bin, and the largest adjacent-bin jump
  z_drift       rise in units of the diffusion drift null,
                Var(dp) = p(1-p)*dt_gen/(2Ne), for Ne in {5000, 10000}
  p_perm        date-permutation null: fraction of permutations of dates across
                the ancient samples producing a rise >= observed
  shape         punctuated (one adjacent jump > 60% of total) vs gradual
  verdict       RISING if p_perm <= 0.05 and z_drift >= 3 in BOTH Ne scenarios;
                falling/flat/oscillating otherwise

Ancients are dated from the date table (date_bp; moderns sit at 0). Samples
without dates or without calls at a locus are dropped per-locus. The ancient
cohort is population-mixed: verdicts must be read with the turnover caveat;
per-population stratification is a refinement.

Usage:
  temporal_traj.py [--table alleles_clusters.tsv] [--dates dates.tsv]
                   [--min-class 20] [--bin-samples 120] [--perms 200]
                   [--gen-years 28] [--loci chr,pos[,...]] [--out trajectories.csv]
  temporal_traj.py --selftest          # synthetic rise, flat and decline cases
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEF_TABLE = ROOT / 'data/derived/old_dups/alleles_clusters.tsv'
DEF_DATES = ROOT / 'data/derived/old_dups/dates.tsv'
DEF_OUT = ROOT / 'data/derived/old_dups/genome_wide/trajectories.csv'


def trajectory(samples, carriers, dates_bp, bin_size, gen_years):
    """samples: ancient sample ids with dates (sorted by date ascending = oldest first).
    carriers: (n_samples x n_classes) counts 0/1/2 per sample.
    dates_bp: matching dates (>0 ancient). Returns per-class per-bin freqs + stats."""
    n = len(samples)
    n_bins = max(2, n // bin_size)
    # equal-count bins over the date-sorted samples
    bin_id = np.concatenate([np.full(bin_size, b) for b in range(n_bins)])[:n]
    if n % bin_size:
        # pad the last bin with the remainder
        bin_id = np.concatenate([bin_id, np.full(n - len(bin_id), n_bins - 1)])
    called = carriers.sum(axis=1)                       # called haps per sample
    W = np.zeros((n_bins, n))
    W[bin_id, np.arange(n)] = 1.0
    bin_called = W @ called
    bin_carriers = W @ carriers                         # (n_bins x n_classes)
    bin_date = W @ dates_bp / np.maximum(W @ np.ones(n), 1e-9)
    with np.errstate(invalid='ignore', divide='ignore'):
        freq = bin_carriers / np.maximum(bin_called[:, None], 1e-9)
    return bin_id, bin_called, bin_date, freq


def rise_stats(bin_called, bin_date, freq, modern_freq, modern_called, gen_years):
    """Per class: oldest-bin freq, modern freq, rise, biggest adjacent jump, shape."""
    out = []
    n_bins, n_classes = freq.shape
    for k in range(n_classes):
        f = freq[:, k]
        f0 = f[0]
        fm = modern_freq[k] if np.isfinite(modern_freq[k]) else np.nan
        rise = fm - f0 if np.isfinite(fm) else np.nan
        adj = [f[i + 1] - f[i] for i in range(n_bins - 1)
               if np.isfinite(f[i]) and np.isfinite(f[i + 1])]
        # include the modern step (last ancient bin -> modern)
        if np.isfinite(f[-1]) and np.isfinite(fm):
            adj.append(fm - f[-1])
        jump = max(adj) if adj else np.nan
        total = (fm - f0) if np.isfinite(fm) else np.nan
        shape = ('punctuated' if (np.isfinite(jump) and np.isfinite(total)
                 and total > 0 and jump >= 0.6 * total) else 'gradual')
        # drift null over [oldest bin, modern]
        dt_years = float(bin_date[0])                       # oldest bin mean date
        dt_gen = dt_years / gen_years
        p0 = f0 if f0 > 0 else 1.0 / max(modern_called, 1)
        sd = np.sqrt(p0 * (1 - p0) * dt_gen / (2 * np.array([5000.0, 10000.0])))
        z = (total / sd) if (np.isfinite(total) and sd[0] > 0) else (np.nan, np.nan)
        out.append(dict(freq_first=round(float(f0), 4) if np.isfinite(f0) else np.nan,
                        freq_modern=round(float(fm), 4) if np.isfinite(fm) else np.nan,
                        rise=round(float(total), 4) if np.isfinite(total) else np.nan,
                        max_jump=round(float(jump), 4) if np.isfinite(jump) else np.nan,
                        shape=shape, z_ne5000=round(float(z[0]), 2) if sd[0] > 0 else np.nan,
                        z_ne10000=round(float(z[1]), 2) if sd[0] > 0 else np.nan))
    return out


def perm_null(samples, carriers, dates_bp, bin_size, modern_freq, n_perm, rng):
    """Empirical null distribution of the total rise: dates are exchangeable, so
    permute WHICH samples fall in which bin (rows of the carrier matrix)."""
    n, n_classes = carriers.shape
    rises = np.zeros((n_perm, n_classes))
    for p in range(n_perm):
        Cp = carriers[rng.permutation(n)]
        _, bin_called, _, freq = trajectory(samples, Cp, dates_bp, bin_size, 28)
        with np.errstate(invalid='ignore'):
            r = modern_freq - freq[0]
        rises[p] = np.where(np.isfinite(r), r, np.nan)
    return rises


def run(table_path, dates_path, args):
    df = pd.read_csv(table_path, sep=None, engine='python')
    for col in ('allele1', 'allele2'):
        df[col] = pd.to_numeric(df[col], errors='coerce')
    if getattr(args, 'modern_pop', None):
        popf = ROOT / 'data/derived/g1k/samples_population.txt'
        keep = set()
        with open(popf) as f:
            next(f)
            for line in f:
                p = line.split()
                if len(p) >= 7 and p[6] == args.modern_pop:
                    keep.add(p[1])
        df = pd.concat([df[~df['sample'].str.match(r'^(HG|NA)\d+$')],
                        df[df['sample'].isin(keep)]])
    dates = pd.read_csv(dates_path, sep=None, engine='python')
    dmap = {}
    for r in dates.itertuples():
        if not pd.isna(r.date_bp):
            dmap[r.sample] = float(r.date_bp)
    loci_sel = None
    if args.loci:
        loci_sel = {tuple(int(x) for x in s.split(':')) for s in args.loci.split(',')}
    rng = np.random.default_rng(args.seed)
    out_rows = []
    for (chrom, pos), g in df.groupby(['chrom', 'pos']):
        if loci_sel and (int(chrom), int(pos)) not in loci_sel:
            continue
        mods = g[g['sample'].str.match(r'^(HG|NA)\d+$')]
        ancs = g[~g['sample'].str.match(r'^(HG|NA)\d+$')]
        # dated + called ancients
        rows = []
        for r in ancs.itertuples():
            d = dmap.get(r.sample)
            if d is None or d <= 0 or (pd.isna(r.allele1) and pd.isna(r.allele2)):
                continue
            rows.append((r.sample, d, r.allele1, r.allele2))
        if getattr(args, 'anc_group_file', None):
            allowed = set()
            with open(args.anc_group_file) as f:
                for line in f:
                    p = line.rstrip('\n').split('\t')
                    if len(p) >= 2 and (not getattr(args, 'group', None)
                                        or p[1] == args.group):
                        allowed.add(p[0])
            rows = [t for t in rows if t[0] in allowed]
        if len(rows) < args.bin_samples * 2:
            continue
        rows.sort(key=lambda t: -t[1])                     # oldest first
        samples = [t[0] for t in rows]
        d_arr = np.array([t[1] for t in rows])
        # class inventory: everything carried by >= min_class haps (ancients+moderns)
        counts = {}
        for _, _, a1, a2 in rows:
            if not pd.isna(a1):
                counts[int(a1)] = counts.get(int(a1), 0) + 1
            if not pd.isna(a2):
                counts[int(a2)] = counts.get(int(a2), 0) + 1
        for r in mods.itertuples():
            for a in (r.allele1, r.allele2):
                if not pd.isna(a):
                    counts[int(a)] = counts.get(int(a), 0) + 1
        classes = sorted(k for k, v in counts.items() if v >= args.min_class)
        if not classes:
            continue
        cidx = {c: i for i, c in enumerate(classes)}
        C = np.zeros((len(rows), len(classes)))
        for i, (_, _, a1, a2) in enumerate(rows):
            if not pd.isna(a1) and int(a1) in cidx:
                C[i, cidx[int(a1)]] += 1
            if not pd.isna(a2) and int(a2) in cidx:
                C[i, cidx[int(a2)]] += 1
        bin_id, bin_called, bin_date, freq = trajectory(samples, C, d_arr,
                                                         args.bin_samples, args.gen_years)
        # modern frequencies
        mC = np.zeros(len(classes))
        m_tot = 0
        for r in mods.itertuples():
            for a in (r.allele1, r.allele2):
                if not pd.isna(a) and int(a) in cidx:
                    mC[cidx[int(a)]] += 1
                    m_tot += 1
        modern_freq = mC / max(m_tot, 1)
        stats = rise_stats(bin_called, bin_date, freq, modern_freq, m_tot, args.gen_years)
        # permutation null on the total rise
        if args.perms:
            pr = perm_null(samples, C, d_arr, args.bin_samples, modern_freq,
                            args.perms, rng)
            obs = np.array([s['rise'] for s in stats], dtype=float)
            with np.errstate(invalid='ignore'):
                ge = np.nansum(pr >= obs[None, :], axis=0)
                nvalid = np.sum(np.isfinite(pr), axis=0)
            pvals = (ge + 1) / (nvalid + 1)
        else:
            pvals = np.full(len(classes), np.nan)
        for i, c in enumerate(classes):
            st = stats[i]
            verdict = 'RISING' if (st['rise'] and st['rise'] > 0
                      and pvals[i] <= 0.05
                      and min(st['z_ne5000'], st['z_ne10000']) >= 3.0) else \
                      ('falling' if (st['rise'] and st['rise'] < -0.05) else 'flat')
            bins_str = ';'.join(f'{int(bin_date[b])}:{round(float(freq[b, i]), 4)}'
                                for b in range(len(bin_date)))
            out_rows.append(dict(chrom=int(chrom), pos=int(pos), allele=c,
                                 group=getattr(args, 'group', 'ALL') or 'ALL',
                                 n_class_haps=counts[c], n_anc=len(rows),
                                 n_bins=len(bin_date),
                                 freqs_by_bin=bins_str, **st,
                                 p_perm=round(float(pvals[i]), 4),
                                 verdict=verdict))
        print(f'chr{chrom}:{pos}: {len(rows)} dated called ancients, '
              f'{len(classes)} classes, {len(out_rows)} rows so far', flush=True)
    res = pd.DataFrame(out_rows)
    res.to_csv(args.out, index=False)
    print(f'wrote {args.out}: {len(res)} trajectories')
    if len(res):
        n_rise = (res['verdict'] == 'RISING').sum()
        print(f'RISING (p_perm<=0.05, z>=3 at Ne 5k and 10k): {n_rise}')
    return res


def selftest():
    """Synthetic check: one rising class, one flat, one falling; permutation p
    must be small for the rising one and large for the flat one."""
    import tempfile
    d = Path(tempfile.mkdtemp())
    rows = ['chrom\tpos\tsample\tallele1\tallele2']
    # 120 ancients at 9000 BP: class 1 at 10%, class 2 at 50%, class 3 at 30%
    # 120 ancients at 1000 BP: class 1 at 70%, class 2 at 50%, class 3 at 5%
    rng = np.random.default_rng(3)
    def emit(tag, n, probs):
        p = np.array(probs, dtype=float); p /= p.sum()
        for i in range(n):
            a1, a2 = rng.choice([1, 2, 3], size=2, p=p)
            rows.append(f'7\t1000\t{tag}{i}\t{a1}\t{a2}')
    emit('OLD', 120, [0.10, 0.20, 0.70])
    emit('NEW', 120, [0.70, 0.20, 0.10])
    for i in range(200):                                        # moderns same as NEW
        a1, a2 = rng.choice([1, 2, 3], size=2, p=[0.70, 0.20, 0.10])
        rows.append(f'7\t1000\tHG{i:05d}\t{a1}\t{a2}')
    (d / 't.tsv').write_text('\n'.join(rows) + '\n')
    drows = ['sample\tdate_bp\tdate_sd\tsource\tculture']
    for i in range(120):
        drows.append(f'OLD{i}\t9000\t100\t14C\tx')
    for i in range(120):
        drows.append(f'NEW{i}\t1000\t100\t14C\tx')
    for i in range(200):
        drows.append(f'HG{i:05d}\t0\t0\tmodern\tmodern')
    (d / 'd.tsv').write_text('\n'.join(drows) + '\n')
    args = argparse.Namespace(table=str(d / 't.tsv'), dates=str(d / 'd.tsv'),
                               out=str(d / 'o.csv'), min_class=20, bin_samples=120,
                               perms=100, gen_years=28, loci=None, seed=1)
    res = run(args.table, args.dates, args)
    r = {int(x['allele']): x for x in res.to_dict('records')}
    ok = True
    ok &= r[1]['verdict'] == 'RISING' and r[1]['rise'] > 0.4
    ok &= r[2]['verdict'] == 'flat' and abs(r[2]['rise']) < 0.15
    ok &= r[3]['verdict'] == 'falling' and r[3]['rise'] < -0.15
    ok &= r[1]['p_perm'] <= 0.1 and r[2]['p_perm'] > 0.2
    print('\nSELFTEST', 'PASSED' if ok else 'FAILED',
          {k: (r[k]['verdict'], r[k]['rise'], r[k]['p_perm']) for k in (1, 2, 3)})
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--table', default=str(DEF_TABLE))
    ap.add_argument('--dates', default=str(DEF_DATES))
    ap.add_argument('--out', default=str(DEF_OUT))
    ap.add_argument('--min-class', type=int, default=20)
    ap.add_argument('--bin-samples', type=int, default=120)
    ap.add_argument('--perms', type=int, default=200)
    ap.add_argument('--gen-years', type=int, default=28)
    ap.add_argument('--loci', default=None)
    ap.add_argument('--modern-pop', default=None,
                    help='restrict moderns to one superpopulation (e.g. EUR) so rises are '
                         'like-for-like against the mostly-European ancients')
    ap.add_argument('--anc-group-file', default=None,
                    help='tsv sample<TAB>group; restrict ancients to one group (--group)')
    ap.add_argument('--group', default=None)
    ap.add_argument('--seed', type=int, default=11)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    run(a.table, a.dates, a)


if __name__ == '__main__':
    main()
