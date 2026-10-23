#!/usr/bin/env python3
"""Erik's tau-free breadth readout: per-class dom_share(tau) profile + crossover tau*.

From the existing multi-tau runs (global_tau200/500/1000 on either substrate):
for each structural class, interpolate the tau at which the dominant background
comes to hold >= 50% of carriers (TH). Confinement is genealogically meaningful
only when tau* < the class's clade age; at tau >= age the background contains
all carriers by construction (docs/loci-analysis-note.md, 2026-10-20 correction).

Usage: tau_star_profiles.py [--gw-dir data/derived/old_dups/genome_wide_erikarg]
"""
import argparse
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAUS = [200.0, 500.0, 1000.0]
TH = 0.5


def tau_star(ds):
    if ds[200.0] >= TH:
        return 200.0, True
    for a, b in zip(TAUS, TAUS[1:]):
        if ds[b] >= TH:
            f = (TH - ds[a]) / (ds[b] - ds[a])
            return math.exp(math.log(a) + f * (math.log(b) - math.log(a))), True
    return 1000.0, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gw-dir', default=str(ROOT / 'data/derived/old_dups/genome_wide_erikarg'))
    a = ap.parse_args()
    cfgs = {}
    for t in TAUS:
        d = {}
        for r in csv.DictReader(open(f'{a.gw_dir}/global_tau{int(t)}/summary.csv')):
            d[(r['chrom'], r['pos'], r['allele'])] = r
        cfgs[t] = d
    common = set(cfgs[200.0]) & set(cfgs[500.0]) & set(cfgs[1000.0])
    rows = []
    for k in common:
        ds = {t: float(cfgs[t][k]['dom_share']) for t in TAUS}
        age = float(cfgs[200.0][k]['age'])
        ts, confined = tau_star(ds)
        vacuous = ts >= age
        rows.append(dict(chrom=k[0], pos=k[1], allele=k[2],
                         n=cfgs[200.0][k]['n'], freq=cfgs[200.0][k]['freq'], age=cfgs[200.0][k]['age'],
                         dom_share_200=ds[200.0], dom_share_500=ds[500.0], dom_share_1000=ds[1000.0],
                         tau_star=round(ts, 1) if confined else '',
                         status=('CONFINED<200' if confined and ts <= 200.0 else
                                 'confined-in-profile' if confined and not vacuous else
                                 'by-construction' if confined and vacuous else
                                 'not-confined<=1000')))
    out = f'{a.gw_dir}/tau_star_profiles.csv'
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    from collections import Counter
    print(f'{len(rows)} classes; status:', dict(Counter(r['status'] for r in rows)))
    print('wrote', out)


if __name__ == '__main__':
    main()
