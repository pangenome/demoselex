#!/usr/bin/env python3
"""Genome-wide run of the delivered sv_phase_pilot (exact single-origin test, Q2).

Builds pilot-format tables from alleles_clusters.tsv (moderns only: the ARG has
no ancients), runs erik's sv_phase_pilot.py verbatim per chromosome on the
delivered full-genome ARG (G1K_ARG_DIR), merges, and cross-tabs against the
tau* profiles.

Known edges handled:
- loci straddling two chunk FILES are degraded to midpoint-only
- ancients (not in the ARG cohort) are dropped from the tables
- TPTE2 (chr13:19.3-19.6 Mb) has no tree coverage on either build: the chunk
  named chr13_10000000_20000000 is actually trimmed to 13,999,701 on both;
  the pilot crashes loudly there (the static analyzer skipped it silently).
  It is excluded and reported.

Usage: see scripts/run_sv_pilot_all.sh for the parallel run, then
       sv_pilot_scan.py --tables <dir> --merge <out.csv>  # merge + cross-tab
"""
import argparse
import csv
import glob
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BED = Path('/moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed')
GWE = ROOT / 'data/derived/old_dups/genome_wide_erikarg'
GW = ROOT / 'data/derived/old_dups/genome_wide'


def build_tables(tables_dir, argdir_glob, table):
    chunks = defaultdict(list)
    for f in glob.glob(argdir_glob):
        m = re.match(r'chr(\d+)_(\d+)_(\d+)\.dated\.trees', f.split('/')[-1])
        if m:
            chunks[int(m.group(1))].append((int(m.group(2)), int(m.group(3))))

    def chunk_of(c, p):
        for s, e in sorted(chunks.get(c, [])):
            if s <= p < e:
                return s
        return None

    bed, seen = {}, set()
    for line in open(BED):
        f = line.split()
        if len(f) >= 4:
            c = f[0].split('#')[-1].replace('chr', '')
            bed[(int(c), (int(f[1]) + int(f[2])) // 2)] = f[3]
    pop = {l.split()[1] for l in open('/moosefs/erikg/demoselex-arg/samples_population.txt')
           if len(l.split()) >= 2 and l.split()[0] != 'FamilyID'}
    rows_by_chrom = defaultdict(list)
    locus_names, cross = {}, 0
    for r in csv.DictReader(open(table), delimiter='\t'):
        if r['sample'] not in pop:
            continue
        c, pos = int(r['chrom']), int(r['pos'])
        s, e = int(r['start']), int(r['end'])
        gene = bed.get((c, pos), f"chr{c}:{pos}")
        if gene not in locus_names:
            if chunk_of(c, s) != chunk_of(c, e) or chunk_of(c, s) is None or chunk_of(c, e) is None:
                cross += 1
                s = e = pos
            if gene in seen:
                gene = f"{gene}@{pos}"
            seen.add(gene)
            locus_names[gene] = (c, s, e)
        c2, s2, e2 = locus_names[gene]
        rows_by_chrom[c2].append((gene, c2, s2, e2, r['sample'], r['allele1'], r['allele2']))
    outdir = Path(tables_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    for c, rows in sorted(rows_by_chrom.items()):
        with open(outdir / f'chr{c}.csv', 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['locus', 'chrom', 'start', 'end', 'sample', 'allele1', 'allele2'])
            w.writerows(rows)
    return cross, len(locus_names)


def merge(tables_dir, glob_pat, out):
    rows = []
    for f in sorted(glob.glob(glob_pat)):
        with open(f) as fh:
            rd = csv.DictReader(fh)
            rows.extend(rd)
    if not rows:
        raise SystemExit('no per-chrom CSVs found')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return rows


def crosstab(rows, out_dir):
    bed = {}
    for line in open(BED):
        f = line.split()
        if len(f) >= 4:
            c = f[0].split('#')[-1].replace('chr', '')
            bed[f[3]] = (c, str((int(f[1]) + int(f[2])) // 2))
    tstar = {}
    for r in csv.DictReader(open(Path(out_dir) / 'tau_star_profiles.csv')):
        tstar[(r['chrom'], r['pos'], r['allele'])] = r
    cls = defaultdict(list)
    for r in rows:
        cls[(r['locus'], r['allele'])].append(r)
    counts = defaultdict(int)
    for k, v in cls.items():
        cp = bed.get(k[0])
        t = tstar.get((cp[0], cp[1], k[1]), {}) if cp else {}
        confined = t.get('status', '') in ('CONFINED<200', 'confined-in-profile')
        counts[(confined, any(x['single_origin'] == '1' for x in v))] += 1
    print(f"classes: {len(cls)}")
    print(f"  tau* confined   & single-origin>=1pos: {counts[(True, True)]}"
          f" | no single-origin: {counts[(True, False)]}")
    print(f"  tau* distributed& single-origin>=1pos: {counts[(False, True)]}"
          f" | no single-origin: {counts[(False, False)]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--build', action='store_true', help='build per-chrom pilot tables')
    ap.add_argument('--table', default=str(ROOT / 'data/derived/old_dups/alleles_clusters.tsv'))
    ap.add_argument('--argdir', default='/moosefs/erikg/demoselex-arg/arg')
    ap.add_argument('--tables', default=str(ROOT / '../pilot_run/tables'))
    ap.add_argument('--merge', help='merge per-chrom CSVs into one file')
    ap.add_argument('--crosstab', help='pilot_genome.csv for the tau* cross-tab')
    a = ap.parse_args()
    if a.build:
        cross, n = build_tables(a.tables, f"{a.argdir}/chr*_*.dated.trees", a.table)
        print(f"built tables for {n} loci ({cross} degraded to midpoint, ancients dropped)")
    if a.merge:
        rows = merge(a.tables, str(Path(a.merge).with_name('chr*.csv').name), a.merge)
        print(f"merged {len(rows)} rows -> {a.merge}")
    if a.crosstab:
        rows = list(csv.DictReader(open(a.crosstab)))
        crosstab(rows, str(GWE))


if __name__ == '__main__':
    main()
