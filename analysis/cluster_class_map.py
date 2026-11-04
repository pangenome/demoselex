#!/usr/bin/env python3
"""Decoder table: MDS cluster -> class label -> carrier count, per locus.

Joins the three levels of the class construction:
  L0a  paths_<locus>.tsv       MDS/DBSCAN cluster of every catalog path (~466/locus)
  L0b  allele_labels_clusters  (path_key -> cluster_id, allele) for ASSIGNED paths
  L0c  alleles_clusters.tsv     class -> carrier haplotypes (moderns)

Output: data/derived/old_dups/cluster_class_map.tsv
  locus, level (cluster|private), cluster_id, n_catalog_paths,
  n_assigned_paths, class_label, n_carrier_haps

Reading it:
- 'cluster' rows: one per MDS cluster; class_label is the renumbering of that
  cluster (by DESCENDING ASSIGNED-path count, not catalog-path count - so the
  biggest catalog cluster is not always class 1: GOLGA6D, LPA, KATNAL2,
  HLA-DQA1, CHST5, C2orf78, ANKRD36C flip this way).
- 'private' rows: one per DBSCAN noise path; each noise path became its OWN
  class (never merged, never NA) and can still recruit many unrelated carriers.
"""
import csv
import glob
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MDS = Path('/moosefs/laurap/old_dups/cosi_like/mds')
AL = ROOT / 'data/derived/old_dups/allele_labels_clusters.tsv'
AC = ROOT / 'data/derived/old_dups/alleles_clusters.tsv'
BED = Path('/moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed')


def main():
    catalog = {}
    for f in MDS.glob('paths_*.tsv'):
        loc = f.stem.replace('paths_', '')
        catalog[loc] = {r['path']: int(r['cluster']) for r in csv.DictReader(open(f), delimiter='\t')}
    asg = defaultdict(dict)
    for r in csv.DictReader(open(AL), delimiter='\t'):
        asg[r['locus']][r['path_key']] = (int(r['cluster_id']), int(r['allele']))
    carriers = defaultdict(Counter)
    for r in csv.DictReader(open(AC), delimiter='\t'):
        if r['sample'][:2] in ('HG', 'NA'):
            for a in (r['allele1'], r['allele2']):
                if a != 'NA':
                    carriers[(r['chrom'], r['pos'])][a] += 1
    bed = defaultdict(list)
    for line in open(BED):
        f = line.split()
        if len(f) >= 4:
            c = f[0].split('#')[-1].replace('chr', '')
            bed[f[3]].append((c, str((int(f[1]) + int(f[2])) // 2)))
    locus_keys = {}
    for loc, keys in bed.items():
        for k in keys:
            if k in carriers:
                locus_keys[loc] = k
                break

    out = ROOT / 'data/derived/old_dups/cluster_class_map.tsv'
    with open(out, 'w') as fh:
        fh.write('locus\tlevel\tcluster_id\tn_catalog_paths\tn_assigned_paths\tclass_label\tn_carrier_haps\n')
        n_real = n_priv = 0
        for loc in sorted(catalog):
            key = locus_keys.get(loc)
            if key is None or loc not in asg:
                continue
            cat, A = catalog[loc], asg[loc]
            cl_info = defaultdict(lambda: [0, None])
            for pk, (cid, al) in A.items():
                if cid >= 0:
                    cl_info[cid][0] += 1
                    cl_info[cid][1] = al
            for cid, v in cl_info.items():
                n_cat = sum(1 for cc in cat.values() if cc == cid)
                nc = carriers[key].get(str(v[1]), 0)
                fh.write(f'{loc}\tcluster\t{cid}\t{n_cat}\t{v[0]}\t{v[1]}\t{nc}\n')
                n_real += 1
            for pk, (cid, al) in sorted(A.items()):
                if cid < 0:
                    nc = carriers[key].get(str(al), 0)
                    fh.write(f'{loc}\tprivate\t-1\t1\t1\t{al}\t{nc}\n')
                    n_priv += 1
    print(f'wrote {out}: {n_real} cluster rows, {n_priv} private rows')


if __name__ == '__main__':
    main()
