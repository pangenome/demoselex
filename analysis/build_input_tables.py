#!/usr/bin/env python3
"""Build the two INPUT CONTRACT tables (docs/loci-analysis-note.md) from the
old_dups genotyping outputs.

Outputs (under data/derived/old_dups/ in the repo):

  alleles.tsv   chrom pos sample allele1 allele2
               One row per (locus, sample), FULL GRID: moderns and ancients,
               NA where not called. Per locus, distinct structural-haplotype
               path IDs are mapped to integer labels 1..k (sorted, so the
               mapping is deterministic). Path identity = 'sample#copy#contig'
               (the ':start-end' interval suffix reported by likegt is the
               graph interval, not a stable identity, and is stripped).

  dates.tsv    sample date_bp date_sd source culture
               From old_dups_metadata-V1_dated.tsv (date_mean_BP, date_sd_BP,
               date_method, cultural_genetic_context); moderns get date 0.

  allele_labels.tsv   locus path_key allele   (the label map, for interpretation)

Sources:
  --bed        svs_refined.hgsvcv3.bed (GRCh38) - gene -> (chrom, pos=midpoint)
  --identity   gene_cn/all_copies_identity.tsv - modern graph paths
               (class GENE rows; haplotype sample#copy, contig sample#copy#ACC)
  --geno-dir   aDNA_real_1x_geno/results/*.geno1.tsv - ancient hap1/hap2 calls
  --metadata   old_dups_metadata-V1_dated.tsv - dates

Moderns: a sample's diplotype at a locus = its own #1 and #2 GENE-class paths
(the same leave-zero-out truth the likegt/COSIGT validation used). This yields
233 moderns (1000G/HGDP IDs, so they sit inside the dated ARG). The remaining
modern cohort (samples without assemblies in the graph) has no diplotype data
yet and is simply absent.

Usage:
  build_input_tables.py [--min-qv 0] [--ancients-only|--moderns-only]
"""
import argparse
import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/old_dups'
OLD_DUPS = Path('/moosefs/laurap/old_dups')


def norm_chrom(tok):
    t = tok.split('#')[-1]
    return int(t[3:]) if t.startswith('chr') else int(t)


def path_key(hap):
    """'HG01934#1#CM088841.1:48425716-48827166' -> 'HG01934#1#CM088841.1'"""
    if pd.isna(hap) or hap in ('', 'NA'):
        return None
    return str(hap).rsplit(':', 1)[0]


def load_bed(bed):
    """gene -> (chrom, pos=midpoint, start, end) (GRCh38)"""
    out = {}
    with open(bed) as f:
        for line in f:
            p = line.split()
            if len(p) < 4 or line.startswith(('#', 'track', 'browser')):
                continue
            gene = p[3]
            if gene in out:
                raise SystemExit(f'duplicate gene {gene} in BED')
            lo, hi = int(p[1]), int(p[2])
            out[gene] = (norm_chrom(p[0]), (lo + hi) // 2, lo, hi)
    return out


# Per-gene names in the CN tables that belong to merged loci in the SV BED. The
# merged-locus diplotype cannot be reconstructed from per-gene CN tables (a
# sample has up to 6 gene copies there, not 2 locus haplotypes), so moderns are
# skipped at these loci; ancient calls define the labels.
GENE_ALIASES = {'GYPA': 'GYP', 'GYPB': 'GYP', 'GYPE': 'GYP'}


def load_modern(identity):
    """gene -> {sample: (path_key_1, path_key_2)} for samples with both copies."""
    df = pd.read_csv(identity, sep='\t')
    df = df[(df['class'] == 'GENE') & df['haplotype'].str.contains('#')]
    df = df[df['haplotype'].str.match(r'^(HG|NA)\d+#\d+$')]
    df = df[~df['gene'].isin(GENE_ALIASES)]          # merged loci: moderns skipped
    df['pk'] = df['contig'].map(path_key)
    modern = {}
    for (gene, hap), g in df.groupby(['gene', 'haplotype'], sort=False):
        sample, copy = hap.split('#')
        if copy not in ('1', '2') or g['pk'].isna().any():
            continue
        pk = g['pk'].iloc[0]
        modern.setdefault(gene, {}).setdefault(sample, {})[int(copy)] = pk
    out = {}
    for gene, samples in modern.items():
        for s, copies in samples.items():
            if 1 in copies and 2 in copies:
                out.setdefault(gene, {})[s] = (copies[1], copies[2])
    return out


def load_ancient(geno_dir, min_qv):
    """locus -> [ (sample, pk1, pk2), ... ]"""
    rows = []
    files = sorted(Path(geno_dir).glob('*.geno1.tsv'))
    if not files:
        raise SystemExit(f'no *.geno1.tsv under {geno_dir}')
    bad = []
    for f in files:
        df = None
        for attempt in (1, 2, 3):                   # moosefs errno-5 is often transient
            try:
                df = pd.read_csv(f, sep='\t')
                break
            except OSError as e:
                if attempt == 3:
                    bad.append(str(f))
                    print(f'WARNING: unreadable, skipped: {f}: {e}', flush=True)
        if df is None:
            continue
        if min_qv:
            df = df[df['graph_qv'] >= min_qv]
        for r in df.itertuples():
            rows.append((r.sample, r.locus, path_key(r.hap1), path_key(r.hap2)))
    if bad:
        print(f'{len(bad)} geno files unreadable after retries; those samples get NA rows')
    anc = {}
    for s, locus, k1, k2 in rows:
        if k1 is None and k2 is None:
            continue
        anc.setdefault(locus, []).append((s, k1, k2))
    return anc, len(files) - len(bad)


def load_modern_geno(path, min_qv):
    """Modern read-based genotypes (geno_3202_master.tsv.gz): locus -> {sample: (pk1, pk2)}"""
    import gzip
    op = gzip.open if str(path).endswith('.gz') else open
    modern = {}
    n_rows = n_qv_drop = 0
    with op(path, 'rt') as f:
        rd = csv.DictReader(f, delimiter='\t')
        for r in rd:
            n_rows += 1
            if min_qv and r.get('graph_qv') and float(r['graph_qv']) < min_qv:
                n_qv_drop += 1
                continue
            k1, k2 = path_key(r.get('hap1')), path_key(r.get('hap2'))
            if k1 is None and k2 is None:
                continue
            modern.setdefault(r['locus'], {})[r['sample']] = (k1, k2)
    print(f'modern genotyping: {n_rows} rows, {n_qv_drop} dropped by min-qv, '
          f'{len(modern)} loci')
    return modern


def load_clusters(mds_dir):
    """locus -> (path_key -> cluster_id) and locus -> QC dict from mds_summary.tsv."""
    import glob as _glob
    import os as _os
    cl_of, qc = {}, {}
    for f in _glob.glob(str(Path(mds_dir) / 'paths_*.tsv')):
        locus = _os.path.basename(f)[len('paths_'):-4]
        d = {}
        for r in csv.DictReader(open(f), delimiter='\t'):
            d[r['path'].rsplit(':', 1)[0]] = int(r['cluster'])
        cl_of[locus] = d
    summ = Path(mds_dir) / 'mds_summary.tsv'
    if summ.exists():
        for r in csv.DictReader(open(summ), delimiter='\t'):
            qc[r['locus']] = r
    return cl_of, qc


def cluster_labels(gene, path_keys, cl_of, qc):
    """Label scheme: clusters by descending size 1..k, then noise paths (and any
    paths missing from the cluster map) get private labels k+1.. in sorted order.
    Returns (path_key -> label, n_clusters, qc_tuple or None)."""
    cl = cl_of.get(gene, {})
    sizes = {}
    for pk in path_keys:
        c = cl.get(pk)
        if c is not None and c >= 0:
            sizes[c] = sizes.get(c, 0) + 1
    ordered = sorted(sizes, key=lambda c: (-sizes[c], c))
    lab = {c: i + 1 for i, c in enumerate(ordered)}
    k = len(ordered)
    nxt = k
    for pk in sorted(path_keys):
        c = cl.get(pk)
        if c is None or c < 0:
            nxt += 1
            lab[pk] = nxt
    relab = {}
    for pk in path_keys:
        c = cl.get(pk)
        relab[pk] = lab[c] if (c is not None and c >= 0) else lab[pk]
    q = qc.get(gene)
    qc_t = (q['k'], q['silhouette'], q['noise_frac']) if q else None
    return relab, k, qc_t


def build(args):
    bed = load_bed(args.bed)
    if getattr(args, 'modern_geno', None):
        modern = load_modern_geno(args.modern_geno, args.min_qv)
    else:
        modern = load_modern(args.identity)
    ancient, n_files = load_ancient(args.geno_dir, args.min_qv)
    loci = sorted(set(modern) | set(ancient))
    no_bed = [g for g in loci if g not in bed]
    if no_bed:
        raise SystemExit(f'loci missing from BED: {no_bed}')
    print(f'{len(loci)} loci: {len(modern)} with modern paths, '
          f'{len(ancient)} with ancient calls ({n_files} ancient files)')

    OUT.mkdir(parents=True, exist_ok=True)
    cl_of, qc = ({}, {})
    if getattr(args, 'mds', None):
        cl_of, qc = load_clusters(args.mds)
        covered = sum(1 for g in loci if g in cl_of)
        print(f'cluster maps: {len(cl_of)} loci, covering {covered}/{len(loci)} of ours')
    # per-locus label map over the union of path keys
    keys, ncl_of, qcinfo = {}, {}, {}
    for gene in loci:
        ks = set()
        for s, (k1, k2) in modern.get(gene, {}).items():
            ks.update(k for k in (k1, k2) if k)
        for s, k1, k2 in ancient.get(gene, []):
            ks.update(k for k in (k1, k2) if k)
        if cl_of:
            relab, k, qc_t = cluster_labels(gene, ks, cl_of, qc)
            keys[gene] = relab
            ncl_of[gene] = k
            qcinfo[gene] = qc_t
        else:
            keys[gene] = {k: i + 1 for i, k in enumerate(sorted(ks))}

    suffix = '_clusters' if cl_of else ''
    n_mod = n_anc = n_na = 0
    all_anc_samples = sorted({t[0] for gene in ancient for t in ancient[gene]})
    with open(OUT / f'alleles{suffix}.tsv', 'w') as f:
        if cl_of:
            f.write('chrom\tpos\tsample\tallele1\tallele2\tpath1\tpath2\t'
                    'cluster_k\tsilhouette\tnoise_frac\tstart\tend\n')
        else:
            f.write('chrom\tpos\tsample\tallele1\tallele2\tstart\tend\n')
        for gene in loci:
            chrom, pos, lo, hi = bed[gene]
            lab = keys[gene]
            m = modern.get(gene, {})
            called = {t[0]: (t[1], t[2]) for t in ancient.get(gene, [])}
            samples = sorted(set(m) | set(called) | set(all_anc_samples))
            kk, sil, nfr = qcinfo.get(gene) or ('NA', 'NA', 'NA')
            for s in samples:
                a1 = a2 = 'NA'
                p1 = p2 = 'NA'
                if s in m:
                    k1, k2 = m[s]
                    a1 = lab[k1] if k1 else 'NA'
                    a2 = lab[k2] if k2 else 'NA'
                    p1, p2 = k1 or 'NA', k2 or 'NA'
                    n_mod += 1
                elif s in called:
                    k1, k2 = called[s]
                    a1 = lab[k1] if k1 else 'NA'
                    a2 = lab[k2] if k2 else 'NA'
                    p1, p2 = k1 or 'NA', k2 or 'NA'
                    n_anc += 1
                else:
                    n_na += 1
                row = f'{chrom}\t{pos}\t{s}\t{a1}\t{a2}'
                if cl_of:
                    row += f'\t{p1}\t{p2}\t{kk}\t{sil}\t{nfr}'
                row += f'\t{lo}\t{hi}'
                f.write(row + '\n')
    print(f'wrote {OUT / f"alleles{suffix}.tsv"}: {n_mod} modern rows, {n_anc} ancient rows, '
          f'{n_na} NA rows')

    with open(OUT / f'allele_labels{suffix}.tsv', 'w') as f:
        if cl_of:
            f.write('locus\tpath_key\tcluster_id\tallele\n')
        else:
            f.write('locus\tpath_key\tallele\n')
        for gene in loci:
            clg = cl_of.get(gene, {})
            for k, a in keys[gene].items():
                cid = clg.get(k, -1)
                row = f'{gene}\t{k}\t{a}'
                if cl_of:
                    row = f'{gene}\t{k}\t{cid}\t{a}'
                f.write(row + '\n')

    md = pd.read_csv(args.metadata, sep='\t')
    md = md[['sample_id', 'date_mean_BP', 'date_sd_BP', 'date_method',
             'cultural_genetic_context']]
    modern_ids = sorted({s for gene in modern for s in modern[gene]})
    extra = pd.DataFrame({'sample_id': modern_ids, 'date_mean_BP': 0,
                          'date_sd_BP': 0, 'date_method': 'modern',
                          'cultural_genetic_context': 'modern'})
    dates = pd.concat([extra, md], ignore_index=True).rename(columns={
        'sample_id': 'sample', 'date_mean_BP': 'date_bp', 'date_sd_BP': 'date_sd',
        'date_method': 'source', 'cultural_genetic_context': 'culture'})
    dates.to_csv(OUT / 'dates.tsv', sep='\t', index=False)
    nd = (dates['date_bp'] > 0).sum()
    print(f'wrote {OUT / "dates.tsv"}: {len(dates)} samples, {nd} with BP dates')
    if getattr(args, 'modern_geno', None):
        print('moderns are READ-BASED genotypes (geno_3202): the full 3202-sample ARG cohort')
    else:
        print('note: moderns carry their true graph diplotypes (leave-zero-out truth); '
              'other moderns are absent until genotyped from reads')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bed', default=str(OLD_DUPS / 'svs_refined.hgsvcv3.bed'))
    ap.add_argument('--identity', default=str(OLD_DUPS / 'gene_cn/all_copies_identity.tsv'))
    ap.add_argument('--geno-dir', default=str(OLD_DUPS / 'aDNA_real_1x_geno/results'))
    ap.add_argument('--metadata', default=str(OLD_DUPS / 'old_dups_metadata-V1_dated.tsv'))
    ap.add_argument('--min-qv', type=float, default=0.0)
    ap.add_argument('--modern-geno',
                    default='/moosefs/raid5/lpignata/old_dups/geno_3202_fa/geno_3202_master.tsv.gz',
                    help='read-based modern genotypes (all 3202); empty string = graph-truth moderns')
    ap.add_argument('--mds', default='/moosefs/laurap/old_dups/cosi_like/mds',
                    help='dir with paths_<locus>.tsv cluster maps; empty string = path-key labels')
    args = ap.parse_args()
    if args.modern_geno == '':
        args.modern_geno = None
    if args.mds == '':
        args.mds = None
    build(args)


if __name__ == '__main__':
    main()
