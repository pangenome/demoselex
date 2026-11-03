#!/usr/bin/env python3
"""Phase link: determine the correspondence between likegt's hap1/hap2 path
assignments and the ARG / VCF phased haplotype space, per (sample, locus).

Method (validated on ZNF146 pilot, 2026-10-25):
- likegt hap numbering is arbitrary per (sample, locus): read-based, no
  genome-wide phasing input. Confirmed: S/SW ~50/50 on het samples.
- The assigned path is a real HPRC donor sequence segment (PanSN name,
  embedded GRCh38 span). Reconstruct its sequence from the per-locus GFA,
  align to the GRCh38 locus window (minimap2 -x asm20 -c), and read its base
  at every ARG site where the sample's two phased haplotypes differ.
- The 2x2 match (path1/path2 x hapA/hapB) gives the orientation: S (likegt
  h1 = VCF hap 1) or SW (swapped), with concordance = 1 - (donor-sample
  divergence at het sites); orientation decided only when the score gap is
  decisive. AMBIG/UND are reported, never guessed.

Inputs: geno_3202_master.tsv.gz, per-locus GFAs, MDS cluster maps,
        the delivered ARG chunks (erik's), GRCh38 FASTA.
Output: one row per (sample, locus) with both haps called:
        sample, locus, chrom, lo, hi,
        cluster_h1, cluster_h2 (likegt numbering, as in alleles_clusters),
        orient (S/SW/AMBIG/UND),
        arg_cluster_hap0, arg_cluster_hap1 (cluster labels on ARG hap 2i / 2i+1),
        n_sites, conc_same, conc_sw (best-diagonal concordances)

Usage: phase_link.py --locus ZNF146 | --all [--jobs N] [--out phase_link.csv]
"""
import argparse
import bisect
import csv
import gzip
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import tskit

MASTER = Path('/moosefs/raid5/lpignata/old_dups/geno_3202_fa/geno_3202_master.tsv.gz')
GFA_DIR = Path('/moosefs/raid5/lpignata/old_dups/gfa_paf_sv')
MDS_DIR = Path('/moosefs/laurap/old_dups/cosi_like/mds')
BED = Path('/moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed')
ARG_DIR = Path('/moosefs/erikg/demoselex-arg/arg')
POPFILE = Path('/moosefs/erikg/demoselex-arg/samples_population.txt')
REF = Path('/moosefs/laurap/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa')
SAMTOOLS = '/home/lpignata/miniconda3/envs/likegt/bin/samtools'
MINIMAP2 = '/home/lpignata/miniconda3/envs/likegt/bin/minimap2'
CHUNK_RE = re.compile(r'chr(\d+)_(\d+)_(\d+)\.dated\.trees')


def chunks_for(chrom, lo, hi):
    out = []
    for f in sorted(ARG_DIR.glob(f'chr{chrom}_*.dated.trees')):
        m = CHUNK_RE.match(f.name)
        if m and int(m.group(2)) <= hi and lo < int(m.group(3)):
            out.append((f, int(m.group(2))))
    return out


def load_arg_sites(chrom, lo, hi, order_idx):
    """ARG variant data over [lo, hi): positions, allele tuples, 6404 x n matrix."""
    sites = []
    for f, start in chunks_for(chrom, lo, hi):
        ts = tskit.load(str(f))
        for v in ts.variants():
            p = start + int(v.position)
            if lo <= p < hi:
                sites.append((p, v.alleles, v.genotypes))
    sites.sort(key=lambda x: x[0])
    if not sites:
        return None
    P = np.array([p for p, a, g in sites])
    A = {p: a for p, a, g in sites}
    G = np.array([g for p, a, g in sites]).T
    return P, A, G


def gfa_paths(needed):
    """Stream-parse the locus GFA once; return {path_name: sequence}."""
    gfa = GFA_DIR / f'{LOCUS}.gfa'
    seqs, out = {}, {}
    with open(gfa) as f:
        for line in f:
            if line.startswith('S\t'):
                p = line.rstrip('\n').split('\t')
                seqs[p[1]] = p[2]
            elif line.startswith('P\t'):
                p = line.rstrip('\n').split('\t')
                if p[1] in needed:
                    s = []
                    for tok in p[2].split(','):
                        tok = tok.strip()
                        if not tok:
                            continue
                        nid, orient = tok[:-1], tok[-1]
                        sq = seqs[nid]
                        if orient == '-':
                            sq = sq[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
                        s.append(sq)
                    out[p[1]] = ''.join(s)
    return out


def paf_blocks(paf_path):
    per = {}
    for line in open(paf_path):
        f = line.rstrip('\n').split('\t')
        qname, qstart, tstart = f[0], int(f[2]), int(f[7])
        cigar = None
        for tag in f[12:]:
            if tag.startswith('cg:Z:'):
                cigar = tag[5:]
        if cigar is None:
            continue
        tpos, qpos = tstart, qstart
        for n, op in re.findall(r'(\d+)([MIDNSH=X])', cigar):
            n = int(n)
            if op in 'M=X':
                per.setdefault(qname, []).append((tpos, qpos, n)); tpos += n; qpos += n
            elif op in 'DN':
                tpos += n
            elif op in 'IS':
                qpos += n
    return per


def run_locus(locus, bed_span, master_rows, cluster_of, order_idx, workdir):
    chrom, lo, hi = bed_span
    rows = []
    # samples where orientation matters: different cluster labels on the two haps
    het, hom = [], []
    for s, h1, h2 in master_rows:
        c1, c2 = cluster_of.get(h1), cluster_of.get(h2)
        if c1 is None and c2 is None:
            continue
        (het if c1 != c2 else hom).append((s, h1, h2, c1, c2))
    if not het:
        return []
    arg = load_arg_sites(chrom, lo, hi, order_idx)
    if arg is None:
        for s, h1, h2, c1, c2 in het:
            rows.append(dict(sample=s, locus=locus, chrom=chrom, lo=lo, hi=hi,
                             cluster_h1=c1, cluster_h2=c2, orient='NO-COVERAGE',
                             arg_cluster_hap0=None, arg_cluster_hap1=None,
                             n_sites=0, n_scored=0, conc_same=None, conc_sw=None))
        return rows
    P, A, G = arg
    # ref window
    with tempfile.TemporaryDirectory(dir=workdir) as td:
        ref_fa = Path(td) / 'ref.fa'
        region = f'chr{chrom}:{lo - 5000}-{hi + 5000}'
        out = subprocess.run([SAMTOOLS, 'faidx', str(REF), region],
                             capture_output=True, text=True)
        if out.returncode != 0:
            raise RuntimeError(f'faidx failed: {out.stderr[:200]}')
        ref_fa.write_text(out.stdout)
        wstart = lo - 5000
        # paths (all assigned paths, even cluster-less ones: their sequence is
        # what the orientation match runs on)
        need = {h for s, h1, h2, c1, c2 in het for h in (h1, h2)}
        pseqs = gfa_paths(need)
        pseqs_u = {k.replace('#', '_'): v for k, v in pseqs.items()}
        paths_fa = Path(td) / 'paths.fa'
        with open(paths_fa, 'w') as f:
            for name, s in pseqs.items():
                f.write(f'>{name.replace("#", "_")}\n{s}\n')
        paf_path = Path(td) / 'paths.paf'
        with open(paf_path, 'w') as f:
            subprocess.run([MINIMAP2, '-x', 'asm20', '-c', '-t', '4',
                            str(ref_fa), str(paths_fa)], stdout=f,
                           stderr=subprocess.DEVNULL)
        blocks = paf_blocks(str(paf_path))
        starts = {q: [b[0] for b in bs] for q, bs in blocks.items()}

        def base_at(qname, tpos):
            bs = blocks.get(qname)
            if not bs:
                return None
            i = bisect.bisect_right(starts[qname], tpos) - 1
            if i < 0:
                return None
            ts_, qs_, n = bs[i]
            off = tpos - ts_
            if off >= n:
                return None
            b = pseqs_u[qname][qs_ + off].upper()
            return b if b in 'ACGT' else None

        for s, h1, h2, c1, c2 in het:
            k = order_idx.get(s)
            if k is None:
                continue
            hA, hB = G[2 * k], G[2 * k + 1]
            inf = [j for j in range(len(P)) if hA[j] != hB[j]]
            p1, p2 = h1.replace('#', '_'), h2.replace('#', '_')

            def conc(path_name, hap):
                m = t = 0
                for j in inf:
                    b = base_at(path_name, int(P[j]) - wstart)
                    if b is None:
                        continue
                    a = A[int(P[j])]
                    if b not in a:
                        continue
                    t += 1
                    if a.index(b) == hap[j]:
                        m += 1
                return (m / t if t else None), t

            c11, t11 = conc(p1, hA)
            c22, t22 = conc(p2, hB)
            c12, t12 = conc(p1, hB)
            c21, t21 = conc(p2, hA)
            n_ok = min(t11 or 0, t22 or 0, t12 or 0, t21 or 0)
            if n_ok < 15 or c11 is None or c22 is None or c12 is None or c21 is None:
                orient, c_same, c_sw = 'UND', None, None
            else:
                c_same, c_sw = (c11 + c22) / 2, (c12 + c21) / 2
                if abs(c_same - c_sw) < 0.10:
                    orient = 'AMBIG'
                elif c_same > c_sw:
                    orient = 'S'
                else:
                    orient = 'SW'
            if orient == 'S':
                a0, a1 = c1, c2
            elif orient == 'SW':
                a0, a1 = c2, c1
            else:
                a0, a1 = None, None
            rows.append(dict(sample=s, locus=locus, chrom=chrom, lo=lo, hi=hi,
                             cluster_h1=c1, cluster_h2=c2, orient=orient,
                             arg_cluster_hap0=a0, arg_cluster_hap1=a1,
                             n_sites=len(inf), n_scored=n_ok,
                             conc_same=round(c_same, 3) if c_same is not None else None,
                             conc_sw=round(c_sw, 3) if c_sw is not None else None))
    return rows


def main():
    global LOCUS
    ap = argparse.ArgumentParser()
    ap.add_argument('--locus')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--out', default='phase_link.csv')
    ap.add_argument('--workdir', default='/tmp/phase_link')
    a = ap.parse_args()
    Path(a.workdir).mkdir(parents=True, exist_ok=True)

    bed = {}
    for line in open(BED):
        f = line.split()
        if len(f) >= 4:
            c = f[0].split('#')[-1].replace('chr', '')
            bed[f[3]] = (int(c), int(f[1]), int(f[2]))

    cluster_of = {}
    for f in MDS_DIR.glob('paths_*.tsv'):
        locus = f.stem.replace('paths_', '')
        for r in csv.DictReader(open(f), delimiter='\t'):
            cluster_of.setdefault(r['path'], r['cluster']) if False else None
            cluster_of[f'{locus}\x1f{r["path"]}'] = r['cluster']
    master = defaultdict(list)
    with gzip.open(MASTER, 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            master[r['locus']].append((r['sample'], r['hap1'], r['hap2']))
    order = [l.split()[1] for l in open(POPFILE)][1:]
    order_idx = {s: i for i, s in enumerate(order)}

    def cluster_for(locus, path):
        return cluster_of.get(f'{locus}\x1f{path}')

    loci = sorted([a.locus]) if a.locus else sorted(bed)
    out_rows = []
    for locus in loci:
        if locus not in bed or locus not in master:
            continue
        LOCUS = locus
        rows_sub = [(s, h1, h2) for s, h1, h2 in master[locus]]
        cl = {}
        for s, h1, h2 in rows_sub:
            cl[h1] = cluster_for(locus, h1)
            cl[h2] = cluster_for(locus, h2)
        het_rows = run_locus(locus, bed[locus], rows_sub, cl, order_idx, a.workdir)
        out_rows.extend(het_rows)
        print(f'{locus}: {len(het_rows)} het rows', file=sys.stderr, flush=True)
    if out_rows:
        with open(a.out, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
            w.writeheader()
            w.writerows(out_rows)
        print(f'wrote {a.out}: {len(out_rows)} rows')


if __name__ == '__main__':
    main()
