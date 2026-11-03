#!/usr/bin/env python3
"""Combined VCFs: 1000G phased SNPs/indels + likegt/HPRC-derived structural
class calls, per locus, with the structural phase in VCF haplotype space.

Per locus <L>, writes <L>.combined.vcf.gz:
- all phased-panel records within [lo-5kb, hi+5kb] (bcftools view), unchanged;
- one biallelic symbolic record per structural class k, at the locus midpoint:
    CHROM  mid  <L>.CL<k>  <refbase>  <CLS:k>  . PASS  LOCUS=<L>;CL=<k>;SRC=likegt-MDS
  GT per sample, in VCF column order:
    phased '1|0' / '0|1' / '1|1' / '0|0' from the phase-link orientation
      ('|.' / '.|' when one haplotype is called and the other NA),
    unphased '0/1' etc. when the orientation is AMBIG/UND (never guessed),
    './.' when the sample has no call at the locus.

Inputs: data/derived/old_dups/alleles_clusters.tsv (labels, likegt order),
        phase_link_run/<LOCUS>.csv (orientations in VCF/ARG space),
        data/derived/g1k/vcf/<chr>.vcf.gz (phased panel),
        GRCh38 FASTA for the reference base at the midpoint.
Usage: build_combined_vcf.py --locus ZNF146 | --all [--jobs N]
"""
import argparse
import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AC = ROOT / 'data/derived/old_dups/alleles_clusters.tsv'
PL = Path('/moosefs/laurap/old_dups/demoselex_exp/phase_link_run')
VCF_DIR = ROOT / 'data/derived/g1k/vcf'
REF = Path('/moosefs/laurap/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa')
BED = Path('/moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed')
BCFTOOLS = None  # not installed; tabix used instead
BGZIP = '/home/lpignata/miniconda3/bin/bgzip'
TABIX = '/home/lpignata/miniconda3/bin/tabix'
SAMTOOLS = '/home/lpignata/miniconda3/envs/likegt/bin/samtools'


def load_locus(locus, bed):
    chrom, lo, hi = bed[locus]
    labels = {}   # sample -> (l1, l2) in likegt order
    classes = set()
    out = subprocess.run(['awk', '-F', '\t', '-v', f'c={chrom}', '-v', f'lo={lo}', '-v', f'hi={hi}',
                          '$1==c && $11==lo && $12==hi', str(AC)],
                         capture_output=True, text=True)
    for line in out.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 7:
            continue
        s, a1, a2 = f[2], f[3], f[4]
        labels[s] = (a1, a2)
        for a in (a1, a2):
            if a != 'NA':
                classes.add(a)
    return chrom, lo, hi, labels, sorted(classes, key=int)


def phase_map(locus):
    m = {}
    f = PL / f'{locus}.csv'
    if f.exists():
        for r in csv.DictReader(open(f)):
            if r['orient'] in ('S', 'SW'):
                m[r['sample']] = (r['orient'], r['arg_cluster_hap0'], r['arg_cluster_hap1'])
            else:
                m[r['sample']] = (r['orient'], None, None)
    return m


def build_one(locus, bed, outdir, samples_order, idx_of):
    chrom, lo, hi, labels, classes = load_locus(locus, bed)
    if not classes:
        return None
    pm = phase_map(locus)
    mid = (lo + hi) // 2
    # reference base at midpoint
    reg = f'chr{chrom}:{mid}-{mid}'
    base = subprocess.run([SAMTOOLS, 'faidx', str(REF), reg],
                          capture_output=True, text=True).stdout.splitlines()
    refbase = base[1].upper() if len(base) > 1 else 'N'
    vcf_in = VCF_DIR / f'1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
    if not vcf_in.exists():
        return None
    region = f'chr{chrom}:{lo - 5000}-{hi + 5000}'
    # stream records; insert class records at the midpoint position
    proc = subprocess.Popen([TABIX, '-h', str(vcf_in), region],
                            stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    out_path = Path(outdir) / f'{locus}.combined.vcf'
    n_ins = 0
    header = ['##INFO=<ID=LOCUS,Number=1,Type=String,Description="locus name">',
              '##INFO=<ID=CL,Number=1,Type=Integer,Description="likegt/MDS structural class label">',
              '##ALT=<ID=CLS,Description="likegt HPRC-path structural class carrier state per phased haplotype">']
    with open(out_path, 'w') as out:
        inserted = False
        for line in proc.stdout:
            if line.startswith('##'):
                out.write(line)
                continue
            if line.startswith('#CHROM'):
                for h in header:
                    out.write(h + '\n')
                out.write(line)
                continue
            pos = int(line.split('\t', 2)[1])
            if not inserted and pos > mid:
                for k in classes:
                    fmt = []
                    for s in samples_order:
                        l1, l2 = labels.get(s, ('NA', 'NA'))
                        pm_s = pm.get(s)
                        if l1 == 'NA' and l2 == 'NA':
                            gt = './.'
                        elif l1 == l2:
                            gt = ('1|1' if l1 == k else '0|0')
                        else:
                            if pm_s and pm_s[0] in ('S', 'SW'):
                                hA, hB = (pm_s[1], pm_s[2])
                                g = []
                                for h in (hA, hB):
                                    if h == 'NA':
                                        g.append('.')
                                    else:
                                        g.append('1' if h == k else '0')
                                gt = g[0] + '|' + g[1]
                            else:
                                a1 = '.' if l1 == 'NA' else ('1' if l1 == k else '0')
                                a2 = '.' if l2 == 'NA' else ('1' if l2 == k else '0')
                                gt = a1 + '/' + a2
                        fmt.append(gt)
                    out.write('\t'.join([f'chr{chrom}', str(mid), f'{locus}.CL{k}',
                                         refbase, f'<CLS:{k}>', '.', 'PASS',
                                         f'LOCUS={locus};CL={k}', 'GT'] + fmt) + '\n')
                    n_ins += 1
                inserted = True
            out.write(line)
        if not inserted:
            for k in classes:
                fmt = []
                for s in samples_order:
                    l1, l2 = labels.get(s, ('NA', 'NA'))
                    pm_s = pm.get(s)
                    if l1 == 'NA' and l2 == 'NA':
                        gt = './.'
                    elif l1 == l2:
                        gt = ('1|1' if l1 == k else '0|0')
                    elif pm_s and pm_s[0] in ('S', 'SW'):
                        g = []
                        for h in (pm_s[1], pm_s[2]):
                            g.append('.' if h == 'NA' else ('1' if h == k else '0'))
                        gt = g[0] + '|' + g[1]
                    else:
                        a1 = '.' if l1 == 'NA' else ('1' if l1 == k else '0')
                        a2 = '.' if l2 == 'NA' else ('1' if l2 == k else '0')
                        gt = a1 + '/' + a2
                    fmt.append(gt)
                out.write('\t'.join([f'chr{chrom}', str(mid), f'{locus}.CL{k}',
                                     refbase, f'<CLS:{k}>', '.', 'PASS',
                                     f'LOCUS={locus};CL={k}', 'GT'] + fmt) + '\n')
                n_ins += 1
    proc.wait()
    # bgzip + index
    subprocess.run([BGZIP, '-f', str(out_path)], check=True)
    subprocess.run([TABIX, '-p', 'vcf', str(out_path) + '.gz'], check=True)
    return n_ins


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--locus')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--outdir', default=str(ROOT / 'data/derived/old_dups/combined_vcf'))
    a = ap.parse_args()
    bed = {}
    for line in open(BED):
        f = line.split()
        if len(f) >= 4:
            c = f[0].split('#')[-1].replace('chr', '')
            bed[f[3]] = (int(c), int(f[1]), int(f[2]))
    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    # VCF sample order = popfile order
    chrom0 = 1
    vcf = VCF_DIR / f'1kGP_high_coverage_Illumina.chr{chrom0}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
    hdr = subprocess.run([TABIX, '-H', str(vcf)], capture_output=True, text=True).stdout
    samples = [l for l in hdr.splitlines() if l.startswith('#CHROM')][0].split('\t')[9:]
    loci = sorted([a.locus]) if a.locus else sorted(bed)
    for locus in loci:
        if locus not in bed:
            continue
        n = build_one(locus, bed, a.outdir, samples, None)
        print(f'{locus}: {n} class records', flush=True)


if __name__ == '__main__':
    main()
