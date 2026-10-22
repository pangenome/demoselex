#!/usr/bin/env python3
"""Loci-only ARG rebuild: build only the 10 Mb dated chunks that a BED of loci touches.

The main build (g1k_phased.py) downloads a whole-chromosome phased VCF, tsinfers the
full chromosome, and tsdates every 10 Mb chunk. For a loci list covering ~6% of the
genome that is 94% wasted work. This script reproduces the *same dated chunks*, with
the same names, grid and coordinate convention, by:

  1. streaming the 1000G high-coverage phased VCF over HTTPS *for a window only*
     (chunk +/- --flank; htslib tabix range requests, nothing stored on disk),
  2. tsinfer on the window (flank = inference context),
  3. keeping the central 10 Mb, trim + simplify, tsdate inside_outside
     (exactly the operations g1k_phased.date applies to a full-chromosome ARG),
  4. dumping data/derived/g1k/chr<N>_<start>_<end>.dated.trees with positions
     relative to the chunk start, as chunk_for() in g1k_allele_analyzer.py expects.

Faithfulness caveat: the main build gives tsinfer whole-chromosome context; this build
gives it (10 + 2*flank) Mb. Local topology in the central chunk has `flank` Mb of
context on each side. Validate with the LCT selftest (g1k_allele_analyzer.py
--selftest) - the LCT chunk is added to every plan by default.

Site filters are copied verbatim from g1k_phased.build: SNVs only, panel AF with
min(AF, 1-AF) >= 0.05, duplicated positions skipped, genotypes = the two phased alleles.

Usage:
  g1k_loci_build.py plan  <loci.bed>              # write the chunk list
  g1k_loci_build.py popfile                      # write samples_population.txt (VCF order)
  g1k_loci_build.py chunk <chrom> <start>         # build one chunk (resumable)
  g1k_loci_build.py all <loci.bed> [--jobs 16]    # plan + build all chunks in parallel
"""
import argparse
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/g1k'
OUT.mkdir(parents=True, exist_ok=True)
URL = ('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/'
       '1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/'
       '1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz')
PED_URL = ('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/'
           '1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt')
POPFILE = OUT / 'samples_population.txt'
CHUNK = 10_000_000
MIN_MAF = 0.05
FLANK = 5_000_000
SELFTEST_CHUNK = (2, 135_000_000)       # LCT rs4988235 lives here; validation anchor
PLAN = OUT / 'loci_chunks.tsv'

# GRCh38 primary-sequence lengths (checked against the VCF header at chunk time)
CLEN = {1: 248_956_422, 2: 242_193_529, 3: 198_295_559, 4: 190_214_555,
        5: 181_538_259, 6: 170_805_979, 7: 159_345_973, 8: 145_138_636,
        9: 138_394_717, 10: 133_797_422, 11: 135_086_622, 12: 133_275_309,
        13: 114_364_328, 14: 107_043_718, 15: 101_991_189, 16: 90_338_345,
        17: 83_257_441, 18: 80_373_285, 19: 58_617_616, 20: 64_444_167,
        21: 46_709_983, 22: 50_818_468}


def norm_chrom(tok):
    """'GRCh38#0#chr7' / 'chr7' / '7' -> 7"""
    t = tok.split('#')[-1]
    t = t[3:] if t.startswith('chr') else t
    return int(t)


def chunks_for_bed(bed, with_selftest=True):
    starts = {}
    with open(bed) as f:
        for line in f:
            p = line.split()
            if len(p) < 3 or line.startswith(('#', 'track', 'browser')):
                continue
            c = norm_chrom(p[0])
            lo, hi = int(p[1]), int(p[2])
            for s in range((lo // CHUNK) * CHUNK, hi, CHUNK):
                starts.setdefault(c, set()).add(s)
    if with_selftest:
        starts.setdefault(SELFTEST_CHUNK[0], set()).add(SELFTEST_CHUNK[1])
    rows = []
    for c in sorted(starts):
        for s in sorted(starts[c]):
            rows.append((c, s, min(s + CHUNK, CLEN[c])))
    return rows


def cmd_plan(args):
    rows = chunks_for_bed(args.bed, with_selftest=not args.no_selftest)
    with open(PLAN, 'w') as f:
        f.write('chrom\tstart\tend\n')
        for c, s, e in rows:
            f.write(f'{c}\t{s}\t{e}\n')
    by = {}
    for c, s, e in rows:
        by[c] = by.get(c, 0) + 1
    print(f'{len(rows)} chunks to build ({sum(CHUNK for _ in rows) / 1e6:.0f} Mb):')
    for c in sorted(by):
        print(f'  chr{c}: {by[c]}')
    print(f'wrote {PLAN}')
    print(f'disk/network: no VCF stored; ~20 Mb region streamed per chunk, '
          f'{len(rows)} windows total')


def cmd_popfile(args):
    """samples_population.txt: the ped table reordered into VCF sample-column order."""
    import cyvcf2
    ped = {}
    order_in_ped = []
    with urllib.request.urlopen(PED_URL) as r:
        header = r.readline().decode().split()
        assert header[1] == 'SampleID', header
        for line in r.read().decode().splitlines():
            p = line.split()
            if len(p) < 7:
                continue
            ped[p[1]] = p
            order_in_ped.append(p[1])
    vcf = cyvcf2.VCF(URL.format(chrom=22))          # header read: range requests only
    vcf_order = list(vcf.samples)
    vcf.close()
    missing = [s for s in vcf_order if s not in ped]
    extra = [s for s in order_in_ped if s not in set(vcf_order)]
    if missing or extra:
        raise SystemExit(f'sample mismatch: {len(missing)} in VCF not in ped, '
                         f'{len(extra)} in ped not in VCF')
    with open(POPFILE, 'w') as f:
        f.write('#FamilyID SampleID FatherID MotherID Sex Population Superpopulation\n')
        for s in vcf_order:
            f.write('\t'.join(ped[s]) + '\n')
    print(f'wrote {POPFILE}: {len(vcf_order)} samples in VCF order')


def build_chunk(chrom, start, flank=FLANK, maf=MIN_MAF):
    end = min(start + CHUNK, CLEN[chrom])
    dst = OUT / f'chr{chrom}_{start}_{end}.dated.trees'
    if dst.exists():
        print(f'chr{chrom}:{start}-{end}: exists, skip', flush=True)
        return dst
    import cyvcf2
    local = (OUT / 'vcf' /
             f'1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz')
    url = str(local) if local.exists() else URL.format(chrom=chrom)
    vcf = cyvcf2.VCF(url)
    clen = dict(zip(vcf.seqnames, vcf.seqlens)).get(f'chr{chrom}')
    if clen and abs(clen - CLEN[chrom]) > 1_000:
        raise SystemExit(f'chr{chrom}: contig length {clen} != expected {CLEN[chrom]}')
    ws = max(0, start - flank)
    we = min(CLEN[chrom], start + CHUNK + flank)
    sd_path = OUT / f'loci_chr{chrom}_{start}.samples'
    import tsinfer as _tsi
    sd = None
    if sd_path.exists():                     # resume a previous interrupted run
        try:
            sd = _tsi.load(str(sd_path))
            if not sd.finalised:
                sd = None
        except Exception:
            sd = None
    if sd is None:
        sd_path.unlink(missing_ok=True)
        sd = _tsi.SampleData(path=str(sd_path), sequence_length=we - ws)
        for _ in vcf.samples:
            sd.add_individual(ploidy=2, time=0.0)
        n_sites = n_seen = 0
        last_pos = -1
        t0 = time.time()
        # capture htslib stderr: any [E:: there means a corrupt/truncated read -
        # delete the SampleData and fail so the scheduler retries from scratch
        import os as _os
        import tempfile as _tf
        errf = _tf.NamedTemporaryFile(prefix='htserr_', dir=str(OUT), delete=False)
        saved = _os.dup(2)
        _os.dup2(errf.fileno(), 2)
        stream_ok = True
        try:
            for rec in vcf(f'chr{chrom}:{ws + 1}-{we}'):
                n_seen += 1
                if rec.POS <= last_pos:
                    continue                       # duplicated position (multiallelic split)
                if len(rec.REF) != 1 or any(len(a) != 1 for a in rec.ALT):
                    continue                                   # SNVs only
                af = rec.INFO.get('AF')
                if isinstance(af, (list, tuple)):
                    af = af[0] if len(af) else None
                if af is None:
                    continue
                if min(af, 1 - af) < maf:
                    continue
                g = np.fromiter((a for gt in rec.genotypes for a in gt[:2]), dtype=np.int8,
                                count=2 * len(vcf.samples))
                sd.add_site(position=rec.POS - ws, genotypes=g,
                            alleles=[rec.REF] + list(rec.ALT))
                last_pos = rec.POS
                n_sites += 1
                if n_sites % 50000 == 0:
                    print(f'  chr{chrom}:{start} {n_sites} sites kept / {n_seen} seen '
                          f'({time.time() - t0:.0f}s)', flush=True)
        finally:
            _os.dup2(saved, 2)
            _os.close(saved)
            errf.flush()
            _os.close(errf.fileno())
        with open(errf.name, 'rb') as fh:
            hts_errors = fh.read()
        Path(errf.name).unlink(missing_ok=True)
        if b'[E::' in hts_errors:
            try:
                sd.close()
            except Exception:
                pass
            sd_path.unlink(missing_ok=True)
            raise SystemExit(f'chr{chrom}:{ws + 1}-{we}: htslib stream error '
                             f'({hts_errors.splitlines()[:1]}); SampleData deleted, retry me')
        sd.finalise()
        print(f'chr{chrom}:{ws}-{we}: {n_sites} sites kept of {n_seen} records, '
              f'SampleData {sd_path.stat().st_size / 2**20:.0f} MB ({time.time() - t0:.0f}s)',
              flush=True)
    n_sites = sd.num_sites
    sd.close()
    vcf.close()
    if n_sites == 0:
        sd_path.unlink(missing_ok=True)
        raise SystemExit(f'chr{chrom}:{start}-{end}: no sites kept')
    t1 = time.time()
    ts = _tsi.infer(_tsi.load(str(sd_path)), path_compression=True,
                   progress_monitor=False)
    print(f'chr{chrom}:{start}: inferred {ts.num_trees} trees, {ts.num_sites} sites, '
          f'{ts.num_samples} samples in {time.time() - t1:.0f}s', flush=True)
    sd_path.unlink(missing_ok=True)
    # central 10 Mb, exactly as g1k_phased.date() carves a chunk from the full ARG
    sub = ts.keep_intervals([[start - ws, min(start + CHUNK, we) - ws]]).trim().simplify()
    if sub.num_sites == 0:
        raise SystemExit(f'chr{chrom}:{start}-{end}: no sites in the central chunk')
    t2 = time.time()
    import tsdate
    dated = tsdate.date(sub, mutation_rate=1.25e-8, population_size=10_000,
                        method='inside_outside', time_units='generations')
    dated.dump(str(dst))
    print(f'chr{chrom}:{start}-{end}: dated {sub.num_sites} sites, '
          f'{dated.num_trees} trees in {time.time() - t2:.0f}s -> {dst.name}',
          flush=True)
    return dst


def cmd_chunk(args):
    build_chunk(int(args.chrom), int(args.start), flank=args.flank)


def cmd_all(args):
    if not PLAN.exists() or args.bed:
        ns = argparse.Namespace(bed=args.bed or None, no_selftest=False)
        cmd_plan(ns)
    rows = []
    with open(PLAN) as f:
        next(f)
        for line in f:
            c, s, e = line.split()
            dst = OUT / f'chr{c}_{s}_{e}.dated.trees'
            if not dst.exists():
                rows.append((int(c), int(s)))
    print(f'{len(rows)} chunks to build, {args.jobs} parallel jobs', flush=True)
    log = OUT / 'loci_build.log'
    t0 = time.time()

    def one(cs):
        c, s = cs
        lf = OUT / f'log_loci_chr{c}_{s}.txt'
        with open(lf, 'w') as f:
            r = subprocess.run([sys.executable, __file__, 'chunk', str(c), str(s),
                                '--flank', str(args.flank)],
                               stdout=f, stderr=subprocess.STDOUT)
        return c, s, r.returncode

    done = failed = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(one, cs) for cs in rows]
        for fut in as_completed(futs):
            c, s, rc = fut.result()
            done += 1
            if rc != 0:
                failed += 1
                print(f'FAILED chr{c}:{s} (rc={rc}, see log_loci_chr{c}_{s}.txt)',
                      flush=True)
            print(f'[{done}/{len(rows)}] chr{c}:{s} rc={rc} '
                  f'({(time.time() - t0) / 60:.1f} min elapsed)', flush=True)
    with open(log, 'a') as f:
        f.write(f'{time.strftime("%F %T")}: {done} finished, {failed} failed, '
                f'{(time.time() - t0) / 60:.1f} min\n')
    print(f'done: {done - failed} ok, {failed} failed, see {log}', flush=True)
    sys.exit(1 if failed else 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest='mode', required=True)
    p = sub.add_parser('plan'); p.add_argument('bed'); p.add_argument('--no-selftest', action='store_true')
    sub.add_parser('popfile')
    p = sub.add_parser('chunk'); p.add_argument('chrom'); p.add_argument('start')
    p.add_argument('--flank', type=int, default=FLANK)
    p = sub.add_parser('all'); p.add_argument('bed', nargs='?')
    p.add_argument('--jobs', type=int, default=16)
    p.add_argument('--flank', type=int, default=FLANK)
    a = ap.parse_args()
    if a.mode == 'plan':
        cmd_plan(a)
    elif a.mode == 'popfile':
        cmd_popfile(a)
    elif a.mode == 'chunk':
        cmd_chunk(a)
    elif a.mode == 'all':
        cmd_all(a)


if __name__ == '__main__':
    main()
