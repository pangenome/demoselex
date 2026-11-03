#!/bin/bash
# Genome-wide build of per-locus combined VCFs (phased 1000G panel records +
# likegt/HPRC structural class records, phased in VCF hap space via phase_link).
BASE=/moosefs/laurap/old_dups/demoselex_exp/demoselex
OUT=$BASE/data/derived/old_dups/combined_vcf
PY=$BASE/.venv/bin/python
mkdir -p "$OUT"
export BASE OUT PY

run_one() {
  loc=$1
  [ -s "$OUT/$loc.combined.vcf.gz" ] && { echo "skip: $loc"; return 0; }
  if "$PY" "$BASE/analysis/build_combined_vcf.py" --locus "$loc" > "$OUT/$loc.build.log" 2>&1; then
    echo "done: $loc ($(gzip -dc "$OUT/$loc.combined.vcf.gz" | grep -c '<CLS') class records)"
  else
    echo "FAILED: $loc (see $OUT/$loc.build.log)"
  fi
}
export -f run_one

cut -f4 /moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed | sort -u > /tmp/vcf_loci.txt
cat /tmp/vcf_loci.txt | xargs -P 10 -I{} bash -c 'run_one {}'
echo "COMBINED VCF ALL DONE"
