#!/bin/bash
# Genome-wide phase-link run: per-locus orientation of likegt hap1/hap2 onto the
# ARG/VCF phased haplotype space (validated on ZNF146 pilot).
BASE=/moosefs/laurap/old_dups/demoselex_exp/demoselex
OUT=/moosefs/laurap/old_dups/demoselex_exp/phase_link_run
PY=$BASE/.venv/bin/python
mkdir -p "$OUT" /tmp/phase_link
export G1K_ARG_DIR=/moosefs/erikg/demoselex-arg/arg

run_one() {
  loc=$1
  [ -s "$OUT/$loc.csv" ] && { echo "skip: $loc"; return 0; }
  if [ ! -f "/moosefs/raid5/lpignata/old_dups/gfa_paf_sv/$loc.gfa" ]; then
    echo "FAILED: $loc (no GFA)"; return 1
  fi
  if "$PY" "$BASE/analysis/phase_link.py" --locus "$loc" --out "$OUT/$loc.csv" \
      --workdir /tmp/phase_link > "$OUT/$loc.log" 2>&1; then
    n=$(tail -n +2 "$OUT/$loc.csv" 2>/dev/null | wc -l)
    echo "done: $loc ($n rows)"
  else
    echo "FAILED: $loc (see $OUT/$loc.log)"
  fi
}
export -f run_one; export BASE OUT PY

cut -f4 /moosefs/laurap/old_dups/svs_refined.hgsvcv3.bed | sort -u > "$OUT/loci.txt"
cat "$OUT/loci.txt" | xargs -P 16 -I{} bash -c 'run_one {}'
echo "PHASE-LINK ALL LOCI DONE"
