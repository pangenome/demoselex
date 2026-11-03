#!/bin/bash
# Genome-wide allele-class analysis on erik's delivered full-genome ARG
# (/moosefs/erikg/demoselex-arg/arg/, 294 whole-genome chunks), using his
# bundle analyzer (identical to branch version + argdir auto-resolution).
# Same 8 configurations as run_genome_wide.sh: global x tau{200,500,1000}
# + per-superpopulation x tau200, all --anchors 9. Reuses the by-chromosome
# table split and superpop ID lists from the first run.
set -u
R=/moosefs/laurap/old_dups/demoselex_exp/demoselex
PY=$R/.venv/bin/python
ANALYZER=/moosefs/erikg/demoselex-arg/g1k_allele_analyzer.py
SPLIT=$R/data/derived/old_dups/genome_wide/by_chr
OUT=$R/data/derived/old_dups/genome_wide_erikarg
JOBS=${JOBS:-22}
ANCH=${ANCH:-9}

mkdir -p "$OUT"

run_config () {
  local name=$1 tau=$2 ids=$3
  mkdir -p "$OUT/$name"
  for c in $(seq 1 22); do
    (
      if [ -n "$ids" ]; then
        awk -F'\t' 'NR==FNR{m[$1]=1; next} FNR==1 || $3 in m' \
            "$SPLIT/$ids.ids" "$SPLIT/chr$c.tsv" > "$OUT/$name/table_chr$c.tsv"
        T="$OUT/$name/table_chr$c.tsv"
      else
        T="$SPLIT/chr$c.tsv"
      fi
      "$PY" "$ANALYZER" --table "$T" --tau "$tau" \
          --anchors "$ANCH" --out "$OUT/$name/chr$c.csv" \
          > "$OUT/$name/chr$c.log" 2>&1
      rm -f "$OUT/$name/table_chr$c.tsv"
    ) &
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
  done
  wait
  : > "$OUT/$name/summary.csv"
  first=1
  for c in $(seq 1 22); do
    f="$OUT/$name/chr$c.csv"
    [ -s "$f" ] || continue
    if [ "$first" = 1 ]; then cat "$f" >> "$OUT/$name/summary.csv"; first=0
    else tail -n +2 "$f" >> "$OUT/$name/summary.csv"; fi
  done
  echo "$(date +%T): $name done: $(($(wc -l < "$OUT/$name/summary.csv") - 1)) allele rows"
}

run_config global_tau200  200  ''
run_config global_tau500  500  ''
run_config global_tau1000 1000 ''
for P in AFR AMR EAS EUR SAS; do
  run_config "${P}_tau200" 200 "$P"
done
echo "$(date +%T): ERIK-ARG ALL CONFIGURATIONS DONE"
