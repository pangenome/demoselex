#!/bin/bash
# Genome-wide structural-allele-class analysis over the loci-only ARG.
#
# Configurations (per the cohort discussion):
#   - global mixed cohort (primary) at tau 200 / 500 / 1000
#   - per-superpopulation (AFR AMR EAS EUR SAS) at tau 200
# All runs use --anchors 9 (multi-anchor aggregation across each locus span).
#
# The script waits until all 167 dated chunks exist, then runs one analyzer
# process per chromosome (22 parallel), merging per-configuration CSVs at the end.
# Ancients are not in the ARG cohort; they are excluded here (temporal layer
# will use them) - each per-superpop table keeps only that superpop's moderns.
#
# Usage: run_genome_wide.sh [--no-wait]   (JOBS=22 ANCH=9 env overrides)
set -u
R=/moosefs/laurap/old_dups/demoselex_exp/demoselex
PY=$R/.venv/bin/python
TAB=$R/data/derived/old_dups/alleles_clusters.tsv
POPF=$R/data/derived/g1k/samples_population.txt
OUT=$R/data/derived/old_dups/genome_wide
JOBS=${JOBS:-22}
ANCH=${ANCH:-9}
WAIT=1
[ "${1:-}" = "--no-wait" ] && WAIT=0

mkdir -p "$OUT/by_chr"

# ---- 0. wait for the full ARG ---------------------------------------------
if [ "$WAIT" = 1 ]; then
  while [ "$(ls $R/data/derived/g1k/*.dated.trees 2>/dev/null | wc -l)" -lt 167 ]; do
    sleep 120
  done
  echo "$(date +%T): all 167 chunks present, starting"
fi

# ---- 1. split the table by chromosome --------------------------------------
if [ ! -s "$OUT/by_chr/chr1.tsv" ]; then
  awk -F'\t' -v d="$OUT/by_chr" '
    NR==1 {h=$0; next}
    !seen[$1]++ {print h > d"/chr"$1".tsv"}
    {print >> d"/chr"$1".tsv"}' "$TAB"
  echo "table split by chromosome: $(ls $OUT/by_chr/chr*.tsv | wc -l) files"
fi

# ---- 2. per-superpopulation modern ID lists --------------------------------
for P in AFR AMR EAS EUR SAS; do
  awk -v p="$P" 'NR>1 && $7==p{print $2}' "$POPF" > "$OUT/by_chr/$P.ids"
done

# ---- 3. one configuration: name, tau, idfile ('' = global) ------------------
run_config () {
  local name=$1 tau=$2 ids=$3
  mkdir -p "$OUT/$name"
  for c in $(seq 1 22); do
    (
      if [ -n "$ids" ]; then
        awk -F'\t' \
            'NR==FNR{m[$1]=1; next} FNR==1 || $3 in m' \
            "$OUT/by_chr/$ids.ids" "$OUT/by_chr/chr$c.tsv" > "$OUT/$name/table_chr$c.tsv"
        T="$OUT/$name/table_chr$c.tsv"
      else
        T="$OUT/by_chr/chr$c.tsv"
      fi
      "$PY" "$R/analysis/g1k_allele_analyzer.py" --table "$T" --tau "$tau" \
          --anchors "$ANCH" --out "$OUT/$name/chr$c.csv" \
          > "$OUT/$name/chr$c.log" 2>&1
      rm -f "$OUT/$name/table_chr$c.tsv"
    ) &
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
  done
  wait
  # merge per-chr CSVs (skip empty chromosomes)
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

# ---- 4. run all configurations ---------------------------------------------
run_config global_tau200  200  ''
run_config global_tau500  500  ''
run_config global_tau1000 1000 ''
for P in AFR AMR EAS EUR SAS; do
  run_config "${P}_tau200" 200 "$P"
done
echo "$(date +%T): ALL CONFIGURATIONS DONE"
