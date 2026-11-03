#!/bin/bash
# Genome-wide sv_phase_pilot run (erik's script verbatim) on the delivered full-genome ARG.
# Tables: pilot_run/tables/chr<N>.csv (moderns only, ARG cohort) -> pilot_run/chr<N>.csv
PY=/moosefs/laurap/old_dups/demoselex_exp/demoselex/.venv/bin/python
BUNDLE=/moosefs/erikg/demoselex-arg
BASE=/moosefs/laurap/old_dups/demoselex_exp/pilot_run
export G1K_ARG_DIR=$BUNDLE/arg
cd "$BUNDLE" || exit 1
run_one() {
  c=$1
  out=$BASE/chr${c}.csv
  [ -s "$out" ] && { echo "chr$c: already done, skip"; return 0; }
  if "$PY" sv_phase_pilot.py --table "$BASE/tables/chr${c}.csv" --out "$out" > "$BASE/chr${c}.log" 2>&1; then
    n=$(tail -n +2 "$out" | wc -l)
    echo "done: chr$c ($n class-rows)"
  else
    echo "FAILED: chr$c (see $BASE/chr${c}.log)"
  fi
}
export -f run_one; export PY BUNDLE BASE
printf '%s\n' $(seq 1 20) 22 | xargs -P 21 -I{} bash -c 'run_one {}'
echo "ALL CHROMOSOMES DONE"
