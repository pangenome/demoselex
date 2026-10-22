#!/bin/bash
# Segmented parallel download of the 22 per-chromosome 1000G phased panel VCFs.
# Each file is split into NSEG byte-range parts fetched concurrently; parts are
# size-verified, concatenated, then the final file size is verified. Resumable
# at part granularity (a part with a matching .done marker is skipped).
set -u
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
V=/moosefs/laurap/old_dups/demoselex_exp/demoselex/data/derived/g1k/vcf
NSEG=8
cd "$V" || exit 1

seg_dl() {
  local c=$1
  local base="1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel"
  local url="$B/$base.vcf.gz"
  local want=$(curl -sI --max-time 60 "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
  [ -n "$want" ] || { echo "chr$c: no content-length"; return 1; }
  if [ -f "$base.vcf.gz" ] && [ "$(stat -c %s "$base.vcf.gz")" = "$want" ]; then
    fetch_tbi "$c" "$base"; echo "chr$c: already complete"; return 0
  fi
  local seg=$(( (want + NSEG - 1) / NSEG ))
  local pids=()
  for i in $(seq 0 $((NSEG - 1))); do
    [ -f "$base.part$i.done" ] && continue
    (
      local lo=$((i * seg)) hi=$(( (i + 1) * seg - 1 ))
      [ $hi -ge $want ] && hi=$((want - 1))
      local size=$((hi - lo + 1))
      for attempt in 1 2 3 4 5 6; do
        curl -sSL --retry 6 --retry-delay 3 -r "$lo-$hi" -o "$base.part$i" "$url"
        if [ "$(stat -c %s "$base.part$i" 2>/dev/null || echo 0)" = "$size" ]; then
          touch "$base.part$i.done"; return 0
        fi
        sleep 3
      done
      echo "chr$c part$i: FAILED ($size bytes wanted)"
      return 1
    ) &
    pids+=($!)
  done
  local rc=0
  for p in "${pids[@]}"; do wait "$p" || rc=1; done
  [ $rc -eq 0 ] || { echo "chr$c: parts failed, will retry next round"; return 1; }
  cat "$base".part* > "$base.vcf.gz.tmp" 2>/dev/null || return 1
  if [ "$(stat -c %s "$base.vcf.gz.tmp")" = "$want" ]; then
    mv "$base.vcf.gz.tmp" "$base.vcf.gz"
    rm -f "$base".part* "$base".part*.done
    fetch_tbi "$c" "$base" || return 1
    echo "chr$c: OK ($((want / 1048576)) MB)"
    return 0
  fi
  echo "chr$c: concat size mismatch, redoing parts"
  rm -f "$base".part* "$base".part*.done "$base.vcf.gz.tmp"
  return 1
}

fetch_tbi() {
  local c=$1 base=$2
  local wt=$(curl -sI --max-time 60 "$B/$base.vcf.gz.tbi" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
  for a in 1 2 3; do
    curl -sSL --retry 6 -o "$base.vcf.gz.tbi" "$B/$base.vcf.gz.tbi"
    [ -n "$wt" ] && [ "$(stat -c %s "$base.vcf.gz.tbi" 2>/dev/null || echo 0)" = "$wt" ] && return 0
    sleep 3
  done
  echo "chr$c: tbi download problem"; return 1
}
export -f seg_dl fetch_tbi
export B V NSEG

# 6 rounds max; each round tries every chromosome not yet complete. Concurrency
# kept at 4 files x 8 parts = 32 connections: EBI's ftp refuses bursts above that.
for round in 1 2 3 4 5 6; do
  todo=()
  for c in $(seq 1 22); do
    base="1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel"
    want=$(curl -sI --max-time 60 "$B/$base.vcf.gz" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
    if [ -f "$base.vcf.gz" ] && [ "$(stat -c %s "$base.vcf.gz")" = "$want" ]; then
      [ -f "$base.vcf.gz.tbi" ] || fetch_tbi "$c" "$base"   # complete VCF, missing index
    else
      todo+=("$c")
    fi
  done
  [ ${#todo[@]} -eq 0 ] && break
  echo "=== round $round: ${#todo[@]} chromosomes to fetch: ${todo[*]}"
  printf '%s\n' "${todo[@]}" | xargs -P 4 -I{} bash -c 'seg_dl {}'
  [ ${#todo[@]} -eq 0 ] || sleep 45
done
n=$(ls *.vcf.gz 2>/dev/null | wc -l)
t=$(ls *.vcf.gz.tbi 2>/dev/null | wc -l)
echo "=== downloads finished: $n/22 VCFs, $t indexes"
[ "$n" = 22 ] && [ "$t" = 22 ]
