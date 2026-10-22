#!/bin/bash
# Stage the allele-analyzer prerequisites for collaborators.
#
#   ./copy_share.sh /moosefs/erikg/demoselex-arg                  # local path
#   ./copy_share.sh erikg@octopus01:/moosefs/erikg/demoselex-arg  # over ssh
#
# Copies the small files (analyzer, sample order, annotations, README) and the 294 dated ARG
# chunks (~58 GiB). Safe to re-run: rsync skips what is already there.
set -u
SRC="$(cd "$(dirname "$0")" && pwd)"            # .../data/share/g1k-arg
ARGDIR="${SRC}/../../../data/derived/g1k"        # repo data/derived/g1k
DEST="${1:-/moosefs/erikg/demoselex-arg}"

SMALL=(README.md g1k_allele_analyzer.py sv_phase_pilot.py sv_anchor_match.py
       samples_population.txt sv265_annotation.bed genomicSuperDups_hg38.txt.gz
       loci-analysis-note.md)

case "$DEST" in
  *:*)  HOST="${DEST%%:*}"; RPATH="${DEST#*:}" ;;
  *)    HOST="";           RPATH="$DEST" ;;
esac

echo "source : ${SRC} (ARG chunks: ${ARGDIR})"
echo "dest   : ${HOST:+${HOST}:}${RPATH}"

remote() { if [ -n "$HOST" ]; then ssh -o BatchMode=yes "$HOST" "$@"; else bash -c "$*"; fi; }

# 0. fail fast if the destination is not writable (a hung FUSE mount blocks everything)
if ! remote "mkdir -p '${RPATH}/arg' && test -w '${RPATH}'" 2>/dev/null; then
  echo "ERROR: cannot create/write ${HOST:+${HOST}:}${RPATH} (mount hung or no permission)." >&2
  exit 1
fi

# 1. small files
rsync -a --info=progress2 "${SMALL[@]/#/${SRC}/}" "${DEST}/" || exit 1

# 2. the ARG chunks
rsync -a --ignore-existing --info=progress2 "${ARGDIR}"/chr*_*_*.dated.trees "${DEST}/arg/" || exit 1

# 3. provenance, count, openness
remote "cd '${RPATH}' && sha256sum ${SMALL[*]} > SHA256SUMS.small 2>/dev/null; \
        ls -1 arg/chr*_*_*.dated.trees 2>/dev/null | wc -l > arg/FILE_COUNT.txt; \
        chmod -R a+rX '${RPATH}'; \
        echo \"chunks: \$(cat arg/FILE_COUNT.txt)\"; du -sh arg 2>/dev/null | tail -1"
echo "done"
