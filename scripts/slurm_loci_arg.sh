#!/bin/bash
# SLURM array: one task = one 10 Mb dated-ARG chunk for the SV loci set.
# Chunk list: data/derived/g1k/slurm_chunks.txt (chrom start), resumable.
#SBATCH --job-name=loci_arg
#SBATCH --partition=workers
#SBATCH --array=0-158%72
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=/moosefs/laurap/old_dups/demoselex_exp/demoselex/data/derived/g1k/slurm/%A_%a.out
set -u
cd /moosefs/laurap/old_dups/demoselex_exp/demoselex
read -r C S <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" data/derived/g1k/slurm_chunks.txt)"
echo "=== chr$C:$S on $(hostname) [task $SLURM_ARRAY_TASK_ID]"
exec .venv/bin/python analysis/g1k_loci_build.py chunk "$C" "$S" --flank 5000000
