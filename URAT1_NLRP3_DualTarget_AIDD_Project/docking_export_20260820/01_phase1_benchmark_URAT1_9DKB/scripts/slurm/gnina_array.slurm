#!/usr/bin/env bash
#SBATCH -J gnina_pi
#SBATCH -o logs/gnina_%A_%a.out
#SBATCH -e logs/gnina_%A_%a.err
#SBATCH --array=0-196%16
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 48:00:00
# Prefer GPU partition; then GNINA_NO_GPU=0 in config.sh
# #SBATCH -p gpu
# #SBATCH --gres=gpu:1

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"
# shellcheck disable=SC1091
source ./config.sh
mkdir -p logs
SHARD=$(printf "shards/shard_%04d.txt" "${SLURM_ARRAY_TASK_ID}")
bash scripts/run_gnina_shard.sh "${SHARD}"
