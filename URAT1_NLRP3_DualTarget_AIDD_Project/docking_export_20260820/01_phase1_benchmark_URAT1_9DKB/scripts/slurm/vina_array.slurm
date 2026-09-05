#!/usr/bin/env bash
#SBATCH -J vina_pi
#SBATCH -o logs/vina_%A_%a.out
#SBATCH -e logs/vina_%A_%a.err
#SBATCH --array=0-196%32
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 24:00:00
# Edit --array after make_shards.py. Edit partition/account for your site.

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"
# shellcheck disable=SC1091
source ./config.sh
mkdir -p logs
SHARD=$(printf "shards/shard_%04d.txt" "${SLURM_ARRAY_TASK_ID}")
bash scripts/run_vina_shard.sh "${SHARD}"
