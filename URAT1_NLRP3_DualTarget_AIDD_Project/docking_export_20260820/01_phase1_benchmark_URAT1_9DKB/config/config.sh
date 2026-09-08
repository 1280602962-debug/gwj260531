#!/usr/bin/env bash
# Copy to config.sh and edit paths for the target server.
# shellcheck disable=SC2034

export PROTOCOL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export INPUT_DIR="${PROTOCOL_ROOT}/inputs"
export WORK_DIR="${PROTOCOL_ROOT}/work"
export RESULT_DIR="${PROTOCOL_ROOT}/results"
export SHARD_DIR="${PROTOCOL_ROOT}/shards"
export LOG_DIR="${PROTOCOL_ROOT}/logs"

# ---- docking params ----
export EXHAUSTIVENESS=32
export NUM_MODES=9
export SEED=42
export CPU_PER_TASK="${CPU_PER_TASK:-4}"
# Cap BLAS/OpenMP/Torch threads so local multi-shard runs do not oversubscribe.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${CPU_PER_TASK}}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-${CPU_PER_TASK}}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-${CPU_PER_TASK}}"
export TORCH_NUM_THREADS="${TORCH_NUM_THREADS:-${CPU_PER_TASK}}"

# Per-molecule wall-clock timeout (seconds). Exceeded → skip & log, continue.
# Typical Vina exh32 ~1 min; 5 min covers outliers without long stalls.
export VINA_TIMEOUT_SEC="${VINA_TIMEOUT_SEC:-300}"      # 5 min
export GNINA_TIMEOUT_SEC="${GNINA_TIMEOUT_SEC:-300}"    # 5 min
export RTM_TIMEOUT_SEC="${RTM_TIMEOUT_SEC:-3600}"      # 1 h for batch rescore chunk

# 9DKB / lesinurad site (same as redock smoke)
export CENTER_X=99.980
export CENTER_Y=102.958
export CENTER_Z=105.657
export SIZE_X=20.000
export SIZE_Y=20.010
export SIZE_Z=20.000

export VINA_BIN="${VINA_BIN:-/usr/bin/vina}"
# Blackwell-capable local build (sm_120); wrapper sets LD_LIBRARY_PATH / OpenBabel.
export GNINA_BIN="${GNINA_BIN:-/home/hww/gwj/NLRP3_URAT1/software/gnina_5090.sh}"
export OBABEL_BIN="${OBABEL_BIN:-obabel}"
export PYTHON_BIN="${PYTHON_BIN:-python3}"
export TIMEOUT_BIN="${TIMEOUT_BIN:-timeout}"  # GNU coreutils

# GPU CNN with local sm_120 build (system /usr/bin/gnina is still incompatible).
export GNINA_NO_GPU="${GNINA_NO_GPU:-0}"
# Comma-separated physical GPUs for round-robin (torch uses CUDA_VISIBLE_DEVICES).
export GNINA_DEVICES="${GNINA_DEVICES:-0}"
export GNINA_DEVICE="${GNINA_DEVICE:-0}"
# Parallel caps: vina uses CPU first; after vina finishes, gnina takes the CPUs.
export VINA_MAX_JOBS="${VINA_MAX_JOBS:-8}"
# Post-vina default (search is still CPU-heavy; CNN rescore uses GPU ~10GB/job).
export GNINA_MAX_JOBS="${GNINA_MAX_JOBS:-6}"

export RTMSCORE_PY="${RTMSCORE_PY:-/home/dakki/repos/RTMScore/example/rtmscore.py}"
export RTMSCORE_MODEL="${RTMSCORE_MODEL:-/home/dakki/repos/RTMScore/trained_models/rtmscore_model1.pth}"
# Prefer explicit env python (conda activate is unreliable in non-interactive launchers).
export RTM_PYTHON_BIN="${RTM_PYTHON_BIN:-/opt/anaconda3/envs/rtmscore/bin/python}"
export RTMSCORE_ENV_ACTIVATE="${RTMSCORE_ENV_ACTIVATE:-source /home/dakki/repos/RTMScore/env.sh}"
export RTM_CHUNK_SIZE="${RTM_CHUNK_SIZE:-50}"  # molecules per incremental RTMScore job
# While vina is running, keep RTM light so docking keeps CPU.
export RTM_OMP_THREADS="${RTM_OMP_THREADS:-2}"
export RTM_POLL_SEC="${RTM_POLL_SEC:-60}"

# On cluster with GPU + cuDNN for gnina, set:
#   export GNINA_NO_GPU=0
# and load CUDA module before running.

export RECEPTOR_PDB="${INPUT_DIR}/9DKB_receptor.pdb"
export RECEPTOR_PDBQT="${INPUT_DIR}/9DKB_receptor.pdbqt"
export AUTOBOX_LIGAND="${INPUT_DIR}/lesinurad_crystal_ref.sdf"
export RTM_POCKET_PDB="${INPUT_DIR}/9DKB_pocket_10.0.pdb"
export POOL_CSV="${INPUT_DIR}/unique_docking_pool.csv"
export TRUE_BENCH_CSV="${INPUT_DIR}/true_decoy_benchmark.csv"
export RANDOM_BENCH_CSV="${INPUT_DIR}/random_decoy_benchmark.csv"
export LIGAND_MANIFEST="${WORK_DIR}/ligand_manifest.csv"
