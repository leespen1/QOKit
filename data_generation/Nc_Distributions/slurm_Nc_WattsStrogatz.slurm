#!/bin/bash --login
#SBATCH --job-name=Nc_WS
#SBATCH --array=0-131
#SBATCH --time=02:00:00
#SBATCH --nodes=1                    # Maximum number of nodes to be allocated
#SBATCH --ntasks-per-node=1          # Maximum number of tasks on each node
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/Nc_WS_%A_%a.out
#SBATCH --error=logs/Nc_WS_%A_%a.err
#SBATCH --constraint=amd20

# Fail fast
set -euo pipefail

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    module purge
    module load CUDA
fi

source ../../../qokitvenv/bin/activate

# -----------------------------
# Parameter grids
# Override from outside, e.g.:
#   NODES="10 15 20" KVALUES="4 6" BETAS="0.1 0.2" sbatch slurm_Nc_WattsStrogatz.slurm
# -----------------------------
NODES=(${NODES:-10 11 12 13 14 15 16 17 18 19 20})           # override: NODES="10 15 20"
KVALUES=(${KVALUES:-4 6 8})                                 # override: KVALUES="4 6"
BETAS=(${BETAS:-0.05 0.10 0.20 0.30})                      # override: BETAS="0.1 0.2"

# -----------------------------
# Fixed parameters
# Override from outside via inline prefix or export:
#   STORE_COSTS=true bash slurm_Nc_WattsStrogatz.slurm
#   export STORE_COSTS=true && sbatch slurm_Nc_WattsStrogatz.slurm
# -----------------------------
SEEDSTART="${SEEDSTART:-0}"
NUMGRAPHS="${NUMGRAPHS:-100}"
FORMAT="${FORMAT:-hdf5}"
BACKEND="${BACKEND:-auto}"
STORE_COSTS="${STORE_COSTS:-false}"      # Store raw cost of each bitstring (-c)
COMPUTE_HOMODIST="${COMPUTE_HOMODIST:-false}" # Compute homogeneous distribution N(c',d,c) via Julia (-d)
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-1}" # Avoid JuliaCall segfault with multithreading

# -----------------------------
# Determine batch ID and count
# -----------------------------
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    BATCH_ID=$(( SLURM_ARRAY_TASK_ID - SLURM_ARRAY_TASK_MIN ))
    NUM_BATCHES=$(( SLURM_ARRAY_TASK_MAX - SLURM_ARRAY_TASK_MIN + 1 ))
elif [[ -n "${2:-}" ]]; then
    # Manual batching: bash slurm_Nc_WattsStrogatz.slurm <batchId> <numBatches>
    BATCH_ID=$1
    NUM_BATCHES=$2
else
    BATCH_ID=0
    NUM_BATCHES=1
fi

echo "=== Watts-Strogatz Parameters ==="
echo "  NODES:            ${NODES[*]}"
echo "  KVALUES:          ${KVALUES[*]}"
echo "  BETAS:            ${BETAS[*]}"
echo "  SEEDSTART:        ${SEEDSTART}"
echo "  NUMGRAPHS:        ${NUMGRAPHS}"
echo "  FORMAT:           ${FORMAT}"
echo "  BACKEND:          ${BACKEND}"
echo "  STORE_COSTS:      ${STORE_COSTS}"
echo "  COMPUTE_HOMODIST: ${COMPUTE_HOMODIST}"
echo "  BATCH_ID:         ${BATCH_ID}"
echo "  NUM_BATCHES:      ${NUM_BATCHES}"
echo "==================================="

OPTIONAL_FLAGS=()
[[ "${STORE_COSTS}" == "true" ]] && OPTIONAL_FLAGS+=(-c)
[[ "${COMPUTE_HOMODIST}" == "true" ]] && OPTIONAL_FLAGS+=(-d)

python Nc_WattsStrogatz.py \
    -n "${NODES[@]}" \
    -k "${KVALUES[@]}" \
    -p "${BETAS[@]}" \
    -s "${SEEDSTART}" -g "${NUMGRAPHS}" \
    -b "${BACKEND}" -f "${FORMAT}" \
    --numBatches "${NUM_BATCHES}" --batchId "${BATCH_ID}" \
    "${OPTIONAL_FLAGS[@]+"${OPTIONAL_FLAGS[@]}"}"
