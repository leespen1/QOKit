#!/bin/bash --login
#SBATCH --job-name=Nc_ER
#SBATCH --array=0-54
#SBATCH --time=02:00:00
#SBATCH --nodes=1                    # Maximum number of nodes to be allocated
#SBATCH --ntasks-per-node=1          # Maximum number of tasks on each node
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/Nc_ER_%A_%a.out
#SBATCH --error=logs/Nc_ER_%A_%a.err
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
# -----------------------------
NODES=(10 11 12 13 14 15 16 17 18 19 20)
#NODES=(10 11 12)
PROBS=(0.1 0.25 0.5 0.75 0.9)

# -----------------------------
# Fixed parameters
# -----------------------------
SEEDSTART="${SEEDSTART:-0}"
NUMGRAPHS="${NUMGRAPHS:-100}"
FORMAT="${FORMAT:-hdf5}"
BACKEND="${BACKEND:-auto}"
STORE_COSTS="${STORE_COSTS:-false}"      # Store raw cost of each bitstring (-c)
COMPUTE_HOMODIST="${COMPUTE_HOMODIST:-false}" # Compute homogeneous distribution N(c',d,c) via Julia (-d)

# -----------------------------
# Determine batch ID and count
# -----------------------------
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    BATCH_ID=$(( SLURM_ARRAY_TASK_ID - SLURM_ARRAY_TASK_MIN ))
    NUM_BATCHES=$(( SLURM_ARRAY_TASK_MAX - SLURM_ARRAY_TASK_MIN + 1 ))
elif [[ -n "${2:-}" ]]; then
    # Manual batching: bash slurm_Nc_ErdosRenyi.slurm <batchId> <numBatches>
    BATCH_ID=$1
    NUM_BATCHES=$2
else
    BATCH_ID=0
    NUM_BATCHES=1
fi

echo "=== Erdos-Renyi Parameters ==="
echo "  NODES:            ${NODES[*]}"
echo "  PROBS:            ${PROBS[*]}"
echo "  SEEDSTART:        ${SEEDSTART}"
echo "  NUMGRAPHS:        ${NUMGRAPHS}"
echo "  FORMAT:           ${FORMAT}"
echo "  BACKEND:          ${BACKEND}"
echo "  STORE_COSTS:      ${STORE_COSTS}"
echo "  COMPUTE_HOMODIST: ${COMPUTE_HOMODIST}"
echo "  BATCH_ID:         ${BATCH_ID}"
echo "  NUM_BATCHES:      ${NUM_BATCHES}"
echo "================================"

OPTIONAL_FLAGS=()
[[ "${STORE_COSTS}" == "true" ]] && OPTIONAL_FLAGS+=(-c)
[[ "${COMPUTE_HOMODIST}" == "true" ]] && OPTIONAL_FLAGS+=(-d)

python Nc_ErdosRenyi.py \
    -n "${NODES[@]}" \
    -p "${PROBS[@]}" \
    -s "${SEEDSTART}" -g "${NUMGRAPHS}" \
    -b "${BACKEND}" -f "${FORMAT}" \
    --numBatches "${NUM_BATCHES}" --batchId "${BATCH_ID}" \
    "${OPTIONAL_FLAGS[@]+"${OPTIONAL_FLAGS[@]}"}"
