#!/bin/bash --login
#SBATCH --job-name=Nc_BA
#SBATCH --array=0-32
#SBATCH --time=02:00:00
#SBATCH --nodes=1                    # Maximum number of nodes to be allocated
#SBATCH --ntasks-per-node=1          # Maximum number of tasks on each node
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/Nc_BA_%A_%a.out
#SBATCH --error=logs/Nc_BA_%A_%a.err
#SBATCH --constraint=amd20
#SBATCH --export=ALL

set -euo pipefail

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    module purge
    module load CUDA
fi

source ../../../qokitvenv/bin/activate

# -----------------------------
# Parameter grids
# Override from outside, e.g.:
#   NODES="10 12 15 20" EDGES_PER_NODE="2 3" sbatch slurm_Nc_BarabasiAlbert.slurm
# -----------------------------
NODES=(${NODES:-10 11 12 13 14 15 16 17 18 19 20})           # override: NODES="10 15 20"
EDGES_PER_NODE=(${EDGES_PER_NODE:-1 2 3})                   # override: EDGES_PER_NODE="2 3"

# -----------------------------
# Fixed parameters
# Override from outside via inline prefix or export:
#   STORE_COSTS=true bash slurm_Nc_BarabasiAlbert.slurm
#   export STORE_COSTS=true && sbatch slurm_Nc_BarabasiAlbert.slurm
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
    # Manual batching: bash slurm_Nc_BarabasiAlbert.slurm <batchId> <numBatches>
    BATCH_ID=$1
    NUM_BATCHES=$2
else
    BATCH_ID=0
    NUM_BATCHES=1
fi

echo "=== Barabasi-Albert Parameters ==="
echo "  NODES:            ${NODES[*]}"
echo "  EDGES_PER_NODE:   ${EDGES_PER_NODE[*]}"
echo "  SEEDSTART:        ${SEEDSTART}"
echo "  NUMGRAPHS:        ${NUMGRAPHS}"
echo "  FORMAT:           ${FORMAT}"
echo "  BACKEND:          ${BACKEND}"
echo "  STORE_COSTS:      ${STORE_COSTS}"
echo "  COMPUTE_HOMODIST: ${COMPUTE_HOMODIST}"
echo "  BATCH_ID:         ${BATCH_ID}"
echo "  NUM_BATCHES:      ${NUM_BATCHES}"
echo "=================================="

OPTIONAL_FLAGS=()
[[ "${STORE_COSTS}" == "true" ]] && OPTIONAL_FLAGS+=(-c)
[[ "${COMPUTE_HOMODIST}" == "true" ]] && OPTIONAL_FLAGS+=(-d)

python Nc_BarabasiAlbert.py \
    -n "${NODES[@]}" \
    -m "${EDGES_PER_NODE[@]}" \
    -s "${SEEDSTART}" -g "${NUMGRAPHS}" \
    -b "${BACKEND}" -f "${FORMAT}" \
    --numBatches "${NUM_BATCHES}" --batchId "${BATCH_ID}" \
    "${OPTIONAL_FLAGS[@]+"${OPTIONAL_FLAGS[@]}"}"
