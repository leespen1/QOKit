#!/bin/bash
#SBATCH --job-name=Nc_WS
#SBATCH --array=0-131
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/Nc_WS_%A_%a.out
#SBATCH --error=logs/Nc_WS_%A_%a.err

set -euo pipefail

# -----------------------------
# Parameter grids
# -----------------------------
NODES=(10 11 12 13 14 15 16 17 18 19 20)
KVALUES=(4 6 8)
BETAS=(0.05 0.10 0.20 0.30)

NUM_NODES=${#NODES[@]}
NUM_K=${#KVALUES[@]}
NUM_BETA=${#BETAS[@]}

# -----------------------------
# Map SLURM_ARRAY_TASK_ID -> (n, k, beta)
# -----------------------------
TASK_ID=${SLURM_ARRAY_TASK_ID}

NODE_INDEX=$(( TASK_ID / (NUM_K * NUM_BETA) ))
REM=$(( TASK_ID % (NUM_K * NUM_BETA) ))
K_INDEX=$(( REM / NUM_BETA ))
BETA_INDEX=$(( REM % NUM_BETA ))

NUMNODES=${NODES[$NODE_INDEX]}
K=${KVALUES[$K_INDEX]}
BETA=${BETAS[$BETA_INDEX]}

# -----------------------------
# Fixed parameters
# -----------------------------
SEEDSTART=0
NUMGRAPHS=100
FORMAT=hdf5
BACKEND=auto

# -----------------------------
# Sanity checks (Watts–Strogatz constraints)
# -----------------------------
if (( K >= NUMNODES )); then
    echo "Invalid WS parameters: k=${K} >= n=${NUMNODES}"
    exit 1
fi

if (( K % 2 != 0 )); then
    echo "Invalid WS parameters: k=${K} must be even"
    exit 1
fi

# -----------------------------
# Run
# -----------------------------
echo "Running N(c) collection (Watts–Strogatz)"
echo "  numNodes            = ${NUMNODES}"
echo "  nearestNeighbors k  = ${K}"
echo "  rewiringProbability = ${BETA}"
echo "  seedStart           = ${SEEDSTART}"
echo "  numGraphs           = ${NUMGRAPHS}"

python Nc_WattsStrogatz.py \
    -n "${NUMNODES}" \
    -k "${K}" \
    -p "${BETA}" \
    -s "${SEEDSTART}" \
    -g "${NUMGRAPHS}" \
    -b "${BACKEND}" \
    -f "${FORMAT}"
