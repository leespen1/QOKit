#!/bin/bash
#SBATCH --job-name=Nc_BA
#SBATCH --array=0-32
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/Nc_BA_%A_%a.out
#SBATCH --error=logs/Nc_BA_%A_%a.err

set -euo pipefail

# -----------------------------
# Parameter grids
# -----------------------------
NODES=(10 11 12 13 14 15 16 17 18 19 20)
EDGES_PER_NODE=(1 2 3)

NUM_NODES=${#NODES[@]}
NUM_M=${#EDGES_PER_NODE[@]}

# -----------------------------
# Map SLURM_ARRAY_TASK_ID -> (n, m)
# -----------------------------
TASK_ID=${SLURM_ARRAY_TASK_ID}

NODE_INDEX=$(( TASK_ID / NUM_M ))
M_INDEX=$(( TASK_ID % NUM_M ))

NUMNODES=${NODES[$NODE_INDEX]}
EDGESPERNODE=${EDGES_PER_NODE[$M_INDEX]}

# -----------------------------
# Fixed parameters
# -----------------------------
SEEDSTART=0
NUMGRAPHS=100
FORMAT=hdf5
BACKEND=auto

# -----------------------------
# Sanity check: BA constraint 1 <= m < n
# -----------------------------
if (( EDGESPERNODE >= NUMNODES )); then
    echo "Invalid BA parameters: m=${EDGESPERNODE} >= n=${NUMNODES}"
    exit 1
fi

# -----------------------------
# Run
# -----------------------------
echo "Running N(c) collection (Barabasi–Albert)"
echo "  numNodes       = ${NUMNODES}"
echo "  edgesPerNode   = ${EDGESPERNODE}"
echo "  seedStart      = ${SEEDSTART}"
echo "  numGraphs      = ${NUMGRAPHS}"

python Nc_BarabasiAlbert.py \
    -n "${NUMNODES}" \
    -m "${EDGESPERNODE}" \
    -s "${SEEDSTART}" \
    -g "${NUMGRAPHS}" \
    -b "${BACKEND}" \
    -f "${FORMAT}"
