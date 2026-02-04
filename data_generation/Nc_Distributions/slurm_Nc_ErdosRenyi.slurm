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

module purge
module load CUDA
source ../../../qokitvenv/bin/activate

# -----------------------------
# Parameter grids
# -----------------------------
NODES=(10 11 12 13 14 15 16 17 18 19 20)
PROBS=(0.1 0.25 0.5 0.75 0.9)

NUM_NODES=${#NODES[@]}
NUM_PROBS=${#PROBS[@]}

# -----------------------------
# Map SLURM_ARRAY_TASK_ID -> (n, p)
# -----------------------------
TASK_ID=${SLURM_ARRAY_TASK_ID}

NODE_INDEX=$(( TASK_ID / NUM_PROBS ))
PROB_INDEX=$(( TASK_ID % NUM_PROBS ))

NUMNODES=${NODES[$NODE_INDEX]}
EDGEPROBABILITY=${PROBS[$PROB_INDEX]}

# -----------------------------
# Fixed parameters
# -----------------------------
SEEDSTART=0
NUMGRAPHS=100
FORMAT=hdf5
BACKEND=auto

# -----------------------------
# Run
# -----------------------------
echo "Running N(c) collection"
echo "  numNodes        = ${NUMNODES}"
echo "  edgeProbability = ${EDGEPROBABILITY}"
echo "  seedStart       = ${SEEDSTART}"
echo "  numGraphs       = ${NUMGRAPHS}"

python Nc_ErdosRenyi.py \
    -n "${NUMNODES}" \
    -p "${EDGEPROBABILITY}" \
    -s "${SEEDSTART}" \
    -g "${NUMGRAPHS}" \
    -b "${BACKEND}" \
    -f "${FORMAT}"
