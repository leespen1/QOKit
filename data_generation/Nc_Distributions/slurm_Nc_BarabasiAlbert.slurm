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
EDGES_PER_NODE=(1 2 3)

NUM_NODES=${#NODES[@]}
NUM_M=${#EDGES_PER_NODE[@]}

# -----------------------------
# Fixed parameters
# -----------------------------
SEEDSTART=0
NUMGRAPHS=100
FORMAT=hdf5
BACKEND=auto

# -----------------------------
# Run a single task by TASK_ID
# -----------------------------
run_task() {
    local TASK_ID=$1

    local NODE_INDEX=$(( TASK_ID / NUM_M ))
    local M_INDEX=$(( TASK_ID % NUM_M ))

    local NUMNODES=${NODES[$NODE_INDEX]}
    local EDGESPERNODE=${EDGES_PER_NODE[$M_INDEX]}

    # Sanity check: BA constraint 1 <= m < n
    if (( EDGESPERNODE >= NUMNODES )); then
        echo "Invalid BA parameters: m=${EDGESPERNODE} >= n=${NUMNODES}"
        return 1
    fi

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
}

# -----------------------------
# Determine task ID(s) to run
# -----------------------------
NUM_TASKS=$(( NUM_NODES * NUM_M ))

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # Running under SLURM job array
    run_task "${SLURM_ARRAY_TASK_ID}"
elif [[ -n "${1:-}" ]]; then
    # Task ID passed as command-line argument
    run_task "$1"
else
    # No SLURM, no argument: run all tasks sequentially
    echo "No SLURM_ARRAY_TASK_ID or argument provided. Running all ${NUM_TASKS} tasks sequentially."
    for (( i=0; i<NUM_TASKS; i++ )); do
        echo "========== Task ${i}/${NUM_TASKS} =========="
        run_task "$i"
    done
fi
