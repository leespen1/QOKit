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
# -----------------------------
NODES=(10 11 12 13 14 15 16 17 18 19 20)
KVALUES=(4 6 8)
BETAS=(0.05 0.10 0.20 0.30)

NUM_NODES=${#NODES[@]}
NUM_K=${#KVALUES[@]}
NUM_BETA=${#BETAS[@]}

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

    local NODE_INDEX=$(( TASK_ID / (NUM_K * NUM_BETA) ))
    local REM=$(( TASK_ID % (NUM_K * NUM_BETA) ))
    local K_INDEX=$(( REM / NUM_BETA ))
    local BETA_INDEX=$(( REM % NUM_BETA ))

    local NUMNODES=${NODES[$NODE_INDEX]}
    local K=${KVALUES[$K_INDEX]}
    local BETA=${BETAS[$BETA_INDEX]}

    # Sanity checks (Watts–Strogatz constraints)
    if (( K >= NUMNODES )); then
        echo "Invalid WS parameters: k=${K} >= n=${NUMNODES}"
        return 1
    fi

    if (( K % 2 != 0 )); then
        echo "Invalid WS parameters: k=${K} must be even"
        return 1
    fi

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
}

# -----------------------------
# Determine task ID(s) to run
# -----------------------------
NUM_TASKS=$(( NUM_NODES * NUM_K * NUM_BETA ))

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
