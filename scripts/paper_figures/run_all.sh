#!/bin/bash
# Run all paper figure scripts sequentially.
#
# Usage:
#   bash run_all.sh              # run all figures
#   bash run_all.sh 2 3 7        # run only figures 2, 3, and 7
#
# Output PNGs land in the repo-root plots/ directory (DrWatson plotsdir()).
# Each figure is run as a separate Julia process so that const declarations
# in one script don't conflict with another.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# The JuliaQAOA project is the repo root, two levels up from scripts/paper_figures.
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

declare -A FIGURES=(
    [2]="figure2_stddev_heatmap.jl"
    [3]="figure3_pearson_correlation.jl"
    [4]="figure4_squared_overlap.jl"
    [5]="figure5_objective_landscapes.jl"
    [6]="figure6_approx_ratio_comparison.jl"
    [7]="figure7_high_depth_performance.jl"
)

# Determine which figures to run
if [ $# -gt 0 ]; then
    TO_RUN=("$@")
else
    TO_RUN=(2 3 4 5 6 7)
fi

JULIA_THREADS="${JULIA_NUM_THREADS:-auto}"
JULIA_ARGS=(--project="$REPO_ROOT" --threads="$JULIA_THREADS")

echo "Julia args: ${JULIA_ARGS[*]}"
echo "Figures to run: ${TO_RUN[*]}"
echo ""

FAILED=()

for num in "${TO_RUN[@]}"; do
    script="${FIGURES[$num]:-}"
    if [ -z "$script" ]; then
        echo "WARNING: No figure $num — skipping"
        continue
    fi

    echo "================================================================"
    echo " Figure $num: $script"
    echo " Started: $(date)"
    echo "================================================================"

    if julia "${JULIA_ARGS[@]}" "$SCRIPT_DIR/$script"; then
        echo "Finished figure $num: $(date)"
    else
        echo "ERROR: figure $num failed (exit $?)"
        FAILED+=("$num")
    fi
    echo ""
done

echo "================================================================"
if [ ${#FAILED[@]} -eq 0 ]; then
    echo "All figures completed successfully."
    echo "Output: $REPO_ROOT/plots/"
else
    echo "Completed with errors. Failed figures: ${FAILED[*]}"
    exit 1
fi
