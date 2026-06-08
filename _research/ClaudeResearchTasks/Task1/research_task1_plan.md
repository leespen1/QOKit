# Research Task 1: Replicating Paper Figures — Plan

## Overview
Reproduce Figures 2–7 from "Parameter-setting heuristic for the quantum alternating operator ansatz" using Julia (CairoMakie for plotting) and existing QOKit/GRIPS infrastructure.

Each script will be placed in `julia/paper_figures/` and will:
- Be self-contained Julia scripts using the JuliaQAOA module
- Use CairoMakie for plotting
- Be heavily commented in a literate programming style
- Support easy customization (graph type, proxy type, parameters)
- Include timestamps

## Figure-by-Figure Plan

### Figure 2: Stddev/Mean Heatmap (paper §5.1)
**Paper params:** G(10, 1/3), 10 instances, c'=7
**Script:** `julia/paper_figures/figure2_stddev_heatmap.jl`
**Approach:**
1. Generate N random ER graphs using Graphs.jl
2. Compute real distribution n(x;d,c) for each via `get_real_distribution_from_costs`
3. Compute homogeneous distribution N(c';d,c) for each
4. Compute mean and stddev across instances
5. Plot stddev/mean heatmap for chosen c' using CairoMakie
6. Support saving all c' values at once
**Quick test:** n=6, 3 instances

### Figure 3: Pearson Correlation + Analytical Comparison (paper §5.1)
**Paper params:** G(10, 1/3), 10 instances
**Script:** `julia/paper_figures/figure3_pearson_correlation.jl`
**Approach:**
1. Generate graphs, compute averaged empirical N(c';d,c)
2. Compute analytical N(c';d,c) from PaperProxy
3. Compute Pearson correlation for each c'
4. Main plot: correlation coefficients vs c', with P(c') overlay
5. Subplots: heatmaps comparing analytical vs empirical for selected c' values
6. Support swapping proxy type for the "analytical" method
**Quick test:** n=6, 3 instances

### Figure 4: Squared Overlap vs QAOA Depth (paper §5.2)
**Paper params:** G(8, 1/2), p=20, linear ramp schedule
**Script:** `julia/paper_figures/figure4_squared_overlap.jl`
**Approach:**
1. Generate one ER graph
2. Compute costs for the graph
3. For each (γ₁, γ_f) pair with fixed β schedule:
   - Run real QAOA layer by layer (full 2^n statevector)
   - Run proxy QAOA layer by layer (using N(c(x);d,c) substitution for all 2^n amplitudes)
   - Compute squared overlap at each layer
4. Plot overlap vs layer for each parameter curve
**Challenge:** Need to implement real QAOA statevector evolution in Julia, or use the Python backend.
The real QAOA uses the FUR simulator from qokit. We can call Python from Julia, or implement the simple evolution directly in Julia since it's just applying diagonal phase gates and X-mixer.
**Decision:** Implement a minimal QAOA statevector simulator in Julia for this figure. The evolution is straightforward:
- Phase gate: multiply each amplitude by exp(-iγc(x))
- X-mixer: apply product of single-qubit X rotations = tensor product of [[cos β, -i sin β], [-i sin β, cos β]]
**Quick test:** n=6, p=5, 2 curves

### Figure 5: Objective Function Landscapes (paper §5.3)
**Paper params:** G(8, 1/2), p=3, sweep γ₃ and β₃ on 30×30 grid
**Script:** `julia/paper_figures/figure5_objective_landscapes.jl`
**Approach:**
1. Generate one ER graph
2. Fix γ₁, γ₂, β₁, β₂ (from paper or random optimization)
3. Sweep γ₃, β₃ on a grid
4. For each (γ₃, β₃):
   - Compute true expectation via full QAOA simulation
   - Compute proxy expectation via proxy
5. Plot side-by-side heatmaps
**Quick test:** n=6, p=2, 15×15 grid

### Figure 6: Approximation Ratio Comparison (paper §6.1)
**Paper params:** Transfer G(9,1/2) → G(20,1/2), p=1,2,3, 10 instances each
**Script:** `julia/paper_figures/figure6_approx_ratio_comparison.jl`
**Approach:**
1. Generate source graphs (small n), optimize QAOA parameters on each
2. Transfer: take median parameters, apply to target graphs
3. Homogeneous: optimize proxy parameters for each target graph
4. Evaluate both via real QAOA on target graphs
5. Compare via box plots
6. Support multiple proxy types
**Challenge:** Real QAOA on n=20 is 2^20 ≈ 1M amplitudes — feasible but slow for many evaluations. For quick test, use smaller n.
**Quick test:** n_source=6, n_target=8, p=1,2, 3 instances

### Figure 7: High Depth Performance (paper §6.2)
**Paper params:** G(20, 1/2), p=4,8,12,16,20, 10 instances, linear ramp, BFGS
**Script:** `julia/paper_figures/figure7_high_depth_performance.jl`
**Approach:**
1. Add linear ramp schedule to JuliaQAOA API (julia/src/linear_ramp.jl)
2. For each p value and graph instance:
   - Optimize 4 linear ramp parameters via proxy
   - Evaluate real QAOA at optimized parameters
3. Box plot of approximation ratios vs p
4. Support multiple proxy types
**Linear ramp API addition:** `linear_ramp(γ₁, γ_f, β₁, β_f, p)` → (γs, βs)
**Quick test:** n=8, p=2,4, 3 instances

## Shared Infrastructure Needs
1. **Graph generation in Julia:** Use Graphs.jl for ER graphs + MaxCut cost computation
2. **Real QAOA simulator in Julia:** Minimal statevector simulator for Figures 4, 5, 6, 7
3. **Linear ramp schedule:** Add to JuliaQAOA module
4. **CairoMakie plotting utilities:** Shared color schemes, layout helpers

## File Structure
```
julia/
  src/
    linear_ramp.jl          # New: linear ramp schedule API
    qaoa_simulation.jl      # New: minimal real QAOA simulator
  paper_figures/
    common.jl               # Shared utilities (graph gen, costs, etc.)
    figure2_stddev_heatmap.jl
    figure3_pearson_correlation.jl
    figure4_squared_overlap.jl
    figure5_objective_landscapes.jl
    figure6_approx_ratio_comparison.jl
    figure7_high_depth_performance.jl
```
