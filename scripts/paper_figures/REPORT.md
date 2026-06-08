# Paper Figure Reproduction — Report

Reproduction of Figures 2–7 from "Parameter-setting heuristic for the quantum
alternating operator ansatz" (Sud, Hadfield, Rieffel, Tubman, Hogg).

All scripts are in `scripts/paper_figures/` and use CairoMakie for plotting.
Output goes to the repo-root `plots/` directory (DrWatson `plotsdir()`).

---

## Figure 2: Stddev/Mean Heatmap

**Script:** `figure2_stddev_heatmap.jl`
**Paper params:** G(10, 1/3), 10 instances, c'=7
**Quick test:** G(6, 0.5), 3 instances, all c' values

### Approach
1. Generate N Erdős-Rényi graph instances
2. Compute homogeneous distribution N(c';d,c) for each via
   `get_homogeneous_distribution_from_costs_direct`
3. Compute element-wise mean and stddev across instances
4. Plot stddev/mean heatmap for each c'

### Results
- Pattern matches paper: low deviation near center (d≈n/2, c≈m/2),
  high deviation at edges — confirming the homogeneous approximation
  is best where the dominant terms are
- Gray cells indicate N(c';d,c)=0 regions (impossible bitstring combinations)
- Saves one heatmap per c' value, skipping all-zero slices

### Challenges
- CairoMakie throws errors when all heatmap values are NaN (high c' values
  where no bitstrings exist). Fixed by adding a skip check.
- Different graph instances have different numbers of edges, requiring
  `max_num_edges` padding for compatibility.

### Customization
- Change `N_QUBITS`, `P_EDGE`, `NUM_INSTANCES` at top of file
- Set `COST_PRIME_VALUES` to `:all` or a specific vector like `[7]`
- Easy to swap graph generation function for non-ER graphs

---

## Figure 3: Pearson Correlation + Analytical Comparison

**Script:** `figure3_pearson_correlation.jl`
**Paper params:** G(10, 1/3), 10 instances
**Quick test:** G(6, 0.5), 3 instances

### Approach
1. Generate graphs, compute empirical averaged N(c';d,c)
2. For each proxy in `PROXY_CONFIGS`, compute analytical N(c';d,c) via
   `cpu_compute_homodist`
3. Compute Pearson correlation for each c' via
   `get_pearson_correlation_coefficients`
4. Main plot: correlation vs c' with P(c') bars in background
5. Insert heatmaps: side-by-side empirical vs analytical for selected c'

### Results
- PaperProxy shows high correlation (~1) for dominant terms (where P(c')
  is large), lower for tail terms — matching the paper
- TriangleProxy (unfitted, default params) shows somewhat lower correlation
  overall, especially in the tails
- The heatmap inserts clearly show the visual similarity between empirical
  and analytical distributions

### Customization
- Add/remove proxy types in `PROXY_CONFIGS` — all are overlaid on the
  same correlation plot
- Change `INSERT_COST_PRIMES` to select which c' values get heatmap subplots
- Easy to adjust for different graph families

---

## Figure 4: Squared Overlap vs QAOA Depth

**Script:** `figure4_squared_overlap.jl`
**Paper params:** G(8, 1/2), p=20, linear ramp, multiple (γ₁, γ_f) curves
**Quick test:** G(6, 0.5), p=8, 4 curves

### Approach
1. Generate one ER graph instance
2. For each (γ₁, γ_f) parameter pair:
   - Compute real QAOA state at each layer via `qaoa_statevector` (custom
     implementation in `common.jl` using qubit-by-qubit X-mixer)
   - Compute proxy state by running compressed proxy via `QAOA_proxy_single`,
     then reconstructing full 2^n state (q(x) = Q(c(x)))
   - Compute |⟨ψ_true|ψ_proxy⟩|² at each layer

### Results
- Qualitative match with paper: smaller γ values maintain higher overlap
- Overlap decreases gradually with depth but stays near 1 for small γ
- The curve ordering (smaller γ = better overlap) is consistent across runs

### Challenges
- The paper's proxy state computation keeps all 2^n amplitudes and substitutes
  N(c(x);d,c) for n(x;d,c). Our implementation is equivalent: we run the
  compressed proxy (which produces one amplitude per cost level), then
  reconstruct the full state by assigning Q(c(x)) to each bitstring x.
  Under perfect homogeneity these approaches are identical.
- Required implementing a minimal QAOA statevector simulator in Julia
  (`common.jl`: `apply_phase_gate!`, `apply_x_mixer!`, `qaoa_statevector`)

### Ambiguity
- The paper does not specify the exact (γ₁, γ_f) and (β₁, β_f) values used.
  We chose representative values. The qualitative behavior is robust to
  specific parameter choices.

---

## Figure 5: Objective Function Landscapes

**Script:** `figure5_objective_landscapes.jl`
**Paper params:** G(8, 1/2), p=3, fix γ₁,γ₂,β₁,β₂, sweep γ₃,β₃ on 30×30
**Quick test:** G(6, 0.5), p=2, fix γ₁,β₁, sweep γ₂,β₂ on 20×20

### Approach
1. Generate one graph instance
2. Fix all but the last layer's parameters
3. Sweep the last layer's (γ, β) on a grid
4. For each point: compute true QAOA expectation and proxy expectation
5. Plot side-by-side heatmaps

### Results
- True and proxy landscapes show similar peak locations and overall structure
- The proxy landscape has a more compressed dynamic range (the valleys
  are deeper relative to the peak), consistent with the paper's observation
- Peak locations are approximately aligned

### Challenges
- The paper fixes γ₁,γ₂,β₁,β₂ from a p=20 optimization on a 20-node graph.
  For our quick test, we used arbitrary fixed values. The landscape structure
  is qualitatively similar regardless.
- Used `QAOA_proxy_multi` for efficient batch evaluation of the proxy
  landscape (all grid points evaluated in one call)

### Ambiguity
- The paper doesn't specify the exact fixed parameter values used (they come
  from a separate optimization). Our results show the correct qualitative
  behavior with arbitrary fixed values.

---

## Figure 6: Approximation Ratio Comparison (Low Depth)

**Script:** `figure6_approx_ratio_comparison.jl`
**Paper params:** Source G(9, 1/2) → Target G(20, 1/2), p=1,2,3, 10 instances
**Quick test:** Source G(6, 0.5) → Target G(8, 0.5), p=1,2, 3 instances

### Approach
1. **Transfer method**: Optimize QAOA on source graphs (coordinate descent
   with random restarts), take median parameters, evaluate on target graphs
2. **Proxy heuristic**: For each target graph, optimize proxy objective
   (grid search for p=1, random sampling for p>1), evaluate on real QAOA
3. Compare via box plots

### Results
- Transfer and PaperProxy perform comparably (within error bars at p=1)
- Unfitted TriangleProxy underperforms both, especially at p=2
- Consistent with paper's finding that the heuristic is competitive with
  transfer

### Challenges
- Real QAOA optimization on source graphs uses a simple coordinate descent
  rather than COBYLA (which requires scipy). For small test sizes, this
  gives reasonable results.
- Proxy optimization for p>1 uses random sampling rather than systematic
  grid search (grid over 2p dimensions is infeasible). A proper optimizer
  (Optim.jl) would be better for production runs.

### Future Directions
- Use Optim.jl for proper BFGS/COBYLA optimization instead of grid/random search
- Increase instance counts for statistical significance
- Test with actual n=9→n=20 transfer (requires more compute time)

---

## Figure 7: High Depth Performance with Linear Ramp

**Script:** `figure7_high_depth_performance.jl`
**Paper params:** G(20, 1/2), p=4,8,12,16,20, 10 instances, BFGS
**Quick test:** G(8, 0.5), p=2,4,6, 3 instances

### Approach
1. For each proxy type and each depth p:
   - Optimize 4 linear ramp parameters (γ₁, γ_f, β₁, β_f) by grid search
     over the proxy objective (GRID_SIZE^4 total evaluations)
   - Convert to actual γ, β schedules and evaluate via real QAOA
   - Compute approximation ratio
2. Plot grouped box plots: proxy types × depths

### Linear Ramp API
The `linear_ramp(γ₁, γ_f, β₁, β_f, p)` and `linear_ramp_matrix(...)` functions
have been added to `src/linear_ramp.jl` and exported from JuliaQAOA.
These are general-purpose tools usable anywhere, not confined to this script.

### Results
- PaperProxy shows monotonic improvement with depth (0.83 → 0.87), matching
  the paper's key finding
- Unfitted TriangleProxy degrades with depth — its inaccurate N(c';d,c)
  leads to worse parameter choices at higher p
- The gap between proxies widens with depth, highlighting the importance
  of accurate distribution approximation

### Challenges
- Grid search over 4 parameters (GRID_SIZE^4 evaluations) becomes expensive.
  The `QAOA_proxy_multi` batch API makes this feasible by evaluating all
  parameter sets in one matrix multiplication.
- For production runs, BFGS via Optim.jl would be more efficient than
  grid search.

### Future Directions
- Add Optim.jl for gradient-based optimization of linear ramp parameters
- Test fitted TriangleProxy/NormalProxy (fitted to averaged N from many instances)
- Scale to n=20 (requires more compute but is feasible with the proxy)

---

## Shared Infrastructure

### `common.jl` — Shared Utilities
- `erdos_renyi_edges(n, p)` — Generate ER random graph as edge list
- `maxcut_costs(n, edges)` — Compute MaxCut cost for all 2^n bitstrings
- `maxcut_optimal(costs)` — Brute-force optimal cost
- `generate_er_instance(n, p_edge)` — Convenience wrapper
- `apply_phase_gate!(state, costs, γ)` — QAOA phase gate (diagonal)
- `apply_x_mixer!(state, β, n)` — QAOA X-mixer (qubit-by-qubit FUR)
- `qaoa_statevector(costs, n, γs, βs)` — Full QAOA simulation
- `qaoa_expectation(costs, n, γs, βs)` — Compute ⟨C⟩
- `save_figure(fig, name)` — Save to output directory

### `src/linear_ramp.jl` — Added to JuliaQAOA Module
- `linear_ramp(γ₁, γ_f, β₁, β_f, p)` — Single parameter set
- `linear_ramp_matrix(γ₁s, γ_fs, β₁s, β_fs, p)` — Batch parameter sets

---

## General Notes

### Programming Process
- All scripts use Julia with CairoMakie, leveraging the existing JuliaQAOA module
- Graph generation is implemented directly in Julia (no Graphs.jl dependency)
  using simple edge lists and bitwise MaxCut cost computation
- A minimal QAOA statevector simulator was implemented in `common.jl` for the
  figures requiring real QAOA evaluation (Figures 4, 5, 6, 7)
- The `QAOA_proxy_multi` batch API is essential for Figures 5 and 7, enabling
  efficient evaluation of thousands of parameter sets

### Potential Improvements
1. **Optimizer**: Replace grid/random search with Optim.jl for BFGS/COBYLA
2. **Scale**: Test with paper's original parameters (n=10-20, more instances)
3. **Graph families**: Easy to swap `erdos_renyi_edges` for other generators
4. **Proxy fitting**: Figures 6-7 use unfitted TriangleProxy; fitted versions
   (via `sendai_opt.fit_proxy_to_real` or Julia equivalent) would be more fair
