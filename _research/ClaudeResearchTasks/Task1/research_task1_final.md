# Research Task 1: Replicating Paper Figures — Final Report

Reproduction of Figures 2–7 from "Parameter-setting heuristic for the quantum
alternating operator ansatz" (Sud, Hadfield, Rieffel, Tubman, Hogg).

All 6 figure scripts are in `julia/paper_figures/`, use CairoMakie for plotting,
and output to `julia/paper_figures/output/`. Each script was run with small
parameters (n=6–8, 3 instances) to verify correctness. All produced valid output.

---

## Deliverables Summary

| Deliverable | Status | Location |
|---|---|---|
| Implementation plan | Done | `research_task1_plan.md` |
| Figure 2 script | Done | `julia/paper_figures/figure2_stddev_heatmap.jl` |
| Figure 3 script | Done | `julia/paper_figures/figure3_pearson_correlation.jl` |
| Figure 4 script | Done | `julia/paper_figures/figure4_squared_overlap.jl` |
| Figure 5 script | Done | `julia/paper_figures/figure5_objective_landscapes.jl` |
| Figure 6 script | Done | `julia/paper_figures/figure6_approx_ratio_comparison.jl` |
| Figure 7 script | Done | `julia/paper_figures/figure7_high_depth_performance.jl` |
| Shared utilities | Done | `julia/paper_figures/common.jl` |
| Linear ramp API | Done | `julia/src/linear_ramp.jl` (exported from JuliaQAOA) |
| CLAUDE.md updates | Done | `CLAUDE.md` § "Paper Figure Reproduction Scripts" |
| Per-figure report | Done | `julia/paper_figures/REPORT.md` |
| Output plots (15 PNGs) | Done | `julia/paper_figures/output/` |

---

## Shared Infrastructure

### `common.jl` — Shared Utilities

All figure scripts `include("common.jl")`, which provides:

- **Graph generation**: `erdos_renyi_edges(n, p; rng)`, `maxcut_costs(n, edges)`,
  `maxcut_optimal(costs)`, `generate_er_instance(n, p_edge; rng)`
- **Real QAOA simulator**: `apply_phase_gate!(state, costs, γ)`,
  `apply_x_mixer!(state, β, n)`, `qaoa_statevector(costs, n, γs, βs)`,
  `qaoa_expectation(costs, n, γs, βs)`
- **Plotting**: `save_figure(fig, name)`, `FIGURE_SIZE`, `HALF_FIGURE_SIZE`

The real QAOA simulator implements the FUR approach (qubit-by-qubit X-mixer
rotation), consistent with the paper's description and the existing Python
implementation in `qokit/fur/`.

### `julia/src/linear_ramp.jl` — General-Purpose API (New)

Added to the JuliaQAOA module (not script-local) per the Figure 7 requirement:

- `linear_ramp(γ₁, γ_f, β₁, β_f, p; pi_units=true)` — generates a single
  linear ramp schedule returning `(γs, βs)` vectors of length `p`
- `linear_ramp_matrix(γ₁s, γ_fs, β₁s, β_fs, p; pi_units=true)` — batch version
  generating `K×p` matrices for use with `QAOA_proxy_multi`

Both are exported from `JuliaQAOA` and available anywhere in the Julia codebase.

---

## Per-Figure Details

### Figure 2: Stddev/Mean Heatmap

**Script**: `julia/paper_figures/figure2_stddev_heatmap.jl`
**Timestamps**: Started 2026-03-17, Finished 2026-03-17
**Output**: `figure2_stddev_heatmap_cprime{0-9}.png` (10 files, ~60 KB each)

#### Approach
1. Generate `NUM_INSTANCES` Erdős-Rényi graphs G(n, p)
2. Compute N(c'; d, c) for each via `get_homogeneous_distribution_from_costs_direct`
3. Compute element-wise mean and stddev across instances using
   `distributions_mean_and_stddev`
4. Plot stddev/mean heatmap for each c' value using CairoMakie

#### Quick-Test Parameters vs Paper
| Parameter | Quick test | Paper |
|---|---|---|
| n (vertices) | 6 | 10 |
| p_edge | 0.5 | 1/3 |
| instances | 3 | 10 |
| c' plotted | all (0–9) | c'=7 |

#### Results
- Pattern matches paper: low coefficient of variation near center (d ≈ n/2,
  c ≈ m/2), higher at edges. This confirms the homogeneous approximation is
  best where the dominant terms of the amplitude sum are.
- Gray cells (NaN) correctly indicate N(c';d,c) = 0 regions (impossible
  combinations of distance and cost).
- All 10 c' values produced valid heatmaps; slices with all-zero data are
  automatically skipped.

#### Challenges
- Different random graph instances have different numbers of edges, requiring
  `max_num_edges` padding for array compatibility.
- CairoMakie throws errors when all heatmap values are NaN (high c' slices
  where no bitstrings exist). Fixed by adding a skip check: `if all(isnan, ratio)`.

#### Customization Knobs
- `N_QUBITS`, `P_EDGE`, `NUM_INSTANCES`, `SEED` — configuration constants at
  top of file
- `COST_PRIME_VALUES` — set to `:all` to save all c' values, or a specific
  vector like `[7]` to match the paper
- To swap graph types: replace the `generate_er_instance` call with a different
  generator (e.g., Barabási-Albert, Watts-Strogatz). The rest of the pipeline
  is graph-agnostic — it only needs `(edges, costs, num_edges, num_vertices)`.

#### Future Directions
- Run with paper parameters (n=10, p=1/3, 10 instances) for quantitative match
- Add non-ER graph generators to `common.jl`

---

### Figure 3: Pearson Correlation + Analytical Comparison

**Script**: `julia/paper_figures/figure3_pearson_correlation.jl`
**Timestamps**: Started 2026-03-17, Finished 2026-03-17
**Output**: `figure3_pearson_correlation.png` (202 KB)

#### Approach
1. Generate graphs, compute empirical averaged N(c';d,c) via
   `get_homogeneous_distribution_from_costs_direct` + `average_distributions`
2. For each proxy in `PROXY_CONFIGS`, compute analytical N(c';d,c) via
   `cpu_compute_homodist`
3. Compute Pearson correlation for each c' via
   `get_pearson_correlation_coefficients`
4. Main plot: correlation vs c' with P(c') bars in background
5. Insert heatmaps: side-by-side empirical vs analytical for selected c' values

#### Quick-Test Parameters vs Paper
| Parameter | Quick test | Paper |
|---|---|---|
| n | 6 | 10 |
| p_edge | 0.5 | 1/3 |
| instances | 3 | 10 |
| proxies | PaperProxy, TriangleProxy | Analytical formula |
| inserts | c'=3, c'=5 | c'=7, c'=13 |

#### Results
- PaperProxy shows high correlation (~1) for dominant terms (where P(c') is
  large), dropping off for tail terms — matching the paper's finding
- TriangleProxy (unfitted, default parameters) shows somewhat lower correlation,
  especially in the tails. This is expected since the default parameters haven't
  been fitted to the actual distribution.
- The heatmap inserts clearly show the visual similarity between empirical and
  analytical distributions for the dominant c' values.

#### Customization Knobs
- `PROXY_CONFIGS` — list of `(label, constructor)` pairs. All proxies are
  overlaid on the same correlation plot with different colors. Easy to add
  NormalProxy or a fitted TriangleProxy.
- `INSERT_COST_PRIMES` — which c' values get side-by-side heatmap subplots
- Graph parameters at top of file

#### Ambiguity
- The paper refers to the "analytical method" without specifying which proxy
  type (since only the PaperProxy formula existed at the time). The script
  treats any proxy's `N_cost_distance_distribution` as an "analytical method"
  to enable comparison.

---

### Figure 4: Squared Overlap vs QAOA Depth

**Script**: `julia/paper_figures/figure4_squared_overlap.jl`
**Timestamps**: Started 2026-03-17, Finished 2026-03-17
**Output**: `figure4_squared_overlap.png` (153 KB)

#### Approach
1. Generate one ER graph instance
2. For each (γ₁, γ_f) parameter pair:
   - Compute real QAOA state at each layer via `qaoa_statevector` with
     `return_intermediates=true`
   - Compute proxy state by running the compressed proxy via
     `QAOA_proxy_single`, then reconstructing the full 2^n state:
     `q(x) = Q(c(x))` for all x
   - Compute |⟨ψ_true|ψ_proxy⟩|² at each layer

The key insight: under the homogeneous approximation, the proxy amplitude for
bitstring x depends only on c(x). So we run the efficient compressed proxy
(which tracks only m+1 amplitudes) and reconstruct the full state by assigning
Q(c(x)) to each bitstring x.

#### Quick-Test Parameters vs Paper
| Parameter | Quick test | Paper |
|---|---|---|
| n | 6 | 8 |
| p (depth) | 8 | 20 |
| curves | 4 (γ₁, γ_f) pairs | multiple |
| β schedule | β₁=0.45, β_f=0.05 | not specified |

#### Results
- Qualitative match with paper: smaller γ values maintain higher overlap
  across all depths
- Overlap decreases gradually with depth but stays near 1.0 for small γ
  (the γ₁=0.05, γ_f=0.15 curve), confirming the proxy is faithful
- The curve ordering (smaller γ → better overlap) is consistent across runs,
  which is the paper's key observation

#### Challenges
- Required implementing a minimal QAOA statevector simulator in Julia
  (`common.jl`: `apply_phase_gate!`, `apply_x_mixer!`, `qaoa_statevector`).
  The X-mixer uses the qubit-by-qubit FUR approach consistent with the paper
  and existing Python code.
- The `proxy_statevector_from_compressed` function handles the reconstruction
  step, including clamping costs to the proxy's valid range.

#### Ambiguity
- **The paper does not specify the exact (γ₁, γ_f) and (β₁, β_f) values
  used in Figure 4.** We chose representative values that demonstrate the
  correct qualitative behavior. The script's `GAMMA_PAIRS` and β constants
  are easy to adjust.

#### Customization Knobs
- `GAMMA_PAIRS` — list of (γ₁, γ_f) pairs for different curves
- `β₁`, `β_f` — shared β schedule
- `PROXY_CONFIG` — which proxy to use for the homogeneous approximation
- `P_DEPTH` — maximum QAOA depth

---

### Figure 5: Objective Function Landscapes

**Script**: `julia/paper_figures/figure5_objective_landscapes.jl`
**Timestamps**: Started 2026-03-17, Finished 2026-03-17
**Output**: `figure5_objective_landscapes.png` (72 KB)

#### Approach
1. Generate one graph instance
2. Fix all but the last layer's parameters (γ₁, …, γ_{p-1}, β₁, …, β_{p-1})
3. Sweep the last layer's (γ_p, β_p) on a grid
4. For the true landscape: evaluate `qaoa_expectation` at each grid point
5. For the proxy landscape: batch-evaluate using `QAOA_proxy_multi` (all grid
   points in one call)
6. Plot side-by-side heatmaps

#### Quick-Test Parameters vs Paper
| Parameter | Quick test | Paper |
|---|---|---|
| n | 6 | 8 |
| p | 2 | 3 |
| grid | 20×20 | 30×30 |
| fixed params | arbitrary | from p=20 optimization |

#### Results
- True and proxy landscapes show similar peak locations and overall structure
- The proxy landscape has a more compressed dynamic range (smoother valleys),
  consistent with the proxy's averaging effect
- Peak locations are approximately aligned, confirming the proxy is useful for
  parameter optimization

#### Challenges
- The paper fixes γ₁,γ₂,β₁,β₂ from a p=20 optimization on a 20-node graph.
  For our quick test, we used arbitrary fixed values. The landscape structure
  is qualitatively similar regardless — the key point is that the proxy and
  true landscapes have correlated features.
- `QAOA_proxy_multi` was essential for efficiency: evaluating 400 grid points
  in a single matrix multiplication rather than 400 separate proxy runs.

#### Ambiguity
- **The paper doesn't specify the exact fixed parameter values** (they come
  from a separate optimization on a larger graph). Our quick-test uses
  `FIXED_GAMMAS = [0.2]`, `FIXED_BETAS = [0.4]`. For a faithful reproduction,
  one would first run a larger optimization to obtain these fixed values.

#### Customization Knobs
- `PROXY_CONFIGS` — list of proxies; each gets its own landscape panel
- `FIXED_GAMMAS`, `FIXED_BETAS` — fixed parameters for early layers
- `GAMMA_RANGE`, `BETA_RANGE` — sweep ranges for the last layer
- `GRID_SIZE` — grid resolution

---

### Figure 6: Approximation Ratio Comparison (Low Depth)

**Script**: `julia/paper_figures/figure6_approx_ratio_comparison.jl`
**Timestamps**: Started 2026-03-17, Finished 2026-03-17
**Output**: `figure6_approx_ratio_comparison.png` (119 KB)

#### Approach
1. **Parameter Transfer**: Optimize QAOA on source graphs (coordinate descent
   with random restarts), take element-wise median of optimal parameters,
   evaluate on target graphs
2. **Homogeneous Heuristic**: For each proxy in `PROXY_CONFIGS` and each target
   graph, optimize the proxy objective (grid search for p=1, random sampling
   for p>1), then evaluate the proxy-optimal parameters on real QAOA
3. Box plot comparison: one panel per depth p, boxes for Transfer and each proxy

#### Quick-Test Parameters vs Paper
| Parameter | Quick test | Paper |
|---|---|---|
| n_source | 6 | 9 |
| n_target | 8 | 20 |
| p values | [1, 2] | [1, 2, 3] |
| source instances | 3 | 10 |
| target instances | 3 | 10 |
| optimizer | coord. descent + random sampling | COBYLA |

#### Results
- Transfer and PaperProxy perform comparably (within error bars at p=1),
  matching the paper's finding that the heuristic is competitive with transfer
- Unfitted TriangleProxy underperforms both, especially at p=2 — expected since
  the default parameters haven't been fitted to the actual distribution
- Even with only 3 instances, the relative ordering of methods is clear

#### Challenges
- Real QAOA optimization on source graphs uses coordinate descent rather than
  COBYLA (which would require Optim.jl or a scipy call). For quick-test sizes,
  this gives reasonable results but is not as thorough as a proper optimizer.
- Proxy optimization for p>1 uses random sampling in the 2p-dimensional
  parameter space (grid search over 2p dimensions is infeasible). The paper
  uses BFGS which would be more efficient.

#### Customization Knobs
- `PROXY_CONFIGS` — add any proxy type; each appears as its own box in the
  comparison. Pre-configured with PaperProxy and TriangleProxy.
- `N_SOURCE`, `N_TARGET`, `P_VALUES`, `NUM_SOURCE_INSTANCES`,
  `NUM_TARGET_INSTANCES` — all configurable at top of file
- `N_RESTARTS` — number of random restarts for source graph optimization

#### Future Directions
- Use Optim.jl for proper BFGS/NelderMead optimization instead of grid/random
  search
- Increase instance counts (10+) for statistical significance
- Test with paper's original sizes: n_source=9, n_target=20
- Add fitted TriangleProxy and NormalProxy for a fairer comparison

---

### Figure 7: High Depth Performance with Linear Ramp

**Script**: `julia/paper_figures/figure7_high_depth_performance.jl`
**Timestamps**: Started 2026-03-17, Finished 2026-03-17
**Output**: `figure7_high_depth_performance.png` (69 KB)

#### Approach
1. For each proxy type and each depth p:
   - Optimize the 4 linear ramp parameters (γ₁, γ_f, β₁, β_f) by exhaustive
     grid search over the proxy objective. Grid of size `GRID_SIZE_PER_DIM^4`
     evaluated in a single `QAOA_proxy_multi` call.
   - Convert optimized parameters to actual γ, β schedules via `linear_ramp`
   - Evaluate via real QAOA and compute approximation ratio
2. Grouped box plot: proxy types × depths

#### Quick-Test Parameters vs Paper
| Parameter | Quick test | Paper |
|---|---|---|
| n | 8 | 20 |
| p values | [2, 4, 6] | [4, 8, 12, 16, 20] |
| instances | 3 | 10 |
| grid per dim | 8 (8^4 = 4096 pts) | BFGS optimizer |

#### Results
- PaperProxy shows monotonic improvement with depth (higher p → higher
  approximation ratio), matching the paper's key finding
- Unfitted TriangleProxy degrades with depth — its inaccurate N(c';d,c)
  leads the optimizer to choose parameters that work well for the proxy
  landscape but poorly for real QAOA
- The gap between proxies widens with depth, highlighting the importance of
  accurate distribution approximation for high-depth QAOA

#### Linear Ramp API
The `linear_ramp` and `linear_ramp_matrix` functions were added to
`julia/src/linear_ramp.jl` and exported from the `JuliaQAOA` module. This
satisfies the task requirement that linear ramp scheduling be a general API,
not confined to this script. The API supports:
- Single parameter set generation: `linear_ramp(γ₁, γ_f, β₁, β_f, p)`
- Batch generation for use with `QAOA_proxy_multi`:
  `linear_ramp_matrix(γ₁s, γ_fs, β₁s, β_fs, p)`

#### Challenges
- Grid search over 4 parameters (8^4 = 4096 evaluations) is expensive but
  feasible thanks to the `QAOA_proxy_multi` batch API, which evaluates all
  parameter sets in one matrix multiplication per layer.
- For production-scale runs (n=20, p=20), BFGS via Optim.jl would be more
  efficient than exhaustive grid search. The grid approach suffices for
  demonstration and correctness verification.

#### Customization Knobs
- `PROXY_CONFIGS` — list of proxies to compare. Each gets its own color in
  the grouped box plot.
- `P_VALUES` — depths to test
- `GRID_SIZE_PER_DIM` — controls optimization resolution (total points =
  GRID_SIZE_PER_DIM^4)
- All graph parameters configurable at top of file

#### Future Directions
- Integrate Optim.jl for gradient-based optimization of linear ramp parameters
- Test fitted TriangleProxy/NormalProxy (fitted to averaged N from many
  instances via `sendai_opt.fit_proxy_to_real` or a Julia equivalent)
- Scale to n=20 with paper parameters for direct quantitative comparison

---

## Verification: Output Files

All 15 output files exist in `julia/paper_figures/output/`:

```
figure2_stddev_heatmap_cprime0.png   (60 KB)
figure2_stddev_heatmap_cprime1.png   (62 KB)
figure2_stddev_heatmap_cprime2.png   (60 KB)
figure2_stddev_heatmap_cprime3.png   (60 KB)
figure2_stddev_heatmap_cprime4.png   (60 KB)
figure2_stddev_heatmap_cprime5.png   (63 KB)
figure2_stddev_heatmap_cprime6.png   (63 KB)
figure2_stddev_heatmap_cprime7.png   (63 KB)
figure2_stddev_heatmap_cprime8.png   (61 KB)
figure2_stddev_heatmap_cprime9.png   (61 KB)
figure3_pearson_correlation.png     (202 KB)
figure4_squared_overlap.png         (153 KB)
figure5_objective_landscapes.png     (72 KB)
figure6_approx_ratio_comparison.png (119 KB)
figure7_high_depth_performance.png   (69 KB)
```

All files dated 2026-03-17 21:48–21:59.

---

## Architecture of the Script Suite

```
julia/
  src/
    JuliaQAOA.jl                    Module orchestrator (exports all)
    linear_ramp.jl                  NEW: linear_ramp, linear_ramp_matrix
    QAOA_proxy.jl                   QAOA_proxy_basic/single/multi, expectation
    cost_distributions.jl           get_homogeneous_distribution_from_costs_direct
    paper_proxy.jl / triangle_proxy.jl / normal_proxy.jl
    utils.jl                        cpu_compute_homodist, allocate_homodist

  paper_figures/
    common.jl                       Graph gen, real QAOA simulator, plotting
    figure2_stddev_heatmap.jl       Stddev/mean heatmaps of N(c';d,c)
    figure3_pearson_correlation.jl  Pearson correlation + heatmap inserts
    figure4_squared_overlap.jl      |⟨ψ_true|ψ_proxy⟩|² vs depth
    figure5_objective_landscapes.jl True vs proxy objective landscapes
    figure6_approx_ratio_comparison.jl  Transfer vs proxy box plots
    figure7_high_depth_performance.jl   Approx ratio vs depth (linear ramp)
    REPORT.md                       Detailed per-figure report
    output/                         15 generated PNG files
```

### Design Pattern

All scripts follow the same pattern:
1. `include("common.jl")` — loads JuliaQAOA module, CairoMakie, shared utilities
2. **CONFIGURATION** block — all tunable constants at top (graph params, proxy
   choice, grid sizes, etc.)
3. **COMPUTATION** block — generate graphs, compute distributions, run proxy/QAOA
4. **PLOTTING** block — create CairoMakie figures, save to output/

Comments are written in a literate programming style: each major step has a
plain-English comment explaining what it does and why.

---

## CLAUDE.md Updates

The following section was added to `CLAUDE.md` under "Paper Figure Reproduction
Scripts":

- Summary of each figure script's purpose, location, and key configuration
  constants
- Description of `common.jl` shared infrastructure
- Description of the `linear_ramp.jl` API addition
- Cross-references to JuliaQAOA module functions used by each script

---

## Known Limitations and Potential Improvements

1. **Optimizer quality**: Figures 6–7 use coordinate descent / grid search
   rather than proper BFGS (which requires Optim.jl). Results are qualitatively
   correct but quantitatively suboptimal.

2. **Scale**: Quick-test parameters (n=6–8, 3 instances) produce correct trends
   but not publication-quality statistics. Paper parameters (n=10–20, 10+
   instances) require more compute time.

3. **Unfitted proxies**: TriangleProxy and NormalProxy are used with default
   (unfitted) parameters in Figures 3, 6, 7. After fitting (via
   `sendai_opt.fit_proxy_to_real` or a Julia equivalent), these proxies would
   perform much better.

4. **Non-ER graphs**: All scripts currently generate Erdős-Rényi graphs.
   Swapping to other families (Barabási-Albert, Watts-Strogatz) requires only
   replacing the `generate_er_instance` call and adding a graph generator to
   `common.jl`.

5. **Exact figure match**: Some paper parameters are not fully specified
   (notably the fixed parameters in Figure 5 and the γ/β values in Figure 4).
   Our choices produce the correct qualitative behavior.
