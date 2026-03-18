# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

QOKit (Quantum Optimization Toolkit) is a Python library for simulating and
benchmarking the Quantum Approximate Optimization Algorithm (QAOA). The `grips/`
directory contains G-RIPS 2024 Sendai research on improving QAOA
parameter-setting heuristics based on the paper "Parameter-setting heuristic for
the quantum alternating operator ansatz." The users of Claude Code in this
project are the Sendai researchers.

**Branches**: Main branch is `grips` (stable). Active development on `gpu` branch.

**Virtual Environment**: At least on Spencer's computer, a virtual environment
(for the python portion of the project) is provided in `../qokitvenv`, which has
the necessary python packages installed.

## Common Commands

```bash
# Install (development mode)
pip install -e .

# Install with GPU support (CUDA 12.x)
pip install -e .[GPU-CUDA12]

# Python-only install (if C compilation fails)
QOKIT_PYTHON_ONLY=1 pip install -e .

# Run all tests with coverage
pytest --cov=qokit --cov-fail-under=75 -rs tests

# Run a single test file
pytest tests/test_qaoa_objective_maxcut.py

# Run GRIPS-specific tests
pytest grips_tests/

# Check formatting
black --check .

# Format code
black .
```

## Architecture

### Simulation Pipeline

The core QAOA simulation uses Fast Unitary Rotation (FUR) algorithm with multiple backends:
- **GPU (CUDA + Numba)**: Fastest, requires `cupy` - `qokit/fur/nbcuda/`
- **C-compiled**: Fast CPU - `qokit/fur/c/`
- **Python fallback**: Reference implementation - `qokit/fur/python/`

`qokit.fur.choose_simulator()` auto-selects the fastest available backend.

### GRIPS Proxy System

The proxy system approximates QAOA state evolution without full quantum simulation:

1. **Real Distribution** (`grips/real_distribution.py`): Computes `n(x; d, c)` - count of bitstrings at Hamming distance `d` with cost `c` from bitstring `x`

2. **Homogeneous Distribution**: Cost-averaged `N(c'; d, c)` enables parameter prediction across graph instances

3. **Proxy Classes** (all in `grips/`):
   - `PaperProxy` (`paper_proxy.py`): Original paper's binomial/multinomial approach
   - `NormalProxy` (`normal_proxy.py`): Multivariate normal approximation
   - `TriangleProxy` (`triangle_proxy.py`): GRIPS contribution - simplified triangle distribution

4. **Interface** (`grips/QAOA_proxy_interface.py`): Unified proxy API with Python (Numba JIT) and Julia backends

### Key Module Relationships

```
grips/QAOA_simulator.py     - Main simulation interface (QAOA_run, get_simulator)
grips/QAOA_proxy_interface.py - Proxy algorithm entry point
  └── paper_proxy.py / normal_proxy.py / triangle_proxy.py
grips/real_distribution.py  - Statistical distribution computation
grips/sendai_opt.py        - Parameter optimization (fit_proxy_to_real)
qokit/fur/__init__.py      - Simulator backend selection
```

### Julia Backend

Julia implementations in `julia/src/` provide high-performance alternatives. Configure with `USE_JULIA=True` in proxy files. Setup via `grips/setup_juliacall.py`.

## Code Style

- **Formatter**: Black with `line-length=160`
- **Python**: 3.10, 3.11
- **Type checking**: pyright
- **License headers**: Apache 2.0, SPDX headers required (checked by `addheader`)

## Key Concepts

- **QAOA**: Quantum Approximate Optimization Algorithm - hybrid quantum-classical algorithm
- **MaxCut**: Graph partitioning problem, primary benchmark in this codebase
- **LABS**: Low Autocorrelation Binary Sequences problem
- **Hamming Distance**: Number of differing bits between two bitstrings
- **FUR**: Fast Unitary Rotation - efficient state evolution algorithm


---

## The Parameter-Setting Paper: Concepts and Code Mapping

The paper "A Parameter Setting Heuristic for the Quantum Alternating Operator
Ansatz" (Sud, Hadfield, Rieffel, Tubman, Hogg) is the theoretical foundation for
the `grips/` code. The full LaTeX source is in
`References/ParameterSettingHeuristicLatexSource/main.tex`. Below is a plain-English
summary of the paper's key ideas, each linked to the specific functions that implement
them.

---

### Idea 1: Perfect Homogeneity (paper §2)

**Plain English**: In a normal QAOA run, each of the 2^n bitstrings gets its own
complex probability amplitude. "Perfect Homogeneity" is the assumption that all
bitstrings with the *same cost value* always have *exactly the same amplitude*. This
means instead of tracking 2^n amplitudes, you only need to track one amplitude per
unique cost value. For MaxCut on m-edge graphs, there are at most m+1 unique costs,
so this can reduce the state description from exponentially large to polynomially
large.

Perfect Homogeneity holds exactly for some symmetric problems (e.g. the Hamming ramp
c(x)~|x|) and is approximately true empirically for random instances of many CSPs,
especially when the QAOA parameters γ are small.

**Code**: This is the core assumption underlying all proxy classes. No single function
implements it, but it motivates the data structure used by every proxy: instead of a
2^n state vector, the proxy tracks a vector of length `num_constraints + 1` (one
amplitude per cost value).

---

### Idea 2: The n(x; d, c) Distribution (paper §2)

**Plain English**: For any bitstring x, `n(x; d, c)` counts how many bitstrings are
simultaneously (a) at Hamming distance d from x, and (b) have cost c. This is the
key quantity needed to propagate QAOA amplitudes: when applying one QAOA layer to
bitstring x, the new amplitude is a sum over all other bitstrings y, weighted by the
mixing operator matrix element (which depends only on Hamming distance d(x,y)) and
the phase factor (which depends only on cost c(y)).

**Code**:
- `grips/real_distribution.py: get_real_distribution()` — computes n(x; d, c) for
  every bitstring x in a graph, as a 3D array indexed by `[x, d, c]`. This is the
  *exact* distribution, computed by brute force in O(2^(2n)) time.
- `grips/real_distribution.py: get_real_distribution_from_costs()` — the core Numba
  JIT-compiled loop. Takes a pre-computed cost array and fills the 3D array.
- `julia/src/cost_distributions.jl: get_real_distribution_from_costs()` — Julia
  equivalent, with multithreading (`@threads`) for speedup.
- `julia/src/cost_distributions.jl: gpu_get_real_distribution_from_costs()` — GPU
  (CUDA) version using a kernel where each thread handles one (x, y) pair.

---

### Idea 3: The Homogeneous Distribution N(c'; d, c) (paper §2)

**Plain English**: Under Perfect Homogeneity, we replace the per-bitstring n(x; d, c)
with a per-*cost* version N(c'; d, c). For a given source cost c', N(c'; d, c) is the
*average* n(x; d, c) over all bitstrings x with cost c(x) = c'. If this average is
a good representative for all individual x, then the proxy faithfully approximates
QAOA. The paper shows empirically (for Erdős-Rényi MaxCut) that the variance of
n(x; d, c) around its cost-class mean is small for the dominant terms of the sum.

There are two ways to compute N(c'; d, c):

**Method A — Empirical averaging** (exact for given instances, O(2^(2n))):
Average n(x; d, c) over all bitstrings x with the same cost. Expensive but exact.
- `grips/real_distribution.py: get_homogeneous_distribution()` — top-level function
  accepting a single graph or list of graphs; returns the averaged N(c'; d, c).
- `grips/real_distribution.py: get_homogeneous_distribution_from_costs()` — low-level
  Numba JIT loop that does the averaging given a pre-computed real distribution.
- `julia/src/cost_distributions.jl: get_homogeneous_distribution_from_costs()` — Julia
  multithreaded version.
- `julia/src/cost_distributions.jl: get_homogeneous_distribution_from_costs_direct()` —
  more memory-efficient Julia version that skips storing the full n(x; d, c).
- `julia/src/cost_distributions.jl: gpu_get_homogeneous_distribution_from_costs_direct()` —
  GPU version using privatized global memory to reduce atomic contention.

**Method B — Analytical formula** (class-level only, O(poly(n)), paper §3):
For random CSPs, N(c'; d, c) can be computed analytically without any graph instances,
using a multinomial formula that depends only on the problem class (n, m, edge
probability). The paper derives this for MaxCut, MaxE3Lin2/Max-k-XOR, and Rand-k-SAT
(paper Eq. 12-16). Implemented in the Proxy classes described below.

---

### Idea 4: Proxy Classes — Implementations of N(c'; d, c) (paper §3)

Each proxy class implements three methods:
- `P_cost_distribution(c')` — probability that a random bitstring has cost c'
- `N_cost_distribution(c')` — expected number of bitstrings with cost c' (= 2^n × P(c'))
- `N_cost_distance_distribution(c', d, c)` — the key N(c'; d, c) value

**PaperProxy** (`grips/paper_proxy.py`, `julia/src/paper_proxy.jl`):
- Directly implements the paper's analytical formula (paper §3.1, Eq. 11-16) for
  MaxCut on Erdős-Rényi graphs.
- `P_cost_distribution`: Binomial distribution (paper Eq. 11).
- `N_cost_distance_distribution`: Multinomial sum (paper Eq. 12-16, involving
  P_both, P_one, P_neither probabilities).
- Parameters: `num_constraints` (number of edges m), `num_qubits` (number of vertices n),
  `prob_edge` (edge probability p_e, default 0.5).
- This is the "ground truth" analytical proxy from the paper; no fitting required.

**TriangleProxy** (`grips/triangle_proxy.py`, `julia/src/triangle_proxy.jl`):
- **GRIPS contribution**: Approximates N(c'; d, c) as a "prism" shape: for each fixed
  distance d, the distribution over costs c is a triangle (piecewise linear function
  with one peak). The peak location slides linearly from c' (at d=0) to m/2 (at
  d=n/2), and the peak height scales with distance.
- Much faster to evaluate than PaperProxy (no multinomial sums).
- Has 4 tunable parameters: `h_tweak_sub` (peak height offset), `hc_tweak_add` (peak
  cost offset), `l_tweak_mul` (left slope), `r_tweak_mul` (right slope).
- `HardCodedTriangleProxy`: Obsolete version with hardcoded default parameters.
- `TriangleProxy` is Numba `@jitclass` compiled for maximum speed.
- Parameters must be fit to real N(c'; d, c) data using `sendai_opt.fit_proxy_to_real()`.

**NormalProxy** (`grips/normal_proxy.py`, `julia/src/normal_proxy.jl`):
- **GRIPS contribution**: Approximates N(c'; d, c) as a 2D multivariate normal
  distribution over (cost c, distance d).
- The covariance matrix is constructed from 3 parameters: `cost_mean`, `cov_1`,
  `cov_2`, rotated by the deviation of c' from the mean.
- Parameters must be fit to real N(c'; d, c) data using `sendai_opt.fit_proxy_to_real()`.

---

### Idea 5: Running the Proxy — Algorithm 1 (paper §2.3)

**Plain English**: Given the N(c'; d, c) distributions, QAOA parameters (γ, β), and
an initial uniform state (all costs get amplitude 1/√2^n), iterate p times:

  Q_l(c') = Σ_{d,c} cos(β)^(n-d) × (-i·sin(β))^d × exp(-iγc) × Q_{l-1}(c) × N(c'; d, c)

This is Eq. (8) in the paper. Instead of the full 2^n quantum state vector, we only
work with a vector of `m+1` complex amplitudes (one per unique cost). Time complexity
O(n × m² × p) vs O(2^n × p) for full simulation.

**Code**:
- `grips/QAOA_proxy_interface.py: QAOA_proxy()` — dispatcher that routes to the
  appropriate backend (Numba JIT, Python, or Julia) based on proxy type.
- `grips/QAOA_proxy_interface.py: compute_amplitude_sum()` — computes one step of the
  inner loop (one new Q_l(c') value). Implements the triple sum over d and c.
- `grips/QAOA_proxy_interface.py: compute_amplitude_sum_njit()` — Numba JIT version.
- `grips/QAOA_proxy_interface.py: QAOA_proxy_njit()` — full Numba JIT proxy run.
- `julia/src/QAOA_proxy.jl: QAOA_proxy_basic()` — readable Julia reference implementation.
- `julia/src/QAOA_proxy.jl: QAOA_proxy_single()` — faster Julia version using BLAS
  matrix-vector multiplication (reshapes the sum into a mat-vec product).
- `julia/src/QAOA_proxy.jl: QAOA_proxy_multi()` — batch Julia version for evaluating
  multiple (γ, β) parameter sets simultaneously via BLAS matrix-matrix multiplication.
- `julia/src/QAOA_proxy.jl: get_β_factors()` / `get_γ_factors()` — helper functions
  that precompute the cos/sin and exp factors for all distances/costs at once.

---

### Idea 6: Computing the Expected Cost from the Proxy (paper §2.2, Eq. 9)

**Plain English**: After running the proxy to get amplitudes Q_p(c'), compute the
expected cost as:

  ⟨C⟩ ≈ Σ_{c'} 2^n × P(c') × |Q_p(c')|² × c'

This is the "homogeneous parameter objective function" — the quantity we maximize
when optimizing γ and β. It replaces the expensive exact quantum expectation value.

**Code**:
- `grips/QAOA_proxy_interface.py: QAOA_proxy_expectation()` — computes this sum given
  the final proxy amplitudes. Dispatches to Numba or Julia versions.
- `grips/QAOA_proxy_interface.py: QAOA_proxy_expectation_njit()` — Numba JIT version.
- `grips/QAOA_proxy_interface.py: QAOA_proxy_expectation_from_gamma_beta()` —
  convenience function that runs the proxy and returns the expectation in one call.
- `julia/src/QAOA_proxy.jl: expectation()` — Julia version. Has two overloads: one
  for a single state vector Q, and one for a matrix Q where each column is a state
  (used with `QAOA_proxy_multi`).

---

### Idea 7: Homogeneous Heuristic for Parameter Setting — Algorithm 3 (paper §4)

**Plain English**: The full parameter-setting strategy:
1. Precompute N(c'; d, c) and P(c') for the problem class.
2. Choose initial parameters (γ_in, β_in), e.g. a linear ramp schedule.
3. Run a classical optimizer (e.g. COBYLA, BFGS) that calls the proxy to evaluate
   the objective function for each candidate (γ, β).
4. The optimizer outputs (γ_out, β_out) which are then used in real QAOA.

Key insight: each proxy evaluation is O(n × m² × p) instead of O(2^n × p), so
high-depth optimization (e.g. p=20) becomes tractable on a laptop.

**Code**:
- `grips/QAOA_proxy_interface.py: QAOA_proxy_optimize_gamma_beta()` — main
  optimization loop. Takes a proxy, initial parameters, and an optimizer method
  (default COBYLA). Returns optimized γ, β, and diagnostics. Dispatches to Numba or
  Julia backends.
- `grips/QAOA_proxy_interface.py: inverse_proxy_objective_function()` — wraps the
  proxy evaluation as a minimization objective (negates because scipy minimizes).

---

### Idea 8: Fitting Proxy Parameters to Real Distributions (GRIPS-specific)

**Plain English**: The TriangleProxy and NormalProxy have free parameters that need
to be tuned so their N(c'; d, c) shape matches the empirically observed distribution.
This is done by minimizing MSE between the proxy's predicted N(c'; d, c) and the
real averaged N(c'; d, c) computed from actual graph instances.

**Code**:
- `grips/sendai_opt.py: fit_proxy_to_real()` — main fitting function. Uses "smart
  random search": perturb randomly, reuse helpful perturbations, shrink step sizes
  after consecutive failures. Optional grid search for initial parameters.
- `grips/real_distribution.py: distribution_mean_squared_error()` — computes MSE
  between a proxy's predicted N(c'; d, c) and a real distribution array.
- `grips/real_distribution.py: get_homogeneous_distribution_from_proxy()` — evaluates
  a proxy's N(c'; d, c) for all cost/distance/cost combinations, returning a 3D array
  for comparison with the real distribution.
- `grips/real_distribution.py: normalize_homodist_slices()` — normalizes each N(c'; :, :)
  slice independently (useful for fair comparison across proxies with different scales).

---

### Idea 9: Validating the Proxy (paper §5)

**Plain English**: The paper provides three empirical checks that the proxy is good:
1. The standard deviation of n(x; d, c) across different x with same cost c' is small
   (relative to the mean) for the dominant terms.
2. The analytically-computed N(c'; d, c) correlates well (Pearson correlation ~1) with
   the empirical average over graph instances, for dominant terms.
3. The squared overlap between the proxy state and the true QAOA state stays high
   across p layers (especially for small γ values and linear ramp schedules).

**Code**:
- `grips/real_distribution.py: distributions_mean_and_stddev()` — computes mean and
  stddev across multiple distribution arrays (used for check 1).
- `grips/real_distribution.py: plot_stddev_div_mean_heatmap()` — plots the coefficient
  of variation heatmap shown in the paper.
- `grips/real_distribution.py: get_pearson_correlation_coefficients()` — computes
  Pearson correlation between two N(c'; d, c) arrays, one per cost c' (used for check 2).
- `julia/src/cost_distributions.jl: get_pearson_correlation_coefficients()` — Julia version.

---

### Idea 10: Cost Distribution P(c') Estimation (paper §2.2, §3; QAOA_proxy_interface.py)

**Plain English**: Beyond the analytical binomial formula in PaperProxy, the GRIPS
team added advanced methods for estimating P(c') from graph instances (useful when
the analytical formula is unavailable or we want instance-specific estimates).

**Code** (in `grips/QAOA_proxy_interface.py`):
- `CostDistribution` class — callable wrapper for P(c') supporting multiple backends:
  `wang_landau` (discrete MCMC-based sampling), `moment_matching` with sub-types
  `gaussian`, `beta`, or `edgeworth` expansion.
- `estimate_P_for_graphs()` — estimates P(c') from a list of graph instances using
  either Monte Carlo sampling + moment matching, or Wang-Landau sampling.
- `_run_wang_landau()` — Wang-Landau MCMC sampler for estimating density of states.

---

## GRIPS Research Goals and Future Directions

The following directions are mentioned in the paper (§7 Discussion) and represent
active or potential research areas for the Sendai team:

1. **Instance-specific N(c'; d, c)**: The paper only uses class-level distributions
   (derived analytically from just n, m, and edge probability, not the specific graph).
   Extending to per-instance distributions (e.g. by Monte Carlo sampling of
   n(x; d, c) for a specific graph) could yield better parameters for individual
   instances. The `get_homogeneous_distribution()` function already supports this
   for small graphs.

2. **New graph types**: The data_generation/ directory includes scripts for
   Barabasi-Albert and Watts-Strogatz graphs in addition to Erdős-Rényi. The
   analytical PaperProxy formulas apply specifically to Erdős-Rényi; other graph
   families require either empirical distributions or new analytical derivations.

3. **Triangle/Normal proxy vs Paper proxy**: The GRIPS triangle and normal proxies
   are faster to evaluate than the paper's multinomial approach. The key research
   question is: after fitting, do these proxies produce parameter schedules that
   perform as well on real QAOA as the PaperProxy? `distribution_mean_squared_error()`
   and `get_pearson_correlation_coefficients()` support this comparison.

4. **Scaling to larger n and p**: The paper demonstrates p=20 on n=20 graphs using
   linear ramp schedules. The Julia backend (`QAOA_proxy_multi`) supports batched
   evaluation of many parameter sets simultaneously, enabling more thorough landscape
   exploration. GPU acceleration in `cost_distributions.jl` enables computing
   N(c'; d, c) for larger graphs.

5. **Linear ramp schedules**: The paper restricts to linear ramps at high depth to
   reduce the parameter space from 2p to 4 parameters. The schedule is:
   γ_ℓ = γ_1 + (γ_f - γ_1) × ℓ/p, same for β. This avoids the curse of dimensionality
   at large p.

6. **Approximation ratio**: The standard evaluation metric. For MaxCut:
   ApxRatio = ⟨C⟩ / c_opt, where c_opt is found by brute force. The `grips/solve_maxcut_exact.py`
   module computes c_opt.

---

## Experimental Findings: P(c') Investigation

The script `grips_examples/investigate_P_distribution.py` tested whether replacing
the paper's analytical P(c') = Binomial(m, 0.5) with the real empirical P(c') from
specific graph instances improves QAOA proxy parameter setting.

### Experiment Setup

Three proxy configurations were compared on G(10, 0.5) Erdős-Rényi MaxCut at p=1:
1. **Paper N + Paper P**: Analytical N(c';d,c) from PaperProxy + Binomial P(c')
2. **Paper N + Real P**: Analytical N(c';d,c) + empirical P(c') from graph instance
3. **Real N + Real P**: Instance-specific N(c';d,c) computed via
   `get_homogeneous_distribution_from_costs_direct()` + empirical P(c')

For each configuration, a (γ, β) grid was swept using `QAOA_proxy_multi()` and
`expectation()`, the proxy-optimal parameters were found, and real QAOA was evaluated
at those parameters via `get_expectation()`.

### Key Results

- **Paper N + Paper P wins consistently** (100/100 seeds over Paper N + Real P;
  10/10 seeds over Real N + Real P).
- Paper N + Real P performs worst (mean approximation ratio deficit ~0.025 vs Paper).
- Real N + Real P is closer to Paper (mean deficit ~0.009) but still worse.

### Interpretation

The analytical proxy's N and P are derived from the same probabilistic model
(Binomial/Multinomial over independent edges), making them **mutually consistent**.
The Q amplitudes computed from analytical N are calibrated to the analytical P's
assumptions. Swapping in real P breaks this consistency. Furthermore, the analytical
formulas act as a **regularizer**: they smooth over instance-specific noise in the
cost landscape, producing a proxy landscape whose optima generalize better.

This means instance-specific distributions are unlikely to beat the paper's approach
for Erdős-Rényi graphs, where the analytical formula is available and accurate.

---

## Current Research Directions

Based on the P(c') investigation findings, the following directions have been
identified for improving on the paper's results. **Spencer is currently focusing
on direction 2.**

### Direction 1: Better Fitted Proxy Shapes

Fit TriangleProxy or NormalProxy to N(c';d,c) data averaged over *many* graph
instances (not just one). Multi-instance averaging provides the same regularization
benefit as the analytical formula. The fitted proxy could then be paired with the
analytical Binomial P(c') for consistency. This leverages the speed advantage of
TriangleProxy while potentially matching PaperProxy's accuracy.

### Direction 2: Non-Erdős-Rényi Graph Families (Current Focus)

**This is the most promising direction.** PaperProxy's analytical formula is derived
specifically for Erdős-Rényi random graphs. For other graph families (Barabási-Albert,
Watts-Strogatz, real-world networks), no analytical formula exists, so the PaperProxy
cannot be used at all. This makes fitted proxies (TriangleProxy, NormalProxy) the
*only* option for these graph types. The data_generation/ directory already includes
scripts for generating Barabási-Albert and Watts-Strogatz graphs. The key question
is whether fitted proxies can achieve good approximation ratios on these non-ER
graph families where PaperProxy is unavailable.

### Direction 3: Higher Depth with Linear Ramp Schedules

The P(c') investigation used p=1. At higher depths (p=5, 10, 20), the proxy
landscape becomes more structured and the gap between proxy and real QAOA may
change. Linear ramp schedules (4 parameters: γ_1, γ_f, β_1, β_f) reduce the
search space and are already supported by the proxy infrastructure.

### Direction 4: Hybrid Proxy Warmstart + Local Refinement

Use the proxy to find a good starting point (γ*, β*), then run a few iterations of
real QAOA optimization (COBYLA/BFGS) starting from that point. This combines the
proxy's cheap global search with real QAOA's accuracy for local refinement.

### Direction 5: Multi-Instance Averaged Homodist

Compute empirical N(c';d,c) from multiple graph instances of the same class,
averaged together. This provides class-level smoothing similar to the analytical
formula but works for any graph family. Combine with a consistent P(c') estimated
from the same set of instances (e.g., via moment matching or Wang-Landau).

---

## Paper Figure Reproduction Scripts

All scripts are in `julia/paper_figures/` with output in `julia/paper_figures/output/`.
Full report: `julia/paper_figures/REPORT.md`. Each script uses CairoMakie for plotting
and the JuliaQAOA module for proxy computations. Configuration constants at the top of
each script make it easy to change graph type, size, proxy, etc.

### Shared Infrastructure

- **`julia/paper_figures/common.jl`**: ER graph generation (`erdos_renyi_edges`),
  MaxCut cost computation (`maxcut_costs`), real QAOA statevector simulator
  (`qaoa_statevector`, `qaoa_expectation`), and plotting utilities.

- **`julia/src/linear_ramp.jl`** (added to JuliaQAOA module): General-purpose
  linear ramp schedule API. `linear_ramp(γ₁, γ_f, β₁, β_f, p)` generates
  parameter vectors; `linear_ramp_matrix(...)` generates batch parameter matrices
  compatible with `QAOA_proxy_multi`.

### Per-Figure Scripts

- **`figure2_stddev_heatmap.jl`**: Stddev/mean heatmap of N(c';d,c) across graph
  instances. Generates one heatmap per c' value. Configurable: graph params,
  instances, which c' to plot.

- **`figure3_pearson_correlation.jl`**: Pearson correlation between empirical and
  proxy N(c';d,c). Main plot with P(c') overlay + insert heatmaps for selected c'.
  Supports comparing multiple proxies on the same plot via `PROXY_CONFIGS`.

- **`figure4_squared_overlap.jl`**: |⟨ψ_true|ψ_proxy⟩|² vs QAOA layer for linear
  ramp schedules. Multiple (γ₁, γ_f) curves. Uses compressed proxy + state
  reconstruction (q(x) = Q(c(x))).

- **`figure5_objective_landscapes.jl`**: Side-by-side heatmaps of true vs proxy
  objective function. Fixes all but last layer's parameters, sweeps on grid.
  Uses `QAOA_proxy_multi` for batch proxy evaluation.

- **`figure6_approx_ratio_comparison.jl`**: Box plots comparing parameter transfer
  (from small source graphs) vs proxy heuristic on target graphs. Supports
  multiple proxy types. Uses coordinate descent for source optimization.

- **`figure7_high_depth_performance.jl`**: Box plots of approximation ratio vs
  QAOA depth using linear ramp schedules optimized via proxy. Demonstrates
  monotonic improvement with depth. Supports multiple proxy types.
