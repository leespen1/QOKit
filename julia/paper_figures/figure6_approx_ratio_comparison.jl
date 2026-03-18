#=
figure6_approx_ratio_comparison.jl — Reproduce Figure 6 from the paper.

Paper Figure 6: Box plot comparing approximation ratios from two approaches:
  1. "Parameter Transfer": Optimize QAOA on small (source) graphs, transfer
     the median parameters to larger (target) graphs
  2. "Homogeneous Heuristic": Optimize proxy parameters directly for each
     target graph

Paper parameters: Source G(9, 1/2), Target G(20, 1/2), p=1,2,3, 10 instances
Quick test:       Source G(6, 0.5), Target G(8, 0.5), p=1,2, 3 instances

Customization:
  - Change PROXY_CONFIGS to compare multiple proxies
  - Adjust source/target graph sizes, p values, number of instances

Started:  2026-03-17
Finished: 2026-03-17
=#

include("common.jl")
using Statistics: median, mean

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_SOURCE = 6           # Source graph size (paper: 9)
const N_TARGET = 8           # Target graph size (paper: 20)
const P_EDGE = 0.5
const NUM_SOURCE_INSTANCES = 3   # (paper: 10)
const NUM_TARGET_INSTANCES = 3   # (paper: 10)
const P_VALUES = [1, 2]         # Depths to test (paper: [1, 2, 3])
const SEED = 42

# Number of random restarts for QAOA optimization on source graphs
const N_RESTARTS = 5

# Proxies to compare (all used for the "Homogeneous Heuristic" approach)
PROXY_CONFIGS = [
    ("PaperProxy", (m, n) -> PaperProxy(m, n, P_EDGE)),
    ("TriangleProxy", (m, n) -> OldTriangleProxy(m, n)),
]

#==============================================================================#
#                    OPTIMIZATION HELPERS                                       #
#==============================================================================#

"""
    optimize_qaoa_random_restarts(costs, n, p, n_restarts; rng)

Find good QAOA parameters by random restart optimization (COBYLA-style
grid search). Returns the best (γs, βs) found and the corresponding ⟨C⟩.
"""
function optimize_qaoa_random_restarts(
    costs::Vector{Float64}, n::Int, p::Int, n_restarts::Int;
    rng=Random.default_rng()
)
    best_exp = -Inf
    best_γs = zeros(p)
    best_βs = zeros(p)

    for _ in 1:n_restarts
        # Random initial parameters in reasonable ranges (QOKit convention)
        γs = rand(rng, p) .* 1.6
        βs = rand(rng, p) .* (π/2)

        # Simple grid refinement around initial point
        current_γs = copy(γs)
        current_βs = copy(βs)
        current_exp = qaoa_expectation(costs, n, current_γs, current_βs)

        # Coordinate descent with shrinking step sizes
        for step_scale in [0.2, 0.1, 0.05, 0.02]
            for param_idx in 1:(2p)
                best_local = current_exp
                best_val = param_idx <= p ? current_γs[param_idx] : current_βs[param_idx - p]

                for delta in [-2, -1, -0.5, 0.5, 1, 2] .* step_scale
                    trial_γs = copy(current_γs)
                    trial_βs = copy(current_βs)
                    if param_idx <= p
                        trial_γs[param_idx] = max(0, best_val + delta)
                    else
                        trial_βs[param_idx - p] = clamp(best_val + delta, 0, π/2)
                    end
                    trial_exp = qaoa_expectation(costs, n, trial_γs, trial_βs)
                    if trial_exp > best_local
                        best_local = trial_exp
                        if param_idx <= p
                            current_γs[param_idx] = trial_γs[param_idx]
                        else
                            current_βs[param_idx - p] = trial_βs[param_idx - p]
                        end
                    end
                end
                current_exp = best_local
            end
        end

        if current_exp > best_exp
            best_exp = current_exp
            best_γs = copy(current_γs)
            best_βs = copy(current_βs)
        end
    end

    return best_γs, best_βs, best_exp
end

"""
    optimize_proxy_grid(proxy, homodist, P_vals, n, p, grid_size)

Optimize proxy parameters by grid search over γ, β.
Returns the best (γs, βs) in radians and the proxy expectation.
"""
function optimize_proxy_grid(
    homodist::AbstractArray{<:Real, 3},
    P_vals::Vector{Float64},
    n::Int, p::Int;
    grid_size::Int=30
)
    m = size(homodist, 1) - 1

    if p == 1
        # Direct 2D grid search for p=1 (QOKit convention)
        γ_range = range(0.02, 1.6, length=grid_size)
        β_range = range(0.01, π/2 - 0.01, length=grid_size)

        # Build parameter matrices for batch evaluation
        K = grid_size^2
        γ_matrix = zeros(K, 1)
        β_matrix = zeros(K, 1)
        idx = 0
        for γ in γ_range, β in β_range
            idx += 1
            γ_matrix[idx, 1] = γ / π  # Convert to pi_units
            β_matrix[idx, 1] = β / π
        end

        Qs = QAOA_proxy_multi(homodist, γ_matrix, β_matrix; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n))

        best_idx = argmax(exps)
        best_γ_pi = γ_matrix[best_idx, 1]
        best_β_pi = β_matrix[best_idx, 1]

        return [best_γ_pi * π], [best_β_pi * π], exps[best_idx]

    else
        # For p>1, use random sampling in the proxy parameter space (QOKit convention)
        K = grid_size^2
        γ_matrix = rand(K, p) .* 1.6 ./ π  # In pi_units
        β_matrix = rand(K, p) .* 0.5       # β ∈ [0, 0.5] in pi_units

        Qs = QAOA_proxy_multi(homodist, γ_matrix, β_matrix; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n))

        best_idx = argmax(exps)
        return γ_matrix[best_idx, :] .* π, β_matrix[best_idx, :] .* π, exps[best_idx]
    end
end

#==============================================================================#
#                          COMPUTATION                                          #
#==============================================================================#

println("=== Figure 6: Approximation Ratio Comparison ===")
println("Source: G($N_SOURCE, $P_EDGE) × $NUM_SOURCE_INSTANCES")
println("Target: G($N_TARGET, $P_EDGE) × $NUM_TARGET_INSTANCES")
println("p values: $P_VALUES")

rng = MersenneTwister(SEED)

# Generate source and target graphs
source_instances = [generate_er_instance(N_SOURCE, P_EDGE; rng) for _ in 1:NUM_SOURCE_INSTANCES]
target_instances = [generate_er_instance(N_TARGET, P_EDGE; rng) for _ in 1:NUM_TARGET_INSTANCES]

# Find optimal costs for target graphs (brute force)
target_optimal = [maxcut_optimal(inst.costs) for inst in target_instances]
println("Target optimal costs: $target_optimal")

# --- Results storage ---
# results[p_idx][method_name] = vector of approximation ratios
all_results = Dict{Int, Dict{String, Vector{Float64}}}()

for (p_idx, p) in enumerate(P_VALUES)
    println("\n--- p = $p ---")
    all_results[p] = Dict{String, Vector{Float64}}()

    # === Method 1: Parameter Transfer ===
    println("  Optimizing source graphs (transfer method)...")
    source_params = map(source_instances) do inst
        γs, βs, exp_val = optimize_qaoa_random_restarts(inst.costs, N_SOURCE, p, N_RESTARTS; rng)
        (γs=γs, βs=βs, exp=exp_val)
    end

    # Take element-wise median of source parameters
    median_γs = [median([sp.γs[ℓ] for sp in source_params]) for ℓ in 1:p]
    median_βs = [median([sp.βs[ℓ] for sp in source_params]) for ℓ in 1:p]
    println("  Median γs: $median_γs")
    println("  Median βs: $median_βs")

    # Evaluate transfer parameters on target graphs
    transfer_ratios = map(enumerate(target_instances)) do (i, inst)
        exp_val = qaoa_expectation(inst.costs, N_TARGET, median_γs, median_βs)
        exp_val / target_optimal[i]
    end
    all_results[p]["Transfer"] = transfer_ratios
    println("  Transfer ratios: $transfer_ratios")

    # === Method 2: Homogeneous Heuristic (for each proxy) ===
    for (label, constructor) in PROXY_CONFIGS
        println("  Optimizing via $label proxy...")
        proxy_ratios = map(enumerate(target_instances)) do (i, inst)
            m_target = inst.num_edges
            proxy = constructor(m_target, N_TARGET)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:m_target]

            # Optimize via proxy
            best_γs, best_βs, proxy_exp = optimize_proxy_grid(
                homodist, P_vals, N_TARGET, p; grid_size=30
            )

            # Evaluate on real QAOA
            real_exp = qaoa_expectation(inst.costs, N_TARGET, best_γs, best_βs)
            real_exp / target_optimal[i]
        end
        all_results[p][label] = proxy_ratios
        println("  $label ratios: $proxy_ratios")
    end
end

#==============================================================================#
#                          PLOTTING                                             #
#==============================================================================#

println("\nPlotting...")
n_methods = 1 + length(PROXY_CONFIGS)  # Transfer + proxies
method_names = vcat(["Transfer"], [cfg[1] for cfg in PROXY_CONFIGS])
colors = [:steelblue, :coral, :mediumseagreen, :orange, :purple]

fig = Figure(size=(300 * length(P_VALUES) + 100, 500))

for (p_idx, p) in enumerate(P_VALUES)
    ax = Axis(fig[1, p_idx],
        xlabel="Method",
        ylabel="Approximation Ratio",
        title="p = $p",
        xticks=(1:n_methods, method_names),
        xticklabelrotation=π/6,
    )

    for (m_idx, method) in enumerate(method_names)
        ratios = all_results[p][method]
        # Box plot via scatter + error bars
        μ = mean(ratios)
        boxplot!(ax, fill(m_idx, length(ratios)), ratios,
            color=colors[mod1(m_idx, length(colors))],
            width=0.6,
        )
    end
end

Label(fig[0, :],
    "Approximation Ratio: Transfer vs Proxy Heuristic\nSource G($N_SOURCE, $P_EDGE) → Target G($N_TARGET, $P_EDGE)",
    fontsize=14, font=:bold)

save_figure(fig, "figure6_approx_ratio_comparison.png")
println("Done!")
