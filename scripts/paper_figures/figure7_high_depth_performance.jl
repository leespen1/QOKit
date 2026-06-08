#=
figure7_high_depth_performance.jl — Reproduce Figure 7 from the paper.

Paper Figure 7: Box plot of approximation ratios vs QAOA depth, using
linear ramp schedules (4 free parameters: γ₁, γ_f, β₁, β_f) optimized
via the homogeneous proxy.

Paper parameters: G(20, 1/2), p=4,8,12,16,20, 10 instances, BFGS optimizer
Quick test:       G(8, 0.5), p=2,4,6, 3 instances, grid search

The linear ramp schedule API (linear_ramp, linear_ramp_matrix) has been
added to the JuliaQAOA module (src/linear_ramp.jl) as a general-
purpose tool, not confined to this script.

Customization:
  - Change PROXY_CONFIGS to compare multiple proxies
  - Adjust graph size, p values, number of instances
  - Adjust GRID_SIZE for proxy optimization resolution

Started:  2026-03-17
Finished: 2026-03-17
=#

include("common.jl")
using Statistics: mean

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_QUBITS = 20
const P_EDGE = 0.5
const NUM_INSTANCES = 10
const P_VALUES = [4, 8, 12, 16, 20]
const SEED = 42

# Grid resolution for proxy optimization over (γ₁, γ_f, β₁, β_f) space
const GRID_SIZE_PER_DIM = 8     # Total grid = GRID_SIZE_PER_DIM^4

# Proxies to compare
PROXY_CONFIGS = [
    ("PaperProxy", (m, n) -> PaperProxy(m, n, P_EDGE)),
    ("TriangleProxy", (m, n) -> OldTriangleProxy(m, n)),
]

#==============================================================================#
#                    LINEAR RAMP PROXY OPTIMIZATION                             #
#==============================================================================#

"""
    optimize_linear_ramp_proxy(homodist, P_vals, n, p; grid_size)

Optimize the 4 linear ramp parameters (γ₁, γ_f, β₁, β_f) by grid search
over the proxy objective function.

Returns best (γ₁, γ_f, β₁, β_f) in pi_units and the proxy expectation.
"""
function optimize_linear_ramp_proxy(
    homodist::AbstractArray{<:Real, 3},
    P_vals::Vector{Float64},
    n::Int, p::Int;
    grid_size::Int=8
)
    # Parameter ranges (in pi_units, since proxy uses pi_units; QOKit convention)
    # γ: [0, 0.6] in pi_units ≈ [0, ~2 radians]
    # β: [0, 0.5] in pi_units = [0, π/2 radians]
    γ₁_range = range(0.02, 0.30, length=grid_size)
    γ_f_range = range(0.10, 0.60, length=grid_size)
    β₁_range = range(0.10, 0.45, length=grid_size)
    β_f_range = range(0.01, 0.20, length=grid_size)

    # Build all parameter combinations
    K = grid_size^4
    γ_matrix = zeros(K, p)
    β_matrix = zeros(K, p)

    idx = 0
    for γ₁ in γ₁_range, γ_f in γ_f_range, β₁ in β₁_range, β_f in β_f_range
        idx += 1
        γs, βs = linear_ramp(γ₁, γ_f, β₁, β_f, p)
        γ_matrix[idx, :] .= γs
        β_matrix[idx, :] .= βs
    end

    # Batch evaluate all parameter sets
    Qs = QAOA_proxy_multi(homodist, γ_matrix, β_matrix; pi_units=true)
    exps = vec(expectation(Qs[end], P_vals, n))

    # Find best
    best_idx = argmax(exps)
    best_exp = exps[best_idx]

    # Recover the linear ramp parameters
    γ₁_idx = ((best_idx - 1) ÷ (grid_size^3)) + 1
    rem = (best_idx - 1) % (grid_size^3)
    γf_idx = (rem ÷ (grid_size^2)) + 1
    rem = rem % (grid_size^2)
    β₁_idx = (rem ÷ grid_size) + 1
    βf_idx = (rem % grid_size) + 1

    best_γ₁ = γ₁_range[γ₁_idx]
    best_γ_f = γ_f_range[γf_idx]
    best_β₁ = β₁_range[β₁_idx]
    best_β_f = β_f_range[βf_idx]

    return (γ₁=best_γ₁, γ_f=best_γ_f, β₁=best_β₁, β_f=best_β_f, proxy_exp=best_exp)
end

#==============================================================================#
#                          COMPUTATION                                          #
#==============================================================================#

println("=== Figure 7: High Depth Performance ===")
println("Parameters: n=$N_QUBITS, p_edge=$P_EDGE, instances=$NUM_INSTANCES")
println("p values: $P_VALUES")
println("Grid: $(GRID_SIZE_PER_DIM)^4 = $(GRID_SIZE_PER_DIM^4) parameter sets per proxy/instance")

rng = MersenneTwister(SEED)
instances = [generate_er_instance(N_QUBITS, P_EDGE; rng) for _ in 1:NUM_INSTANCES]
optimal_costs = [maxcut_optimal(inst.costs) for inst in instances]
println("Optimal costs: $optimal_costs")

# Results: results[proxy_label][p] = vector of approx ratios
all_results = Dict{String, Dict{Int, Vector{Float64}}}()

for (label, constructor) in PROXY_CONFIGS
    println("\n--- Proxy: $label ---")
    all_results[label] = Dict{Int, Vector{Float64}}()

    for p in P_VALUES
        println("  p = $p:")
        ratios = Float64[]

        for (i, inst) in enumerate(instances)
            # Build proxy and homodist
            m = inst.num_edges
            proxy = constructor(m, N_QUBITS)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:m]

            # Optimize linear ramp parameters via proxy
            result = optimize_linear_ramp_proxy(
                homodist, P_vals, N_QUBITS, p;
                grid_size=GRID_SIZE_PER_DIM
            )

            # Generate the actual parameter vectors (convert from pi_units to radians)
            γs, βs = linear_ramp(result.γ₁, result.γ_f, result.β₁, result.β_f, p)
            γs_rad = γs .* π
            βs_rad = βs .* π

            # Evaluate via real QAOA (GPU for large N_QUBITS)
            real_exp = qaoa_expectation_device(inst.costs, N_QUBITS, γs_rad, βs_rad)
            ratio = real_exp / optimal_costs[i]
            push!(ratios, ratio)

            println("    Instance $i: ratio = $(round(ratio, digits=4))")
        end

        all_results[label][p] = ratios
    end
end

#==============================================================================#
#                          PLOTTING                                             #
#==============================================================================#

println("\nPlotting...")
n_proxies = length(PROXY_CONFIGS)
proxy_labels = [cfg[1] for cfg in PROXY_CONFIGS]
colors = [:steelblue, :coral, :mediumseagreen, :orange, :purple]

fig = Figure(size=(800, 500))
ax = Axis(fig[1, 1],
    xlabel="QAOA Depth p",
    ylabel="Approximation Ratio",
    title="Proxy-Optimized Linear Ramp Performance\nG($N_QUBITS, $P_EDGE), $NUM_INSTANCES instances",
    xticks=P_VALUES,
)

# Width and offset for grouped box plots
total_width = 0.6
single_width = total_width / n_proxies

for (proxy_idx, label) in enumerate(proxy_labels)
    # Offset each proxy's boxes so they don't overlap
    offset = (proxy_idx - (n_proxies + 1) / 2) * single_width

    positions = Float64[]
    values = Float64[]

    for p in P_VALUES
        ratios = all_results[label][p]
        for r in ratios
            push!(positions, p + offset)
            push!(values, r)
        end
    end

    boxplot!(ax, positions, values,
        color=colors[mod1(proxy_idx, length(colors))],
        width=single_width * 0.8,
        label=label,
    )
end

axislegend(ax, position=:lt)

save_figure(fig, "figure7_high_depth_performance.png")
println("Done!")
