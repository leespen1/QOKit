#=
figure5_objective_landscapes.jl — Reproduce Figure 5 from the paper.

Paper Figure 5: Side-by-side heatmaps of the "true" parameter objective
function (from full statevector simulation) and the "homogeneous proxy"
objective function. At p=3, the first two layers' parameters (γ₁,γ₂,β₁,β₂)
are fixed, and γ₃, β₃ are swept on a grid.

Paper parameters: G(8, 1/2), p=3, 30×30 grid
Quick test:       G(6, 0.5), p=2, 20×20 grid (fix γ₁,β₁, sweep γ₂,β₂)

Customization:
  - Change PROXY_CONFIGS to compare multiple proxies
  - Adjust GRID_SIZE, fixed parameters, sweep ranges

Started:  2026-03-17
Finished: 2026-03-17
=#

include("common.jl")

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_QUBITS = 6
const P_EDGE = 0.5
const P_DEPTH = 2            # Total QAOA depth (paper: 3)
const GRID_SIZE = 20         # Grid resolution (paper: 30)
const SEED = 42

# Fixed parameters for layers 1 to P_DEPTH-1 (in radians, QOKit convention)
# For p=2: fix γ₁, β₁; sweep γ₂, β₂
# For p=3: fix γ₁, γ₂, β₁, β₂; sweep γ₃, β₃
const FIXED_GAMMAS = [0.4]      # Length = P_DEPTH - 1
const FIXED_BETAS = [0.4]       # Length = P_DEPTH - 1

# Sweep range for the last layer's parameters (radians, QOKit convention)
const GAMMA_RANGE = range(0.0, 2.0, length=GRID_SIZE)
const BETA_RANGE = range(0.0, π/2, length=GRID_SIZE)

# Proxies to compare
PROXY_CONFIGS = [
    ("PaperProxy", (m, n) -> PaperProxy(m, n, P_EDGE)),
]

#==============================================================================#
#                          COMPUTATION                                          #
#==============================================================================#

println("=== Figure 5: Objective Function Landscapes ===")
println("Parameters: n=$N_QUBITS, p_edge=$P_EDGE, p=$P_DEPTH, grid=$(GRID_SIZE)×$(GRID_SIZE)")

# Generate one graph instance
rng = MersenneTwister(SEED)
inst = generate_er_instance(N_QUBITS, P_EDGE; rng)
m = inst.num_edges
println("Graph: $(inst.num_vertices) vertices, $m edges")

# --- True QAOA landscape ---
println("Computing true QAOA landscape...")
true_landscape = zeros(GRID_SIZE, GRID_SIZE)
for (i, γ_last) in enumerate(GAMMA_RANGE)
    for (j, β_last) in enumerate(BETA_RANGE)
        γs = vcat(FIXED_GAMMAS, γ_last)
        βs = vcat(FIXED_BETAS, β_last)
        true_landscape[i, j] = qaoa_expectation(inst.costs, N_QUBITS, γs, βs)
    end
end
println("  True landscape range: [$(minimum(true_landscape)), $(maximum(true_landscape))]")

# --- Proxy landscapes ---
println("Computing proxy landscapes...")
proxy_landscapes = map(PROXY_CONFIGS) do (label, constructor)
    proxy = constructor(m, N_QUBITS)
    homodist = cpu_compute_homodist(proxy)
    P_vals = [P_cost_distribution(proxy, c) for c in 0:m]

    landscape = zeros(GRID_SIZE, GRID_SIZE)

    # Build all (γ, β) parameter sets for batch evaluation
    # QAOA_proxy_multi expects K×p matrices
    K = GRID_SIZE * GRID_SIZE
    γ_matrix = zeros(K, P_DEPTH)
    β_matrix = zeros(K, P_DEPTH)

    idx = 0
    for (i, γ_last) in enumerate(GAMMA_RANGE)
        for (j, β_last) in enumerate(BETA_RANGE)
            idx += 1
            # Fill fixed layers (convert from radians to pi_units)
            for ℓ in 1:(P_DEPTH-1)
                γ_matrix[idx, ℓ] = FIXED_GAMMAS[ℓ] / π
                β_matrix[idx, ℓ] = FIXED_BETAS[ℓ] / π
            end
            # Swept layer
            γ_matrix[idx, P_DEPTH] = γ_last / π
            β_matrix[idx, P_DEPTH] = β_last / π
        end
    end

    # Run batch proxy evaluation
    Qs = QAOA_proxy_multi(homodist, γ_matrix, β_matrix; pi_units=true)
    Q_final = Qs[end]  # (m+1) × K matrix

    # Compute expectations
    expectations = expectation(Q_final, P_vals, N_QUBITS)

    # Reshape into grid
    landscape = reshape(vec(expectations), GRID_SIZE, GRID_SIZE)

    println("  $label landscape range: [$(minimum(landscape)), $(maximum(landscape))]")
    (label=label, landscape=landscape)
end

#==============================================================================#
#                          PLOTTING                                             #
#==============================================================================#

println("Plotting...")
n_panels = 1 + length(proxy_landscapes)
fig = Figure(size=(400 * n_panels + 100, 400))

# True QAOA landscape
ax1 = Axis(fig[1, 1],
    xlabel="γ_$P_DEPTH",
    ylabel="β_$P_DEPTH",
    title="True QAOA ⟨C⟩",
)
hm1 = heatmap!(ax1, collect(GAMMA_RANGE), collect(BETA_RANGE), true_landscape,
    colormap=:viridis)

# Proxy landscapes
for (i, res) in enumerate(proxy_landscapes)
    ax = Axis(fig[1, 1 + i],
        xlabel="γ_$P_DEPTH",
        ylabel="β_$P_DEPTH",
        title="$(res.label) ⟨C⟩",
    )
    hm = heatmap!(ax, collect(GAMMA_RANGE), collect(BETA_RANGE), res.landscape,
        colormap=:viridis)
    if i == length(proxy_landscapes)
        Colorbar(fig[1, 2 + i], hm, label="Expected Cost ⟨C⟩")
    end
end

# Add supertitle
Label(fig[0, :], "Objective Function Landscape — G($N_QUBITS, $P_EDGE), p=$P_DEPTH",
    fontsize=16, font=:bold)

save_figure(fig, "figure5_objective_landscapes.png")
println("Done!")
