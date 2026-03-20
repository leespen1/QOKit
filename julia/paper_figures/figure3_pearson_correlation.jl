#=
figure3_pearson_correlation.jl — Reproduce Figure 3 from the paper.

Paper Figure 3: Pearson correlation coefficients between the empirical
(instance-averaged) N(c';d,c) and the analytical proxy N(c';d,c), plotted
for each c'. Includes P(c') overlay showing dominant terms, and insert
heatmaps comparing the distributions for selected c' values.

Paper parameters: G(10, 1/3), 10 instances
Quick test:       G(6, 0.5), 3 instances

Customization:
  - Swap the "analytical method" (proxy) by changing PROXY_CONFIGS
  - Change INSERT_COST_PRIMES to select which c' values get heatmap subplots
  - Compare multiple proxies on the same correlation plot

Started:  2026-03-17
Finished: 2026-03-17
=#

include("common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

# Graph parameters
const N_QUBITS = 10
const P_EDGE = 1/3
const NUM_INSTANCES = 10
const SEED = 42

# Which proxy types to compare? Each entry is (label, proxy_constructor)
# The constructor takes (num_edges, num_qubits) and returns a proxy.
PROXY_CONFIGS = [
    ("PaperProxy", (m, n) -> PaperProxy(m, n, P_EDGE)),
    ("TriangleProxy", (m, n) -> OldTriangleProxy(m, n)),
]

# Which c' values to show as insert heatmaps (paper uses c'=7, c'=13)
const INSERT_COST_PRIMES = [7, 13]  # Paper uses c'=7 and c'=13

#==============================================================================#
#                          COMPUTATION                                          #
#==============================================================================#

println("=== Figure 3: Pearson Correlation ===")
println("Parameters: n=$N_QUBITS, p_edge=$P_EDGE, instances=$NUM_INSTANCES")

# Step 1: Generate graphs and compute empirical averaged N(c';d,c)
rng = MersenneTwister(SEED)
instances = [generate_er_instance(N_QUBITS, P_EDGE; rng) for _ in 1:NUM_INSTANCES]
max_edges = maximum(inst.num_edges for inst in instances)

println("Computing empirical homogeneous distributions...")
homodists = map(instances) do inst
    get_homogeneous_distribution_from_costs_direct(
        inst.costs, inst.num_edges, inst.num_vertices;
        max_num_edges=max_edges
    )
end
empirical_homodist = average_distributions(homodists)
num_costs = size(empirical_homodist, 1)
m = num_costs - 1  # max edges for the padded distribution

# Step 2: For each proxy, compute analytical N(c';d,c) and Pearson correlations
println("Computing proxy distributions and correlations...")
proxy_results = map(PROXY_CONFIGS) do (label, constructor)
    proxy = constructor(m, N_QUBITS)
    proxy_homodist = cpu_compute_homodist(proxy)

    # Pad to match empirical
    proxy_homodist, emp = pad_to_match(proxy_homodist, empirical_homodist)

    # Pearson correlation for each c'
    correlations = get_pearson_correlation_coefficients(proxy_homodist, emp)

    # P(c') values
    P_values = [P_cost_distribution(proxy, c) for c in 0:m]

    (label=label, proxy=proxy, homodist=proxy_homodist, correlations=correlations, P=P_values)
end

#==============================================================================#
#                          PLOTTING                                             #
#==============================================================================#

println("Plotting...")

# Count total subplots: 1 main + 2 per insert c' (one per each c', showing
# empirical vs first proxy)
n_inserts = length(INSERT_COST_PRIMES)

fig = Figure(size=(1200, 400 + 300 * ceil(Int, n_inserts / 2)))

# --- Main plot: Pearson correlation coefficients vs c' ---
ax_main = Axis(fig[1, 1:2],
    xlabel="Cost c'",
    ylabel="Pearson Correlation / P(c')",
    title="Pearson Correlation: Empirical vs Analytical N(c';d,c)\nG($N_QUBITS, $P_EDGE), $NUM_INSTANCES instances",
)

# Plot correlation for each proxy
colors = [:blue, :red, :green, :orange, :purple]
for (i, res) in enumerate(proxy_results)
    # Replace NaN correlations with 0 for plotting
    corr_clean = replace(res.correlations, NaN => 0.0)
    scatterlines!(ax_main, 0:m, corr_clean,
        label="Corr: $(res.label)",
        color=colors[mod1(i, length(colors))],
        markersize=8,
    )
end

# Overlay P(c') (scaled to fit on same axes)
P_vals = proxy_results[1].P
P_max = maximum(P_vals)
P_scaled = P_vals ./ P_max  # Scale to [0, 1] range
barplot!(ax_main, collect(0:m), P_scaled,
    color=(:gray, 0.3),
    label="P(c') (scaled)",
)
axislegend(ax_main, position=:lb)

# --- Insert heatmaps: empirical vs analytical for selected c' values ---
for (idx, c_prime) in enumerate(INSERT_COST_PRIMES)
    if c_prime >= num_costs
        println("Skipping insert c'=$c_prime (out of range)")
        continue
    end

    row = 1 + ceil(Int, idx / 2)
    col = 2 - (idx % 2)  # alternating columns

    # Empirical slice
    emp_slice = empirical_homodist[c_prime + 1, :, :]
    # Analytical slice (first proxy)
    ana_slice = proxy_results[1].homodist[c_prime + 1, :, :]

    # Side-by-side heatmaps
    ax_emp = Axis(fig[row, 2*col - 1],
        xlabel="Cost c", ylabel="Distance d",
        title="Empirical, c'=$c_prime",
    )
    ax_ana = Axis(fig[row, 2*col],
        xlabel="Cost c", ylabel="Distance d",
        title="$(proxy_results[1].label), c'=$c_prime",
    )

    # Use same color range for comparison
    vmax = max(maximum(emp_slice), maximum(ana_slice))
    hm1 = heatmap!(ax_emp, 0:m, 0:N_QUBITS, emp_slice, colorrange=(0, vmax), colormap=:viridis)
    hm2 = heatmap!(ax_ana, 0:m, 0:N_QUBITS, ana_slice, colorrange=(0, vmax), colormap=:viridis)
    Colorbar(fig[row, 2*col + 1], hm2, label="N(c';d,c)")
end

save_figure(fig, "figure3_pearson_correlation.png")
println("Done!")
