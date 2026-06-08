#=
figure2_stddev_heatmap.jl — Reproduce Figure 2 from the paper.

Paper Figure 2: Heat map of stddev/mean of N(c'; d, c) across multiple
graph instances, for a fixed c'. Shows that the homogeneous approximation
(replacing per-bitstring n(x;d,c) with cost-averaged N(c';d,c)) is good
where it matters: near the center of the (d,c) distribution.

Paper parameters: G(10, 1/3), 10 instances, c' = 7
Quick test:       G(6, 0.5), 3 instances

Started:  2026-03-17
Finished: 2026-03-17
=#

include("common.jl")

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

# Graph parameters — change these to explore different graph families
const N_QUBITS = 10          # Number of vertices
const P_EDGE = 1/3           # Edge probability
const NUM_INSTANCES = 10     # Number of graph instances
const SEED = 42              # Random seed for reproducibility

# Which c' values to plot? Set to :all to save all, or a vector like [7]
const COST_PRIME_VALUES = :all

#==============================================================================#
#                          COMPUTATION                                          #
#==============================================================================#

println("=== Figure 2: Stddev/Mean Heatmap ===")
println("Parameters: n=$N_QUBITS, p_edge=$P_EDGE, instances=$NUM_INSTANCES")

# Step 1: Generate random graph instances and compute their homogeneous
#          distributions N(c'; d, c) individually.
rng = MersenneTwister(SEED)
instances = [generate_er_instance(N_QUBITS, P_EDGE; rng) for _ in 1:NUM_INSTANCES]

# Find the maximum number of edges across instances (for padding)
max_edges = maximum(inst.num_edges for inst in instances)
println("Edge counts: ", [inst.num_edges for inst in instances])
println("Max edges: $max_edges")

# Step 2: Compute homogeneous distribution for each instance.
#          Each is a 3D array N(c'; d, c) of shape (m+1, n+1, m+1).
println("Computing homogeneous distributions...")
homodists = map(instances) do inst
    get_homogeneous_distribution_from_costs_direct(
        inst.costs, inst.num_edges, inst.num_vertices;
        max_num_edges=max_edges
    )
end

# Step 3: Compute element-wise mean and stddev across instances.
mean_dist, stddev_dist = distributions_mean_and_stddev(homodists)
println("Distribution shape: $(size(mean_dist))")

#==============================================================================#
#                          PLOTTING                                             #
#==============================================================================#

# Determine which c' values to plot
num_costs = size(mean_dist, 1)
c_prime_list = if COST_PRIME_VALUES == :all
    collect(0:(num_costs - 1))
else
    COST_PRIME_VALUES
end

for c_prime in c_prime_list
    # Extract the 2D slice for this c': mean[c'+1, d, c] and stddev[c'+1, d, c]
    # Axes: x = cost c (0:m), y = Hamming distance d (0:n)
    mean_slice = mean_dist[c_prime + 1, :, :]    # (n+1) × (m+1)
    stddev_slice = stddev_dist[c_prime + 1, :, :] # (n+1) × (m+1)

    # Compute stddev/mean ratio; handle zeros (where N=0, ratio is undefined)
    ratio = similar(mean_slice)
    for idx in eachindex(mean_slice)
        if mean_slice[idx] > 0
            ratio[idx] = stddev_slice[idx] / mean_slice[idx]
        else
            ratio[idx] = NaN  # Will appear as gray in plot
        end
    end

    # Skip if all values are NaN (no bitstrings have this cost)
    if all(isnan, ratio)
        println("Skipping c'=$c_prime (no data)")
        continue
    end

    # Create figure
    fig = Figure(size=FIGURE_SIZE)
    ax = Axis(fig[1, 1],
        xlabel="Cost c",
        ylabel="Hamming Distance d",
        title="Stddev/Mean of N(c'=$c_prime; d, c)\nG($N_QUBITS, $P_EDGE), $NUM_INSTANCES instances",
    )

    # Plot heatmap: ratio is (n+1, m+1), rows=distance, cols=cost
    hm = heatmap!(ax, 0:(num_costs-1), 0:N_QUBITS, ratio,
        colormap=:viridis,
        nan_color=:gray80,  # Gray for undefined (N=0) cells
    )
    Colorbar(fig[1, 2], hm, label="Stddev / Mean")

    # Save
    save_figure(fig, "figure2_stddev_heatmap_cprime$(c_prime).png")
end

println("Done! Saved $(length(c_prime_list)) heatmap(s).")
