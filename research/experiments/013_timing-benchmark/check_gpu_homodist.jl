# Correctness gate for the grid-strided GPU homodist kernel: the GPU result
# must match the CPU direct computation exactly before any timing is trusted.
using JuliaQAOA
using CUDA
using Random: MersenneTwister

CUDA.functional() || error("CUDA not functional; this check requires a GPU node")

for n in (10, 12)
    edges = erdos_renyi_edges(n, 0.5; rng=MersenneTwister(1))
    m = length(edges)
    costs = maxcut_costs(n, edges)
    cpu = get_homogeneous_distribution_from_costs_direct(costs, m, n)
    gpu = Array(gpu_get_homogeneous_distribution_from_costs_direct(costs, m, n))
    err = maximum(abs.(cpu .- gpu))
    err < 1e-12 || error("GPU homodist mismatch at n=$n: max abs error $err")
    println("GPU homodist matches CPU at n=$n (max abs error $err)")
end
println("check_gpu_homodist passed")
