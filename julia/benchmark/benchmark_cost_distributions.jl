using JuliaQAOA, BenchmarkTools, Random, CUDA

println("Benchmarking cost distribution generation")
println("Julia threads: ", Base.Threads.nthreads())
has_gpu = CUDA.has_cuda_gpu()
if has_gpu
    println("GPU: ", CUDA.name(CUDA.device()))
else
    println("GPU: none detected")
end
println()

function make_random_costs(num_vertices, num_edges; seed=42)
    rng = MersenneTwister(seed)
    num_bitstrings = 2^num_vertices
    return Float64.(rand(rng, 0:num_edges, num_bitstrings))
end

# Use @belapsed for small n, single @elapsed for large n
function timed_run(f, n_threshold=14)
    return f()  # just call it, caller handles timing
end

# Warmup with small problem
warmup_costs = make_random_costs(4, 6)
get_real_distribution_from_costs(warmup_costs, 6, 4)
get_homogeneous_distribution_from_costs_direct(warmup_costs, 6, 4)
if has_gpu
    gpu_get_homogeneous_distribution_from_costs_direct(warmup_costs, 6, 4)
    CUDA.synchronize()
end

# --- CPU benchmarks ---
println("="^60)
println("CPU: get_real_distribution_from_costs")
println("="^60)
for n in 10:2:18
    num_edges = div(n * (n - 1), 4)  # ~p=0.5 Erdos-Renyi
    costs = make_random_costs(n, num_edges)
    if n <= 14
        t = @belapsed get_real_distribution_from_costs($costs, $num_edges, $n)
    else
        t = @elapsed get_real_distribution_from_costs(costs, num_edges, n)
    end
    pairs = 2.0^(2n)
    println("  n=$n  edges=$num_edges  time=$(round(t, sigdigits=3))s  pairs/sec=$(round(pairs/t, sigdigits=3))")
end

println()
println("="^60)
println("CPU: get_homogeneous_distribution_from_costs_direct")
println("="^60)
for n in 10:2:18
    num_edges = div(n * (n - 1), 4)
    costs = make_random_costs(n, num_edges)
    if n <= 14
        t = @belapsed get_homogeneous_distribution_from_costs_direct($costs, $num_edges, $n)
    else
        t = @elapsed get_homogeneous_distribution_from_costs_direct(costs, num_edges, n)
    end
    pairs = 2.0^(2n)
    println("  n=$n  edges=$num_edges  time=$(round(t, sigdigits=3))s  pairs/sec=$(round(pairs/t, sigdigits=3))")
end

# --- GPU benchmarks ---
if has_gpu
    println()
    println("="^60)
    println("GPU: gpu_get_homogeneous_distribution_from_costs_direct")
    println("="^60)
    for n in 10:2:24
        num_edges = div(n * (n - 1), 4)
        costs = make_random_costs(n, num_edges)
        # Warmup for this size
        gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, n)
        CUDA.synchronize()
        if n <= 18
            t = @belapsed begin
                gpu_get_homogeneous_distribution_from_costs_direct($costs, $num_edges, $n)
                CUDA.synchronize()
            end
        else
            t = CUDA.@elapsed begin
                gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, n)
            end
        end
        pairs = 2.0^(2n)
        println("  n=$n  edges=$num_edges  time=$(round(t, sigdigits=3))s  pairs/sec=$(round(pairs/t, sigdigits=3))")
    end

    # --- GPU vs CPU comparison ---
    println()
    println("="^60)
    println("GPU vs CPU speedup: get_homogeneous_distribution_from_costs_direct")
    println("="^60)
    for n in 10:2:18
        num_edges = div(n * (n - 1), 4)
        costs = make_random_costs(n, num_edges)
        # Warmup GPU
        gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, n)
        CUDA.synchronize()

        if n <= 14
            t_cpu = @belapsed get_homogeneous_distribution_from_costs_direct($costs, $num_edges, $n)
            t_gpu = @belapsed begin
                gpu_get_homogeneous_distribution_from_costs_direct($costs, $num_edges, $n)
                CUDA.synchronize()
            end
        else
            t_cpu = @elapsed get_homogeneous_distribution_from_costs_direct(costs, num_edges, n)
            GC.gc()
            CUDA.synchronize()
            t_gpu = CUDA.@elapsed gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, n)
        end
        speedup = t_cpu / t_gpu
        println("  n=$n  CPU=$(round(t_cpu, sigdigits=3))s  GPU=$(round(t_gpu, sigdigits=3))s  speedup=$(round(speedup, sigdigits=3))x")
    end
else
    println("\nSkipping GPU benchmarks (no GPU detected)")
end

println("\nDone.")
