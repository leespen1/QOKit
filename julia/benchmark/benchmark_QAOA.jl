using JuliaQAOA, BenchmarkTools, CUDA
import Random: MersenneTwister

suite = BenchmarkGroup()

## Long version
#nm_pairs = [(n=5, m=10), (n=10, m=50), (n=20, m=210), (n=30, m=450)]
#p_vals = [1, 2, 4, 8]
#num_param_sets_vals = [2, 3]

# Short version
nm_pairs = [(n=5, m=10), (n=10, m=50)]
p_vals = [1, 2]
num_param_sets_vals = [2, 3]

for pair in nm_pairs
    n, m = pair
    N = rand(MersenneTwister(0), Float32, 1+m, 1+n, 1+m)
    for p in p_vals
        params_tup = (n=n, m=m, p=p, num_param_sets=1)
        this_suite = BenchmarkGroup()
        suite[string(params_tup)] = this_suite

        # Single set of params
        γs = rand(MersenneTwister(1), Float32, p)
        βs = rand(MersenneTwister(2), Float32, p)

        this_suite["basic"] = @benchmarkable QAOA_proxy_basic($N, $γs, $βs)
        this_suite["single, no blas"] = @benchmarkable QAOA_proxy_single($N, $γs, $βs, blas=false)
        this_suite["single, blas"] = @benchmarkable QAOA_proxy_single($N, $γs, $βs, blas=true)
        this_suite["multi, no blas"] = @benchmarkable QAOA_proxy_multi($N, _expand($γs), _expand($βs), blas=false)
        this_suite["multi, blas"] = @benchmarkable QAOA_proxy_multi($N, _expand($γs), _expand($βs), blas=true)

        # Multiple sets of params
        for num_param_sets in num_param_sets_vals
            params_tup = (n=n, m=m, p=p, num_param_sets=num_param_sets)
            this_suite = BenchmarkGroup()
            suite[string(params_tup)] = this_suite

            γs = rand(MersenneTwister(1), Float32, num_param_sets, p)
            βs = rand(MersenneTwister(2), Float32, num_param_sets, p)

            this_suite["basic"] = @benchmarkable [QAOA_proxy_basic($N, $γs[i,:], $βs[i,:]) for i in 1:$num_param_sets]
            this_suite["single, no blas"] = @benchmarkable [QAOA_proxy_single($N, $γs[i,:], $βs[i,:], blas=false) for i in 1:$num_param_sets]
            this_suite["single, blas"]    = @benchmarkable [QAOA_proxy_single($N, $γs[i,:], $βs[i,:], blas=true) for i in 1:$num_param_sets]
            this_suite["multi, no blas"] = @benchmarkable QAOA_proxy_multi($N, $γs, $βs, blas=false)
            this_suite["multi, blas"]    = @benchmarkable QAOA_proxy_multi($N, $γs, $βs, blas=true)
        end
    end
end


println("Tuning ...")
tune_ret = tune!(suite)
println("Finished tuning!")
println("Running ...")
run_ret = run(suite, verbose=true)
println("Finished running!")
