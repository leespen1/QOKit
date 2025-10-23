using JuliaQAOA, ProfileView, BenchmarkTools, Random
println("Using this many threads: ", Base.Threads.nthreads())


num_qubits = 20
num_constraints = div(num_qubits^2, 2)
proxy = HardCodedTriangleProxy(num_constraints, num_qubits)
@btime homodist = JuliaQAOA.efficient_order_homodist(proxy);
@btime simple_homodist = JuliaQAOA.simple_homodist(proxy);

sampled_homodist = JuliaQAOA.allocate_homodist(proxy)
sampled_homodist .= rand(MersenneTwister(0), size(sampled_homodist))
mses_cpu = JuliaQAOA.cpu_triangle_proxy_sweep_homodist([proxy, proxy], sampled_homodist)

profile = false
if profile
    #@profview JuliaQAOA.efficient_order_homodist(small_proxy)
    t1 = time()
    @profview JuliaQAOA.efficient_order_homodist(proxy)
    t2 = time()
    println("Profiling took $(t2 - t1) seconds.")
end
