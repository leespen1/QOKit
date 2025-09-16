using JuliaQAOA, ProfileView, BenchmarkTools
println("Using this many threads: ", Base.Threads.nthreads())

small_proxy = HardCodedTriangleProxy(5, 4)

num_qubits = 20
num_constraints = div(num_qubits^2, 2)
proxy = HardCodedTriangleProxy(num_constraints, num_qubits)
@btime homodist = JuliaQAOA.efficient_order_homodist(proxy);
@btime simple_homodist = JuliaQAOA.simple_homodist(proxy);

profile = false
if profile
    #@profview JuliaQAOA.efficient_order_homodist(small_proxy)
    t1 = time()
    @profview JuliaQAOA.efficient_order_homodist(proxy)
    t2 = time()
    println("Profiling took $(t2 - t1) seconds.")
end
