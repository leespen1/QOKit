using JuliaQAOA, Random
println("Using this many threads: ", Base.Threads.nthreads())
p = 5
gammas = rand(MersenneTwister(0), p)
betas = rand(MersenneTwister(1), p)

# Do small example to compile functions
small_proxy = HardCodedTriangleProxy(5, 4)
qaoa_proxy_circuit(small_proxy, gammas, betas)
QAOA_proxy(small_proxy, gammas, betas)

num_qubits = 60
num_constraints = div(num_qubits^2, 2)
proxy = HardCodedTriangleProxy(num_constraints, num_qubits)
#println("Running old version")
#@time old_result = QAOA_proxy(proxy, gammas, betas)
println("Running new version")
@time new_result = qaoa_proxy_circuit(proxy, gammas, betas)

print_matrices = false
if print_matrices
    println("Old")
    display(old_result)
    println("New")
    display(new_result)
end
#println("maximum(abs, new_result - old_result) = ", maximum(abs, new_result - old_result))

#@profview old_result = QAOA_proxy(proxy, gammas, betas)
@profview new_result = qaoa_proxy_circuit(proxy, gammas, betas)