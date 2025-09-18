using JuliaQAOA, BenchmarkTools, CUDA

# Get compilation out of the way with a small example
println("#"^40)
println("Benchmarking GPU vs CPU")
println("\nFor a graph with N qubits/vertices, we will assume N^2/2 edges/constraints/")
println("#"^40)

println("-"^40)
println("Homogeneous distribution computation for a single proxy")
println("-"^40)

function benchmark_gpu_compute_homodist(proxy)
    CUDA.@sync homodist = JuliaQAOA.gpu_compute_homodist(proxy) 
    return homodist
end


for num_qubits in 5:5:30
    num_constraints = div(num_qubits^2, 2)
    proxy = HardCodedTriangleProxy(num_constraints, num_qubits)
    homodist_size = (1+num_qubits)*(1+num_constraints)^2
    println("$num_qubits qubits, $num_constraints constraints, problem size $homodist_size")
    println("CPU")
    @btime JuliaQAOA.cpu_compute_homodist($proxy)
    println("GPU")
    @btime benchmark_gpu_compute_homodist($proxy)
    println()
end

println("Finished!")
