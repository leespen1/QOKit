using DrWatson
@quickactivate "JuliaQAOA"
using JuliaQAOA, BenchmarkTools, CUDA, Random

do_single_proxy = false
do_multi_proxy = true

# Get compilation out of the way with a small example
println("#"^40)
println("Benchmarking GPU vs CPU")
println("\nFor a graph with N qubits/vertices, we will assume N^2/2 edges/constraints/")
println("#"^40)

if do_single_proxy 
    println("-"^40)
    println("Homogeneous distribution computation for a single proxy")
    println("-"^40)

    function benchmark_gpu_compute_homodist(proxy)
        CUDA.@sync homodist = JuliaQAOA.gpu_compute_homodist(proxy) 
        return homodist
    end


    for num_qubits in 5:5:30
        num_constraints = div(num_qubits^2, 2)
        homodist_size = (1+num_qubits)*(1+num_constraints)^2

        println("$num_qubits qubits, $num_constraints constraints, problem size $homodist_size")

        println("Hardcoded Triangle Proxy")
        hard_proxy = HardCodedTriangleProxy(num_constraints, num_qubits)
        println("CPU")
        @btime JuliaQAOA.cpu_compute_homodist($hard_proxy)
        println("GPU")
        @btime benchmark_gpu_compute_homodist($hard_proxy)

        println("\nParameterized Triangle Proxy")
        soft_proxy = TriangleProxy(num_constraints, num_qubits)
        println("CPU")
        @btime JuliaQAOA.cpu_compute_homodist($soft_proxy)
        println("GPU")
        @btime benchmark_gpu_compute_homodist($soft_proxy)
        println("-----\n")
    end
end

if do_multi_proxy
    println("-"^40)
    println("MSE computation for many proxies")
    println("-"^40)

    function benchmark_gpu_compute_mse(proxies)
        CUDA.@sync mses = JuliaQAOA.gpu_multi_proxy_mse(proxies) 
        return mses
    end


    for num_qubits in (5, 10)
        for num_proxies in (1, 10, 100, 1000, 10000)
            num_constraints = div(num_qubits^2, 2)
            homodist_size = (1+num_qubits)*(1+num_constraints)^2
            # Dummy homogeneous distribution
            sampled_homodist = JuliaQAOA.allocate_homodist(num_constraints, num_qubits)
            println("$num_proxies proxies, $num_qubits qubits, $num_constraints constraints, problem size $homodist_size")

            println("Hardcoded Triangle Proxy")
            hard_proxies = [HardCodedTriangleProxy(num_constraints, num_qubits)
                            for i in 1:num_proxies]
            println("CPU")
            @btime JuliaQAOA.cpu_multi_proxy_mse($hard_proxies, $sampled_homodist)
            println("GPU")
            @btime JuliaQAOA.gpu_multi_proxy_mse($hard_proxies, $sampled_homodist)

            println("\nParameterized Triangle Proxy")
            soft_proxies = [TriangleProxy(num_constraints, num_qubits, rand(MersenneTwister(i), 4)...)
                            for i in 1:num_proxies]
            println("CPU")
            @btime JuliaQAOA.cpu_multi_proxy_mse($soft_proxies, $sampled_homodist)
            println("GPU")
            @btime JuliaQAOA.gpu_multi_proxy_mse($soft_proxies, $sampled_homodist)
            println("-----\n")
        end
    end
end


println("Finished!")
