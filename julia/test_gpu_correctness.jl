using JuliaQAOA, Test, Random

@testset "CPU and GPU homogeneous distributions agree for a single proxy" begin
    num_qubits = 10
    num_constraints = div(num_qubits^2, 2)
    proxy_hard = HardCodedTriangleProxy(num_constraints, num_qubits)
    proxy_soft = TriangleProxy(num_constraints, num_qubits)
    proxies = [proxy_hard, proxy_soft]
    for proxy in proxies
        cpu_homodist = JuliaQAOA.cpu_compute_homodist(proxy)
        gpu_homodist = JuliaQAOA.gpu_compute_homodist(proxy) |> Array

        @test all(isapprox.(cpu_homodist, gpu_homodist, rtol=1e-15, atol=1e-12))
        println("Maximum disagreement: ", maximum(abs, cpu_homodist .- gpu_homodist))
    end
end

@testset "CPU and GPU mean squared-error values agree for multiple proxies." begin
    num_qubits = 10
    num_constraints = div(num_qubits^2, 2)
    sampled_homodist = JuliaQAOA.allocate_homodist(num_constraints, num_qubits)

    @testset "Hardcoded Triangle" begin
        proxy_hard = HardCodedTriangleProxy(num_constraints, num_qubits)
        proxies = [proxy_hard, proxy_hard]
        cpu_mses = JuliaQAOA.cpu_multi_proxy_mse(proxies, sampled_homodist)
        gpu_mses = JuliaQAOA.gpu_multi_proxy_mse(proxies, sampled_homodist) |> Array |> vec
        @test all(isapprox.(cpu_mses, gpu_mses, rtol=1e-15, atol=1e-12))
        println("Maximum disagreement: ", maximum(abs, cpu_mses .- gpu_mses))
    end


    @testset "Parameterized Triangle" begin
        params_vec = [[0.0, 0, 1, 1], [5, 3, 0.9, 1.2]]
        proxies = [TriangleProxy(num_constraints, num_qubits, params...)
                   for params in params_vec]
        cpu_mses = JuliaQAOA.cpu_multi_proxy_mse(proxies, sampled_homodist)
        gpu_mses = JuliaQAOA.gpu_multi_proxy_mse(proxies, sampled_homodist) |> Array |> vec
        @test all(isapprox.(cpu_mses, gpu_mses, rtol=1e-15, atol=1e-12))
        println("Maximum disagreement: ", maximum(abs, cpu_mses .- gpu_mses))
    end
    
end

