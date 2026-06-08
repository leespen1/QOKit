using JuliaQAOA, Test
using CUDA

if CUDA.has_cuda_gpu()
    @testset "CPU and GPU homogeneous distributions agree for a single proxy" begin
        num_qubits = 10
        num_constraints = div(num_qubits^2, 2)
        proxy_hard = HardCodedTriangleProxy(num_constraints, num_qubits)
        proxy_soft = IntuitiveTriangleProxy(num_constraints, num_qubits)
        proxy_soft2 = IntuitiveTriangleProxy(num_constraints, num_qubits, 2, 0.4, 0.3, 0.2)
        proxies = [proxy_hard, proxy_soft, proxy_soft2]
        for proxy in proxies
            cpu_homodist = JuliaQAOA.cpu_compute_homodist(proxy)
            gpu_homodist = JuliaQAOA.gpu_compute_homodist(proxy) |> Array
            @test all(isapprox.(cpu_homodist, gpu_homodist, rtol=1e-15, atol=1e-12))
            println("Maximum disagreement: ", maximum(abs, cpu_homodist .- gpu_homodist))
        end
    end
end

@testset "CPU and GPU mean squared-error values agree for multiple proxies." begin
    num_qubits = 10
    num_constraints = div(num_qubits^2, 2)
    # Dummy homogeneous distribution
    sampled_homodist = JuliaQAOA.allocate_homodist(num_constraints, num_qubits)
    @show maximum(sampled_homodist)

    @testset "Hardcoded Triangle" begin
        proxy_hard = HardCodedTriangleProxy(num_constraints, num_qubits)
        proxies = [proxy_hard, proxy_hard]
        cpu_mses = JuliaQAOA.cpu_multi_proxy_mse(proxies, sampled_homodist)
        gpu_mses = JuliaQAOA.gpu_multi_proxy_mse(proxies, sampled_homodist) |> Array |> vec
        @test all(isapprox.(cpu_mses, gpu_mses, rtol=1e-15, atol=1e-12))
        display(cpu_mses)
        display(gpu_mses)
        println("Maximum disagreement: ", maximum(abs, cpu_mses .- gpu_mses))
    end


    @testset "Parameterized Triangle" begin
        params_vec = [[1.0, 0.5, 0.25, 0.25], [2, 0.4, 0.3, 0.2]]
        proxies = [IntuitiveTriangleProxy(num_constraints, num_qubits, params...)
                   for params in params_vec]
        cpu_mses = JuliaQAOA.cpu_multi_proxy_mse(proxies, sampled_homodist)
        gpu_mses = JuliaQAOA.gpu_multi_proxy_mse(proxies, sampled_homodist) |> Array |> vec
        display(cpu_mses)
        display(gpu_mses)
        @test all(isapprox.(cpu_mses, gpu_mses, rtol=1e-15, atol=1e-12))
        println("Maximum disagreement: ", maximum(abs, cpu_mses .- gpu_mses))
    end
    
end
