using DrWatson
@quickactivate "JuliaQAOA"
using JuliaQAOA, BenchmarkTools, Random, CUDA

function benchmark_orig(N, gammas, betas)
    return [JuliaQAOA.QAOA_proxy(N, gammas[:,i], betas[:,i]) for i in 1:size(gammas, 2)]
end

function benchmark_matvec(N, gammas, betas; use_BLAS=true)
    return [JuliaQAOA.QAOA_proxy_matvec(N, gammas[:,i], betas[:,i], use_BLAS=use_BLAS) for i in 1:size(gammas, 2)]
end

function benchmark_matmat(N, gammas, betas; use_BLAS=true)
    return JuliaQAOA.QAOA_proxy_matmat(N, gammas, betas, use_BLAS=use_BLAS)
end

function benchmark_matmat_gpu(N, gammas, betas; use_BLAS=true)
    @assert (isa(N, CuArray) && isa(gammas, CuArray) && isa(betas, CuArray)) "Arrays must be CuArrays"
    CUDA.@sync begin
        result = JuliaQAOA.QAOA_proxy_matmat(N, gammas, betas, use_BLAS=use_BLAS)
    end
    return result
end

println('='^40, "\nBenchmarks\n", '='^40, '\n')

println("Small Sample\n", '-'^40, '\n')
N_small = rand(50,8,50)
gammas_small = rand(4,10)
betas_small = rand(4,10)

println("Original version")
@btime benchmark_orig($N_small, $gammas_small, $betas_small)
println("Matvec version with BLAS")
@btime benchmark_matvec($N_small, $gammas_small, $betas_small)
println("Matvec version without BLAS")
@btime benchmark_matvec($N_small, $gammas_small, $betas_small, use_BLAS=$false)
println("Matmat version with BLAS")
@btime benchmark_matmat($N_small, $gammas_small, $betas_small)
println("Matmat version without BLAS")
@btime benchmark_matmat($N_small, $gammas_small, $betas_small, use_BLAS=$false)
if CUDA.has_cuda_gpu()
    N_small_gpu = CuArray(N_small)
    gammas_small_gpu = CuArray(gammas_small)
    betas_small_gpu = CuArray(betas_small)
    println("GPU Matmat version with BLAS")
    @btime benchmark_matmat_gpu($N_small_gpu, $gammas_small_gpu, $betas_small_gpu)
    println("GPU Matmat version without BLAS")
    @btime benchmark_matmat_gpu($N_small_gpu, $gammas_small_gpu, $betas_small_gpu, use_BLAS=$false)
end


println('-'^40, "\nLarge Sample\n", '-'^40, '\n')
N_large = rand(190,20,190)
gammas_large = rand(4,100)
betas_large = rand(4,100)

println("Original version")
@btime benchmark_orig($N_large, $gammas_large, $betas_large)
println("Matvec version with BLAS")
@btime benchmark_matvec($N_large, $gammas_large, $betas_large)
println("Matvec version without BLAS")
@btime benchmark_matvec($N_large, $gammas_large, $betas_large, use_BLAS=$false)
println("Matmat version with BLAS")
@btime benchmark_matmat($N_large, $gammas_large, $betas_large)
println("Matmat version without BLAS")
@btime benchmark_matmat($N_large, $gammas_large, $betas_large, use_BLAS=$false)
if CUDA.has_cuda_gpu()
    N_large_gpu = CuArray(N_large)
    gammas_large_gpu = CuArray(gammas_large)
    betas_large_gpu = CuArray(betas_large)
    println("GPU Matmat version with BLAS")
    @btime benchmark_matmat_gpu($N_large_gpu, $gammas_large_gpu, $betas_large_gpu)
    println("GPU Matmat version without BLAS")
    @btime benchmark_matmat_gpu($N_large_gpu, $gammas_large_gpu, $betas_large_gpu, use_BLAS=$false)
end
