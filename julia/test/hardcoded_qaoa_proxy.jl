using JuliaQAOA, Test, BenchmarkTools
using LinearAlgebra: norm

benchmark = true
use_gpu = true

# Test that the QAOA proxy agrees with a hardcoded/manufactured solution
N = zeros(2,2,2)
N[1,:,:] .= [1 3
             2 4]
N[2,:,:] .= [5 7
             6 8]

gammas = [pi]
betas = [3pi/4]


output_vec = JuliaQAOA.QAOA_proxy(N, gammas, betas)
output_vec_matvec_version = JuliaQAOA.QAOA_proxy_matvec(N, gammas, betas)
output_vec_matmat_version = JuliaQAOA.QAOA_proxy_matmat(N, gammas, betas)
@show output_vec
@show output_vec_matvec_version
@show output_vec_matmat_version
@testset "Hardcoded Example" begin
    @test isapprox(output_vec, [1+1im, 1+1im], rtol=1e-15)
end

# Test that original and matrix versions agree for larger, random example

N2 = rand(50,8,50)
gammas2 = rand(4)
betas2 = rand(4)
output_vec2 = JuliaQAOA.QAOA_proxy(N2, gammas2, betas2)
output_vec_matvec_version2 = JuliaQAOA.QAOA_proxy_matvec(N2, gammas2, betas2)
output_vec_matmat_version2 = JuliaQAOA.QAOA_proxy_matmat(N2, gammas2, betas2)
@show output_vec2[1]
@show output_vec_matvec_version2[1]
@show output_vec_matmat_version2[1]
@show norm(output_vec2 .- output_vec_matvec_version2)
@show maximum(abs, output_vec2 .- output_vec_matvec_version2)
@show norm(output_vec_matmat_version2 .- output_vec_matvec_version2)
@show maximum(abs, output_vec_matmat_version2 .- output_vec_matvec_version2)

@testset "Random example, single set of parameters" begin
    @test isapprox(output_vec2, output_vec_matvec_version2, rtol=1e-14)
    @test isapprox(output_vec_matmat_version2, output_vec_matvec_version2, rtol=1e-14)
end

# Check agreement of orignal method and matmat  for multiple betas
N3 = rand(50,8,50)
gammas3 = rand(4,3)
betas3 = rand(4,3)

output_vecs3 = [JuliaQAOA.QAOA_proxy(N3, gammas3[:,i], betas3[:,i]) for i in 1:size(gammas3, 2)]
output_mat3 = hcat(output_vecs3...)

output_mat_matmat3 = JuliaQAOA.QAOA_proxy_matmat(N3, gammas3, betas3)
@testset "Random example, multiple sets of parameters" begin
    @test isapprox(output_mat3, output_mat_matmat3, rtol=1e-14)
end

function benchmark_orig(N, gammas, betas)
    return [JuliaQAOA.QAOA_proxy(N, gammas[:,i], betas[:,i]) for i in 1:size(gammas, 2)]
end

function benchmark_matvec(N, gammas, betas; use_BLAS=true)
    return [JuliaQAOA.QAOA_proxy_matvec(N, gammas[:,i], betas[:,i], use_BLAS=use_BLAS) for i in 1:size(gammas, 2)]
end

function benchmark_matmat(N, gammas, betas; use_BLAS=true)
    return JuliaQAOA.QAOA_proxy_matmat(N, gammas, betas, use_BLAS=use_BLAS)
end

if benchmark
    println('='^40, "\nBenchmarks\n", '='^40, '\n')

    println("Small Sample\n", '-'^40, '\n')
    N4 = rand(50,8,50)
    gammas4 = rand(4,10)
    betas4 = rand(4,10)

    println("Original version")
    @btime benchmark_orig($N4, $gammas4, $betas4)
    println("Matvec version with BLAS")
    @btime benchmark_matvec($N4, $gammas4, $betas4)
    println("Matvec version without BLAS")
    @btime benchmark_matvec($N4, $gammas4, $betas4, use_BLAS=$false)
    println("Matmat version with BLAS")
    @btime benchmark_matmat($N4, $gammas4, $betas4)
    println("Matmat version without BLAS")
    @btime benchmark_matmat($N4, $gammas4, $betas4, use_BLAS=$false)


    println('-'^40, "\nLarge Sample\n", '-'^40, '\n')
    N5 = rand(190,20,190)
    gammas5 = rand(4,100)
    betas5 = rand(4,100)

    println("Original version")
    @btime benchmark_orig($N5, $gammas5, $betas5)
    println("Matvec version with BLAS")
    @btime benchmark_matvec($N5, $gammas5, $betas5)
    println("Matvec version without BLAS")
    @btime benchmark_matvec($N5, $gammas5, $betas5, use_BLAS=$false)
    println("Matmat version with BLAS")
    @btime benchmark_matmat($N5, $gammas5, $betas5)
    println("Matmat version without BLAS")
    @btime benchmark_matmat($N5, $gammas5, $betas5, use_BLAS=$false)



    nothing
end
