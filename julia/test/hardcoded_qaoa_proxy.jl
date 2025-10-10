using JuliaQAOA
using LinearAlgebra: norm
# Test that the QAOA proxy agrees with a hardcoded/manufactured solution
N = zeros(2,2,2)
N[1,:,:] .= [1 3
             2 4]
N[2,:,:] .= [5 7
             6 8]

gammas = [pi]
betas = [3pi/4]


output_vec = JuliaQAOA.QAOA_proxy(N, gammas, betas)
output_vec_matrix_version = JuliaQAOA.QAOA_proxy_matrix(N, gammas, betas)
@show output_vec
@show output_vec_matrix_version
@assert isapprox(output_vec, [1+1im, 1+1im], rtol=1e-15)

# Test that original and matrix versions agree for larger, random example

N2 = rand(50,8,50)
gammas2 = rand(4)
betas2 = rand(4)
output_vec2 = JuliaQAOA.QAOA_proxy(N2, gammas2, betas2)
output_vec_matrix_version2 = JuliaQAOA.QAOA_proxy_matrix(N2, gammas2, betas2)
@show output_vec2[1]
@show output_vec_matrix_version2[1]
@show norm(output_vec2 .- output_vec_matrix_version2)
@show maximum(abs, output_vec2 .- output_vec_matrix_version2)

@assert isapprox(output_vec2, output_vec_matrix_version2, rtol=1e-14)
