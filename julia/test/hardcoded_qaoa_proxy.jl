# Test that the QAOA proxy agrees with a hardcoded/manufactured solution
using JuliaQAOA
N = zeros(2,2,2)
N[1,:,:] .= [1 3
             2 4]
N[2,:,:] .= [5 7
             6 8]

gammas = [pi]
betas = [3pi/4]

output_vec = JuliaQAOA.QAOA_proxy(N, gammas, betas)

@assert isapprox(output_vec, [1+1im, 1+1im], rtol=1e-15)
