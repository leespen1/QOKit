using Test

@testset "QAOA Proxy" begin
    include("QAOA_proxy.jl")
end
@testset "Distribution Generation" begin
    include("test_gpu_distribution_generation.jl")
end

