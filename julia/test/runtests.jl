using Test

@testset "QAOA Proxy" begin
    include("test_QAOA.jl")
end
@testset "Distribution Generation" begin
    include("test_gpu_distribution_generation.jl")
end
@testset "Real Distribution" begin
    include("test_real_distribution.jl")
end

