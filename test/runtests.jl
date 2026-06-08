include(joinpath(@__DIR__, "..", "find_python.jl"))

using Test

@testset "QAOA Proxy" begin
    include("test_QAOA.jl")
end
@testset "Distribution Generation" begin
    include("test_gpu_distribution_generation.jl")
end
@testset "Real Distribution" begin
    include("test_cost_distribution.jl")
end
@testset "QAOA Simulation (analytical)" begin
    include("test_qaoa_analytical.jl")
end
@testset "QAOA Simulation (vs QOKit Python)" begin
    include("test_qaoa_simulation.jl")
end
@testset "GPU QAOA Simulation (vs CPU)" begin
    include("test_qaoa_simulation_gpu.jl")
end

