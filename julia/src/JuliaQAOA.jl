module JuliaQAOA

using Distributions: pdf, Normal, MvNormal, Binomial, Multinomial
using StaticArrays: @SVector, @SMatrix, SVector
using PythonCall: PyArray, pyconvert
using LinearAlgebra: mul!
using Base.Threads: @threads, nthreads, threadid
using CUDA
using ProgressMeter: @showprogress
using InteractiveUtils
using Statistics: mean, std, cor
#using LoopVectorization: @turbo, vmap, vmapreduce, indices

# Functions that will be made available when I call "using JuliaQAOA"
export P_cost_distribution, N_cost_distribution, N_cost_distance_distribution
export NormalProxy, PaperProxy, TriangleProxy, HardCodedTriangleProxy, IntuitiveTriangleProxy, OldTriangleProxy
#export compute_amplitude_sum, QAOA_proxy, QAOA_proxy_expectation
#export inverse_proxy_objective_function, QAOA_proxy_expectation


include("utils.jl") # Has abstract type definitions! Include this first!
include("QAOA_proxy.jl")
export _expand, get_β_factors, get_γ_factors, QAOA_proxy_basic, QAOA_proxy_single, QAOA_proxy_multi, expectation
include("triangle_proxy.jl")
include("normal_proxy.jl")
include("paper_proxy.jl")
include("cost_distributions.jl")
include("linear_ramp.jl")
export linear_ramp, linear_ramp_matrix
export cpu_compute_homodist, gpu_compute_homodist, allocate_homodist
export cpu_multi_proxy_mse, gpu_multi_proxy_mse, sum_squared_error
include("qaoa_simulation.jl")
export maxcut_costs, apply_phase_gate!, apply_x_mixer!, qaoa_statevector, qaoa_expectation
include("qaoa_simulation_gpu.jl")
export gpu_apply_phase_gate!, gpu_apply_x_mixer!, gpu_qaoa_statevector, gpu_qaoa_expectation, gpu_maxcut_costs

end # Module

#=
To make the exported functions available in python, assuming you are using a
script in the grips_examples directory, add these lines to the top of the
script:

from juliacall import Main as jl
jl.seval('''
using Pkg
Pkg.activate(joinpath(@__DIR__, "../julia"))
Pkg.instantiate()
using JuliaQAOA
''')

Then, for example, you can do:

paper_proxy = jl.PaperProxy(num_constraints, num_qubits, prob_edge)
probability_of_cost = jl.P_cost_distribution(paper_proxy, cost)

=#
