module JuliaQAOA

using Distributions: pdf, Normal, MvNormal, Binomial, Multinomial
using StaticArrays: @SVector, @SMatrix, SVector
using PythonCall: PyArray, pyconvert
using LinearAlgebra: mul!
using Base.Threads: @threads

# Functions that will be made available when I call "using JuliaQAOA"
export P_cost_distribution, N_cost_distribution, N_cost_distance_distribution
export NormalProxy, PaperProxy, TriangleProxy, HardCodedTriangleProxy
export compute_amplitude_sum, QAOA_proxy, QAOA_proxy_expectation
export inverse_proxy_objective_function
export qaoa_proxy_circuit

include("QAOA_proxy_interface.jl")
include("QAOA_proxy_transfer_matrix.jl")
include("normal_proxy.jl")
include("paper_proxy.jl")
include("triangle_proxy.jl")
include("utils.jl")

end

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
