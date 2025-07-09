

#%% imports 
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import qokit.maxcut as mc
from grips.QAOA_simulator import get_expectation, get_simulator
from grips.QAOA_proxy_interface import QAOA_proxy, QAOA_proxy_expectation
import grips.triangle_proxy as tpr
import grips.paper_proxy as ppr
import grips.normal_proxy as npr
import os  

from juliacall import Main as jl
jl.seval('''
using Pkg
Pkg.activate(joinpath(@__DIR__, "../julia"))
Pkg.instantiate()
using JuliaQAOA
''')
#%%
edge_probability = 0.5
num_constraints = 10 # Number of constraints
num_qubits = 10  # Number of qubits
gammas = np.linspace(0, np.pi, 10)  # Gamma values for QAOA
betas = np.linspace(0, np.pi, 10)  # Beta values for QAOA


paper_proxy = ppr.PaperProxy(num_constraints, num_qubits, edge_probability)
triangle_proxy = tpr.TriangleProxy(num_constraints, num_qubits)
paper_results = QAOA_proxy(paper_proxy, gammas, betas)
triangle_results = QAOA_proxy(triangle_proxy, gammas, betas)

final_amplitudes = triangle_results[-1]
expectation = QAOA_proxy_expectation(triangle_proxy, final_amplitudes)
print("Expectation value:", expectation)



