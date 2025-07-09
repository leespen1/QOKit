#%%  imports 
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
import grips.real_distribution as rd 

from juliacall import Main as jl
jl.seval('''
using Pkg
Pkg.activate(joinpath(@__DIR__, "../julia"))
Pkg.instantiate()
using JuliaQAOA
''')

# %%
num_edges = 7
edge_probability = 0.5
graphs = [nx.erdos_renyi_graph(num_edges, edge_probability) for _ in range(5)]

dist = rd.get_homogeneous_distribution(graphs)
print(dist)
dist.shape

# %%

