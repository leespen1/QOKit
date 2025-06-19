import numpy as np
import time
import scipy
import typing
import os
import warnings
from juliacall import Main as jl
dir_path = os.path.dirname(os.path.realpath(__file__))
jl.seval('using Pkg')
jl.seval('Pkg.activate(joinpath(@__DIR__, "..", "julia"))')
jl.seval('Pkg.instantiate()')
jl.seval(f'include("{dir_path}/../julia/QAOA_proxy_peak_distance.jl")')

"""
(written in Julia)
"""
def QAOA_proxy_peak_distance(
    p: int,
    gamma: np.ndarray,
    beta: np.ndarray,
    num_constraints: int,
    num_qubits: int,
    h_peak: float,
    distance: int,
    terms_to_drop_in_expectation: int = 0,
):
    return jl.QAOA_proxy_peak_distance(p, gamma, beta, num_constraints, num_qubits, h_peak, distance, terms_to_drop_in_expectation)
