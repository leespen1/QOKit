"""
Required arguments:
- num_qubits: int
- N_cost_mean: float
- N_cov_1: float
- N_cov_2: float
"""
import numpy as np
import math
import typing
import time
import scipy
from scipy.stats import multivariate_normal



###
### The following functions are defined, but QAOA_paper_proxy now uses the julia
### implementations instead of the python ones. So most of the functions in this
### file will never be called. But they can be called using QAOA_paper_proxy_python
###

# P(c') from paper
def prob_cost_norm(cost: int, prob_cost_mean: float, prob_cost_cov: float) -> float:
    return multivariate_normal.pdf(cost, mean = prob_cost_mean, cov = prob_cost_cov)


# N(c') from paper
def make_N_cost_distribution_normal(num_qubits: int, N_cost_mean: float, N_cov_1: float, N_cov_2: float)
    def N_cost_distribution(cost: int, num_qubits: int) -> float:
        scale = 1 << num_qubits
        return prob_cost_norm(cost, num_qubits/2, num_qubits/4) * scale

    return N_cost_distribution


# N(c'; d, c) from paper
def make_N_cost_distance_distribution_normal(num_qubits: int, N_cost_mean: float, N_cov_1: float, N_cov_2: float):

    def N_cost_distance_distribution(cost_1: int, cost_2: int, distance: int) -> float:
        N_distance_mean = num_qubits / 2

        P = np.matrix([[cost_1 - N_cost_mean, N_distance_mean], [-N_distance_mean, cost_1 - N_cost_mean]])
        P_inv = scipy.linalg.inv(P)
        cov_mat = P@np.matrix([[N_cov_1, 0], [0, N_cov_2]])@P_inv
        cov_mat[0, 1] = cov_mat[1, 0] # cov_mat must be symmetric and is prone to floating point error
        return multivariate_normal([N_cost_mean, N_distance_mean], cov_mat).pdf([cost_2, distance])*(1 << num_qubits), cov_mat

    return N_cost_distance_distribution
