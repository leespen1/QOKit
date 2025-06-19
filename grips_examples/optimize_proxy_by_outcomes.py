#%% imports 
import os
import matplotlib.pyplot as plt
import networkx as nx
import qokit.maxcut as mc
import numpy as np
from grips.QAOA_paper_proxy import QAOA_paper_proxy
from grips.QAOA_simulator import get_simulator, get_expectation, get_result
from grips.plot_utils import plot_heat_map
dir_path = os.path.dirname(os.path.realpath(__file__))
from juliacall import Main as jl
jl.seval('using Pkg')
jl.seval('Pkg.activate(joinpath(@__DIR__, "..", "julia"))')
jl.seval('Pkg.instantiate()')
jl.seval(f'include("{dir_path}/../julia/QAOA_proxy_peak_distance.jl")')


'''This is a replacement for QAOA proxy, with the peak and distance as arguments.'''
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

def collect_parameter_data_proxy_peak_distance(
    num_constraints: int,
    num_qubits: int,
    gammas: np.ndarray,
    betas: np.ndarray,
    h_peak: float,
    distance: int
) -> np.ndarray:
    expectations = np.zeros([len(gammas), len(betas)])
    for i in range(len(gammas)):
        for j in range(len(betas)):
            gamma = np.array([gammas[i]])
            beta = np.array([betas[j]])
            _, expectations[i][j] = QAOA_proxy_peak_distance(
                1, gamma, beta, num_constraints, num_qubits, h_peak, distance, 0
            )
    return expectations


#%% Here starts the real stuff 
# Define parameter ranges
gammas = np.linspace(0, np.pi, 40)
betas = np.linspace(0, np.pi, 40)

#these are multipliers for the peak and distance
#setting them to 1 uses the previous default values
M_peak = 0 
M_distance = 0

# Probabilities for the Erdos-Renyi graph generation
for p in [0.5]: # probability for the Erdos-Renyi graph generation
    for N in range(2, 10): # range of nodes, very small for testing
        G = nx.erdos_renyi_graph(N, p, seed=18) # generate graphs
        M = G.number_of_edges()

        # Define paths
        base_path = f"data_for_Expectation_Heatmaps/Erdős_Rényi/ER_p={p}"
        paper_proxy_path = os.path.join(base_path, f"Paper_Proxy_N={N}_M={M}.npz")
        peak_distance_proxy_path = os.path.join(base_path, f"Peak_Distance_Proxy_N={N}_M={M}.npz")

        # Create directories if they do not exist
        os.makedirs(base_path, exist_ok=True)

        if os.path.exists(paper_proxy_path):
            data = np.load(paper_proxy_path)
            expectation_proxies = data['expectations']
            gammas = data['gammas']
            betas = data['betas']
        else:
            expectation_proxies = collect_parameter_data_proxy(M, N, gammas, betas)
            np.savez(paper_proxy_path, expectations=expectation_proxies, gammas=gammas, betas=betas)

        # Set h_peak and distance values here
        num_qubits = N
        h_peak = float(2**(num_qubits - 4)*M_peak)
        distance = int(num_qubits // 2 * M_distance)

        if os.path.exists(peak_distance_proxy_path):
            data = np.load(peak_distance_proxy_path)
            expectation_proxies_peak_distance = data['expectations']
            gammas = data['gammas']
            betas = data['betas']
        else:
            expectation_proxies_peak_distance = collect_parameter_data_proxy_peak_distance(
                M, N, gammas, betas, h_peak, distance
            )
            np.savez(peak_distance_proxy_path, expectations=expectation_proxies_peak_distance, gammas=gammas, betas=betas)

        # Define image save paths
        img_base_path = f"Expectation_Heatmaps/Erdős_Rényi/ER_p={p}"
        paper_proxy_img_path = os.path.join(img_base_path, f"paper_proxy_N={N}_M={M}.png")
        peak_distance_proxy_img_path = os.path.join(img_base_path, f"peak_distance_proxy_N={N}_M={M}.png")

        # Create directories for images if they do not exist
        os.makedirs(img_base_path, exist_ok=True)

        # Generate and save heatmaps
        _ = plot_heat_map(gammas, betas, expectation_proxies, f"Expectation Proxies from Paper (N={N},M={M})", "Gamma", "Beta")
        plt.savefig(paper_proxy_img_path)
        _ = plot_heat_map(gammas, betas, expectation_proxies_peak_distance, f"Expectation Proxies (Peak Distance) (N={N},M={M})", "Gamma", "Beta")
        plt.savefig(peak_distance_proxy_img_path)

        plt.show()

# %%
