
'''TO DO: 
Modify this to also specify the type of proxy to use! 
Add paper proxy!!'''


#%% imports
print("Importing python packages/functions ...")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import qokit.maxcut as mc
import os
import grips
import pandas as pd
from grips import (
    get_simulator, QAOA_run, QAOA_proxy, QAOA_proxy_expectation, QAOA_proxy_optimize_gamma_beta,
    get_homogeneous_distribution,get_homogeneous_distribution_from_proxy,
    TriangleProxy, NormalProxy, PaperProxy,
    inverse_objective_function, get_expectation,  
    maxcut, maxcut_approx_ratio, spsa_for_scipy,
    plot_distribution_lines_all, fit_proxy_to_real, 
    pad_and_stack
)
from grips import distribution_mean_squared_error as mse_dist_loss
from scipy.optimize import minimize
from scipy.optimize import dual_annealing
print("Finished importing python packages/functions!")

#%%  Julia imports
print("Importing Julia functions ...")
from juliacall import Main as jl
jl.seval('''
using Pkg
Pkg.activate(joinpath(@__DIR__, "../julia"))
Pkg.instantiate()
using JuliaQAOA
''')
print("Finished importing Julia functions!")

# %% Set up a single graph, get its statistical homogeneous distribution
# 3 nodes, edge probability 0.3, seed=4 results in a graph with 2 edges

def fit_multiple_graphs(num_graphs, num_nodes, edge_probability = 0.5, graph_type = 'erdos_renyi',\
                         ws_num_neighbors = 2):

    print("\nSetting up graphs, getting homogeneous distributions ...")
    # num_nodes = 10
    # num_graphs = 100
    # edge_probability = 0.5
    graphs = []
    homodists = []

    assert graph_type in ['erdos_renyi', 'barabasi_albert', 'watts_strogatz'],\
        "graph_type must be one of 'erdos_renyi', 'barabasi_albert', or 'watts_strogatz'"
    
    if graph_type == 'watts_strogatz':
        print('Using edge_probability as rewiring probability for Watts-Strogatz graph.')
        assert ws_num_neighbors >= 2, "For Watts-Strogatz graph, ws_num_neighbors must be at least 2."

    import random
    for i in range(num_graphs):
        if graph_type == 'barabasi_albert':
            # Barabasi-Albert graph requires a parameter for number of edges to attach from a new node to existing nodes
            m = max(1, int(edge_probability * (num_nodes - 1)))  # Ensure at least one edge
            graph = nx.barabasi_albert_graph(num_nodes, m)
        elif graph_type == 'watts_strogatz':
            # Watts-Strogatz graph requires parameters for number of nearest neighbors and rewiring probability
            k = max(ws_num_neighbors, int(edge_probability * (num_nodes - 1)))  # Ensure at least two neighbors
            p = edge_probability  # Rewiring probability
            graph = nx.watts_strogatz_graph(num_nodes, k, p)
        else:
            graph = nx.erdos_renyi_graph(num_nodes, edge_probability)
        #approx ratios later break if we have graphs with no edges
        #add a random edge if the graph doesn't have one
        if graph.number_of_edges() == 0:
            nodes = list(graph.nodes)
            if len(nodes) >= 2:
                u, v = random.sample(nodes, 2)
                graph.add_edge(u, v)
        homodist = get_homogeneous_distribution(graph)
        graphs.append(graph)
        homodists.append(homodist)



    #Q: is the scaling correct here with np.mean()?
    homodists = pad_and_stack(homodists)
    homodist = np.mean(homodists, axis=0)

    print("Finished setting up graph and getting homogeneous distribution!")


    #% initial proxy values needed
    initial_params = [0, 0, 1, 1]  
    max_num_edges = max(g.number_of_edges() for g in graphs)
    num_constraints = max_num_edges # Number of constraints -- +1 here caused error previously!
    num_qubits = num_nodes  # Number of qubits

    # proxy = TriangleProxy(num_constraints, num_qubits, *initial_params)
    # initial_triangle_homodist = get_homogeneous_distribution_from_proxy(proxy)
    # plot_distribution_lines_all(initial_triangle_homodist, f"Initial Triangle Proxy Distribution {initial_params}")

    # % Define mean-squared error loss function for fitting triangle proxy to statistical homogeneous distribution
    print("\nDefining triangle proxy loss function ... ")

    # def mse_dist_loss(params, homodist, num_constraints=0):
    #     h_tweak_sub, hc_tweak_add, l_tweak_mul, r_tweak_mul = params
    #     num_constraints = max(homodist.shape[0] - 1, num_constraints)
    #     num_qubits = homodist.shape[1] - 1

    #     proxy = TriangleProxy(
    #         num_constraints=num_constraints,
    #         num_qubits=num_qubits,
    #         h_tweak_sub=h_tweak_sub,
    #         hc_tweak_add=hc_tweak_add,
    #         l_tweak_mul=l_tweak_mul,
    #         r_tweak_mul=r_tweak_mul
    #     )
    #     return grips.distribution_mean_squared_error(proxy, homodist)

    print("Finished defining triangle proxy loss function!")

    #% Fit triangle proxy to homogeneous distribution
    print("\nFitting triangle proxy to homogeneous distribution ...")

    # Initial guess and bounds
    # use_small_bounds = False

    ##PK: Using good bounds from data 
    bounds = [
        (0, num_nodes**2/3),   # h_tweak_sub >= 0
        (-10, 10),   # hc_tweak_add can be small positive/negative
        (0.005, 2),   # l_tweak_mul > 0
        (0.05, 2),   # r_tweak_mul > 0
    ]

    # # result = minimize(
    # #     mse_dist_loss,
    # #     initial_params,
    # #     args=(homodist, num_constraints, num_qubits),
    # #     method='Nelder-Mead',  
    # #     bounds=bounds,
    # #     options={'maxiter': 10000, 'epsilon': 0.0001}
    # # )

    # bounds[0] = (0,max_num_edges) #need non-None bounds 

    # small_bounds = [
    #     (0, 3),   # h_tweak_sub >= 0
    #     (-3, 3),   # hc_tweak_add can be small positive/negative
    #     (0.3, 5),   # l_tweak_mul > 0
    #     (0.3, 5),   # r_tweak_mul > 0
    # ]


    # if use_small_bounds:
    #     bounds = small_bounds

    # result = dual_annealing(
    #     mse_dist_loss,
    #     bounds=bounds,
    #     args=(homodist,),
    #     maxiter=100,
    #     #maxiter=50000  #this is almost definitely too many its but works for now :) 
    # )

    # fitted_params = result.x

    #note: fit_proxy_to_real is modifying the proxy in-place, 
    #so need a separate initial proxy here. 
    opt_init_proxy = TriangleProxy(
        num_constraints=num_constraints,
        num_qubits=num_qubits,
        h_tweak_sub=initial_params[0],
        hc_tweak_add=initial_params[1],
        l_tweak_mul=initial_params[2],
        r_tweak_mul=initial_params[3]
    )

    fitted_params, fitmse1 = fit_proxy_to_real(opt_init_proxy, homodist,\
                                        initial_params, bounds,\
                                            max_iter = 10000, fail_til_shrink = 50, 
                                            fail_til_end = 100, 
                                            grid_size_start= 9)
    print("Fitted parameters:", fitted_params)
    print("Fitted mse was ", fitmse1, "\n\n")

    fitted_proxy = TriangleProxy(
        num_constraints=num_constraints,
        num_qubits=num_qubits,
        h_tweak_sub=fitted_params[0],
        hc_tweak_add=fitted_params[1],
        l_tweak_mul=fitted_params[2],
        r_tweak_mul=fitted_params[3]
    )
    fitted_proxy_homodist = get_homogeneous_distribution_from_proxy(fitted_proxy)
    #plot_distribution_lines_all(fitted_proxy_homodist, f"Fitted Proxy Homogeneous Distribution for Random Graph {fitted_params}")

    print("Finished fitting triangle proxy to homogeneous distribution!")

    # % Compare MSE for initial and fitted parameters
    print("\nComparing MSE loss function for initial and fitted parameters ...")
    initial_proxy = TriangleProxy(
        num_constraints=num_constraints,
        num_qubits=num_qubits,
        h_tweak_sub=initial_params[0],
        hc_tweak_add=initial_params[1],
        l_tweak_mul=initial_params[2],
        r_tweak_mul=initial_params[3]
    )
    # Compare MSE for initial and fitted parameters
    initial_mse = mse_dist_loss(initial_proxy, homodist)
    fitted_mse = mse_dist_loss(fitted_proxy, homodist)
    print(f"Initial MSE: {initial_mse}")
    print(f"Fitted MSE: {fitted_mse}")
    print("Finished comparing MSE loss function for initial and fitted parameters!")


    #% Run QAOA with fitted proxy
    print("\nRunning QAOA with fitted proxy (but not tuned gamma/beta)...")
    gammas = np.linspace(0, np.pi, 10)  # Gamma values for QAOA
    betas = np.linspace(0, np.pi, 10)  # Beta values for QAOA
    fitted_triangle_results = QAOA_proxy(fitted_proxy, gammas, betas)
    # print("Fitted Proxy Results:", fitted_triangle_results)

    final_amplitudes = fitted_triangle_results[-1]
    expectation = QAOA_proxy_expectation(fitted_proxy, final_amplitudes)
    print("Expectation value:", expectation)
    print("Finished running QAOA with fitted proxy!")

    # % Run QAOA with initial, unfitted proxy
    print("\nRunning QAOA with initial, unfitted proxy ...")
    # Compute results with initial parameters


    initial_triangle_results = QAOA_proxy(initial_proxy, gammas, betas)
    # print("Initial Proxy Results:", initial_triangle_results)

    # initial_final_amplitudes = initial_triangle_results[-1]
    # initial_expectation = QAOA_proxy_expectation(initial_proxy, initial_final_amplitudes)
    # print("Initial Expectation value:", initial_expectation)

    # # Compare expectations
    # print(f"Initial Expectation: {initial_expectation}")
    # print(f"Fitted Expectation: {expectation}")
    print("Finished running QAOA with initial, unfitted proxy!")


    # Comparing QAOA results for initial versus fitted proxy
    '''
    -defining initial gammas and betas
    -finding best gamma and beta according to initial versus fitted proxy
    -using these to run QAOA and comparing expectations
    '''

    gamma_0 = np.array([0.1])
    beta_0 = np.array([0.1])
    init_result = QAOA_proxy_optimize_gamma_beta(initial_proxy, gamma_0, beta_0, optimizer_method = 'Nelder-Mead', optimizer_options={'maxiter': 1})
    gamma_init = init_result["gamma"]
    beta_init = init_result["beta"]



    fitted_result = QAOA_proxy_optimize_gamma_beta(fitted_proxy, gamma_0, beta_0, optimizer_method = 'Nelder-Mead', 
                                                optimizer_options = {'maxiter': 10000, 'epsilon': 0.00001})
    gamma_fitted = fitted_result["gamma"]
    beta_fitted = fitted_result["beta"]

    graph = nx.erdos_renyi_graph(num_nodes, edge_probability)


    #% get QAOA expectations of inverse objective function for initial and fitted proxies
    initial_expectations = []
    fitted_expectations = []
    initial_approx_ratios = []
    fitted_approx_ratios = []

    for graph in graphs:
        ising_model = mc.get_maxcut_terms(graph)
        N = graph.number_of_nodes()
        sim = get_simulator(N, ising_model)
        p = 1
        mixer = "x"
        expectations = []
        overlaps = []

        '''
        we want the QAOA expectations of the fitted gammas and bets from init proxy versus fitted proxy, 
        not using them as starter values to optimize. What is currently commented out is doing the latter, 
        which we don't want, but may want to look at later too. 
        '''

        # initres = QAOA_run(
        #     ising_model=ising_model,
        #     N=N,
        #     p=p,
        #     init_gamma=gamma_init,
        #     init_beta=beta_init,
        #     optimizer_method="COBYLA",
        #     optimizer_options=None,
        #     mixer=mixer,
        #     expectations=expectations,
        #     overlaps=overlaps
        # )
        inv_obj = inverse_objective_function(ising_model, N, p, mixer, None, None)
        initres = inv_obj(np.hstack([gamma_init, beta_init]))
        # print("Initial Proxy Result (negative expectation):", initres)

        # For fitted parameters
        expectations = []
        overlaps = []
        # fitres = QAOA_run(
        #     ising_model=ising_model,
        #     N=N,
        #     p=p,
        #     init_gamma=gamma_fitted,
        #     init_beta=beta_fitted,
        #     optimizer_method='COBYLA',
        #     optimizer_options=None,
        #     mixer=mixer,
        #     expectations=expectations,
        #     overlaps=overlaps
        # )
        inv_obj_fit = inverse_objective_function(ising_model, N, p, mixer, None, None)
        fitres = inv_obj_fit(np.hstack([gamma_fitted, beta_fitted]))
        # print("Fitted Proxy Result (negative expectation):", fitres)

        # Store the expectations (note: inverse_objective_function returns -expectation)
        initial_expectations.append(-initres)
        fitted_expectations.append(-fitres)

        initial_approx_ratios.append(maxcut_approx_ratio(graph, -initres))
        fitted_approx_ratios.append(maxcut_approx_ratio(graph, -fitres))


    # Print the mean expectations with 2 decimals
    print("\n\nInitial Proxy Mean Expectation: {:.2f}".format(np.mean(initial_expectations)))
    print("Fitted Proxy Mean Expectation : {:.2f}".format(np.mean(fitted_expectations)))

    print("\n\nInitial Proxy Mean Approx Ratio: {:.2f}".format(np.mean(initial_approx_ratios)))
    print("Fitted Proxy Mean Approx Ratio : {:.2f}".format(np.mean(fitted_approx_ratios)))

    print(f"\n\nInitial Proxy MSE: {initial_mse:.2f}")
    print(f"Fitted  Proxy MSE: {fitted_mse:.2f}")
    # %

    print(f'With {num_nodes} nodes, {num_graphs} graphs, edge probability {edge_probability}, \n'
        f'Fitted parameters were {[round(x,2) for x in fitted_params]}')
    

    results_dict = {
        'num_nodes': [num_nodes],
        'num_graphs': [num_graphs],
        'graph_type': [graph_type],
        'edge_probability': [edge_probability],
        'fitted_params': [fitted_params],
        'mean_initial_expectation': [np.mean(initial_expectations)],
        'mean_fitted_expectation': [np.mean(fitted_expectations)],
        'mean_initial_approx_ratio': [np.mean(initial_approx_ratios)],
        'mean_fitted_approx_ratio': [np.mean(fitted_approx_ratios)],
        'initial_mse': [initial_mse],
        'fitted_mse': [fitted_mse]
    }
    if graph_type == 'watts_strogatz':
        results_dict['ws_num_neighbors'] = [ws_num_neighbors]

    return pd.DataFrame(results_dict)

if __name__ == "__main__":
    # Run the fitting process
    num_graphs = 50
    edge_probability = 0.5
    startnodes = 4
    endnodes = 10
    graph_type = 'barabasi_albert'  # Change to 'barabasi_albert' or 'watts_strogatz' as needed

    print(f"\n\n---------------------\nFitting for {startnodes} nodes...")
    df = fit_multiple_graphs(num_graphs=num_graphs, num_nodes=startnodes,\
                            edge_probability=edge_probability, graph_type=graph_type)
        
    for num_nodes in range(startnodes+1, endnodes+1):
        print(f"\n\n---------------------\nFitting for {num_nodes} nodes...")
        df_new = fit_multiple_graphs(num_graphs=num_graphs, num_nodes=num_nodes,\
                                    edge_probability=edge_probability, graph_type=graph_type)
        df = pd.concat([df, df_new], ignore_index=True)
        
    print(df)

    # Save the DataFrame to a CSV file
    output_file = f"fitting_results_{graph_type}_{edge_probability}edgeprob_{startnodes}to{endnodes}nodes_{num_graphs}graphs.csv"
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
