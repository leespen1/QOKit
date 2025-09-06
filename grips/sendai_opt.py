import numpy as np
import grips

# def mse_dist_loss_direct(proxy, realdist, num_constraints, num_qubits):
    
#     predicted = np.zeros_like(realdist)
    
#     # Loop over all cost_1, distance, cost_2
#     #note: this is parallelizable
#     for cost_1 in range(num_constraints+1):
#         for distance in range(num_qubits+1):
#             for cost_2 in range(num_constraints+1):
#                 predicted[cost_2, distance, cost_1] = proxy.N_cost_distance_distribution(cost_1, distance, cost_2)
    
#     # Normalize both to sum to 1 (or same scale)
#     predicted_sum = predicted.sum()
#     if predicted_sum == 0:
#         # Return a large MSE if the predicted distribution is all zeros
#         # to penalize these parameters in the optimization.
#         return 10000.0  # Or some other large, constant error value

#     predicted /= predicted_sum
#     realdist_norm = realdist / realdist.sum()
    
#     # Compute MSE
#     mse = np.mean((predicted - realdist_norm)**2)

#     return mse



def fit_proxy_to_real(proxy, realdist, init_params, bounds, optimizer='smart random search', max_iter=1000,\
                        fail_til_shrink = 4, fail_til_end = 50,\
                        SD_vs_paramrange = 0.2):
    """
    Fits a proxy distribution to a given distribution realdist.
    Currently uses a smart random search 
    (perturb randomly, re-use helpful perturbations, shrink if we fail consecutively)

    This *should* be easy to reuse across different proxies. 

    Args:
        proxy (QAOA_proxy): The proxy model to fit.
        realdist (np.ndarray): The target real distribution.
        init_params (np.ndarray): Initial parameters for the proxy.
        bounds (list): A list of [lower, upper] bounds for each parameter.
        num_constraints (int): The number of constraints.
        num_qubits (int): The number of qubits.
        optimizer: The optimization algorithm to use. Defaults to 'smart random search'.
        max_iter = 1000: The maximum number of iterations for the optimizer. Defaults to 1000.
        fail_til_shrink = 4: Number of failures to improve til we shrink SDs of perturbs
        fail_til_end = 50: Number of failures to improve until we terminate. 
        SD_vs_paramrange = 0.1: Ratio (SD of perturbation)/(width of param bounds) for each perturbation
        
    Returns:
        tuple: A tuple containing:
            - np.ndarray: The best parameters found.
            - float: The final mean squared error loss.
    """
    if optimizer == 'smart random search':
        current_params = np.array(init_params)
        
        # Set initial proxy parameters
        proxy.set_params(current_params)
        current_mse = grips.distribution_mean_squared_error(proxy, realdist)
        
        bounds = np.array(bounds)
        param_ranges = bounds[:, 1] - bounds[:, 0]
        sds = param_ranges*SD_vs_paramrange #initial SDs of perturbations
        
        consecutive_failures = 0
        consecutive_failures_aftershrink = 0
        
        # Generate initial perturbation
        perturbation = np.random.normal(0, sds)

        for i in range(max_iter):
            perturbed_params = current_params + perturbation
            
            # Clip parameters to be within bounds
            perturbed_params = np.clip(perturbed_params, bounds[:, 0], bounds[:, 1])
            
            proxy.set_params(perturbed_params)
            new_mse = grips.distribution_mean_squared_error(proxy, realdist)
            if new_mse < current_mse:
                current_params = perturbed_params
                current_mse = new_mse
                consecutive_failures = 0
                consecutive_failures_aftershrink = 0
                # Re-use the same perturbation if it helped
            else:
                consecutive_failures += 1
                consecutive_failures_aftershrink += 1
                # Generate a new random perturbation
                perturbation = np.random.normal(0, sds)

            #if we failed consecutively, shrink SDs of the perturbations
            if consecutive_failures >= fail_til_shrink:
                sds *= (2.0 / 3.0) 
                consecutive_failures = 0
                # Generate a new random perturbation with the new sds
                perturbation = np.random.normal(0, sds)

            #stop if we fail too many times even after shrinking
            if consecutive_failures_aftershrink >= fail_til_end: 
                break
        return current_params, current_mse
    else:
        raise NotImplementedError(f"Optimizer '{optimizer}' is not implemented.")

