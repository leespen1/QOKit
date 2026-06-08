import sys, os, numpy as np
sys.path.append(os.path.abspath("..")) # Allows us to import from grips and qokit directories
from grips import HardCodedTriangleProxy
from grips.QAOA_proxy_interface import QAOA_proxy, QAOA_proxy_expectation, QAOA_proxy_optimize_gamma_beta, jl
# The following imports and seval statements make Julia proxy functions available

num_constraints = 10
num_qubits = 6
triangle_proxy_julia = jl.HardCodedTriangleProxy(num_constraints, num_qubits)
triangle_proxy_python = HardCodedTriangleProxy(num_constraints, num_qubits)

np.random.seed(42)
gammas = np.random.rand(5)
betas = np.random.rand(5)

julia_proxy_amplitudes = QAOA_proxy(triangle_proxy_julia, gammas, betas) 
python_proxy_amplitudes = QAOA_proxy(triangle_proxy_python, gammas, betas) 
print(f"Julia proxy amplitudes: {julia_proxy_amplitudes}")
print(f"Python proxy amplitudes: {python_proxy_amplitudes}")
print("\n")


test_amplitudes = julia_proxy_amplitudes[-1]
julia_proxy_expectation = QAOA_proxy_expectation(triangle_proxy_julia, test_amplitudes)
python_proxy_expectation = QAOA_proxy_expectation(triangle_proxy_python, test_amplitudes)
print(f"Julia proxy expectation: {julia_proxy_expectation}")
print(f"Python proxy expectation: {python_proxy_expectation}")
print("\n")


return_dict_julia = QAOA_proxy_optimize_gamma_beta(
    triangle_proxy_julia, gammas, betas, optimizer_method="Nelder-Mead"
)
return_dict_python = QAOA_proxy_optimize_gamma_beta(
    triangle_proxy_python, gammas, betas, optimizer_method="Nelder-Mead"
)
print(f"Julia return dict: {return_dict_julia}")
print(f"Python return dict: {return_dict_python}")
print("\n")
