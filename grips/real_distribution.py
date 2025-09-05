################################################################################
#
# This file contains function definitions for taking a graph and computing the
# real distribution of costs n(x; d, c), and the cost-averaged homogeneous
# distribution N(c'; d, c).
#
################################################################################
import networkx as nx, qokit.maxcut as mc, numpy as np, typing
import numba
from numba import njit
from grips.QAOA_simulator import QAOA_run, get_simulator
import matplotlib.pyplot as plt

def get_homogeneous_distribution(graph):
    """
    Given a graph, compute the real distribution n(x; d, c) (as a 3D array),
    then compute N(c'; d, c) by averaging n(x; d, c) over each x with cost c'.
    """
    costs = get_costs(graph)
    num_edges = graph.number_of_edges()
    num_vertices = graph.number_of_nodes()
    real_distribution = get_real_distribution_from_costs(costs, num_edges, num_vertices)
    return get_homogeneous_distribution_from_costs(costs, real_distribution)



def get_real_distribution(graph):
    """
    Given a graph, compute the real distribution n(x; d, c) from the paper.

    n(x; d, c) gives the number of bitstrings that have cost c and are Hamming
    distance c from bitstring x.

    The distribution is implemented as an 3D array:
    - Index 1 chooses the bitstring, e.g. 000, 001, 010, 011, etc
    - Index 2 chooses the Hamming distance, e.g. 0, 1, 2, ..., num_qubits
    - Index 3 chooses the cost, e.g. 0, 1, 2, ... 

    This file is for taking a particular graph, and obtaining the real distribution
    n(x; d, c) from the paper. That is, the number of bitstrings with cost c that
    are Hamming distance d from the bitstring x.
    """

    costs = get_costs(graph)
    num_edges = graph.number_of_edges()
    num_vertices = graph.number_of_nodes()
    return get_real_distribution_from_costs(costs, num_edges, num_vertices)


def get_costs(graph):
    """
    Given a graph, get the cost associated with each bitstring (for the maxcut
    problem), as a 1D array.
    """
    num_vertices = graph.number_of_nodes()

    # Get the cost associated with each bitstring, as a 1D array
    ising_model = mc.get_maxcut_terms(graph)
    sim = get_simulator(num_vertices, ising_model)
    costs = sim.get_cost_diagonal()

    return costs



@njit
def get_real_distribution_from_costs(costs, num_edges, num_vertices):
    """
    Given a the set of costs associated with each bitstring for a graph
    partition, the number of edges in the graph, and the number of vertices in
    the graph, compute the real distribution n(x; d, c) from the paper.

    n(x; d, c) gives the number of bitstrings that have cost c and are Hamming
    distance c from bitstring x.

    The distribution is implemented as an 3D array:
    - Index 1 chooses the bitstring, e.g. 000, 001, 010, 011, etc
    - Index 2 chooses the Hamming distance, e.g. 0, 1, 2, ..., num_qubits
    - Index 3 chooses the cost, e.g. 0, 1, 2, ... 

    This file is for taking a particular graph, and obtaining the real distribution
    n(x; d, c) from the paper. That is, the number of bitstrings with cost c that
    are Hamming distance d from the bitstring x.
    """

    num_bitstrings = 2 ** num_vertices
    num_distances = 1 + num_vertices
    num_costs = 1 + num_edges

    assert len(costs.shape) == 1
    assert len(costs) == num_bitstrings

    n_distribution = np.zeros((num_bitstrings, num_distances, num_costs))

    for x in range(num_bitstrings): # 0:num_bitstrings-1
        for y in range(num_bitstrings):
            d = hamming_distance(x, y)
            cost_y = int(costs[y]) # Convert from float to int
            n_distribution[x, d, cost_y] += 1

    return n_distribution



@njit
def get_homogeneous_distribution_from_costs(costs, real_distribution):
    """
    Given the costs associated with each bitstring, and the real distribution
    n(x; d, c) (as a 3D array), compute N(c'; d, c) by averaging n(x; d, c)
    over each x with cost c'.
    """
    assert len(real_distribution.shape) == 3
    num_bitstrings = real_distribution.shape[0]
    num_distances = real_distribution.shape[1]
    num_costs = real_distribution.shape[2]

    homogeneous_distribution = np.zeros((num_costs, num_distances, num_costs))
    num_cost_occurences = np.zeros(num_costs)

    # Sum over bitstrings with cost c'
    for bitstring_i in range(num_bitstrings):
        bitstring_cost = int(costs[bitstring_i])
        num_cost_occurences[bitstring_cost] += 1

        homogeneous_distribution[bitstring_cost, :, :] += real_distribution[bitstring_i, :, :]

    # Take average over number of bitstrings with cost c'
    for bitstring_cost_int in range(num_costs):
        if num_cost_occurences[bitstring_cost_int] != 0:
            homogeneous_distribution[bitstring_cost_int, :, :] /= num_cost_occurences[bitstring_cost_int]
        # If num_cost_occurences[bitstring_cost_int] != 0,
        # then homogeneous_distribution[bitstring_cost_int, :, :] is already 0

    return homogeneous_distribution



@njit
def hamming_distance(bitstring1: int, bitstring2: int):
    """
    Compute the Hamming distance between two bitstrings, given as integers.

    The Hamming distance between bitstrings x and y is the number of bits which
    are different. E.g:
    - d(000, 111) = 3
    - d(100, 111) = 2
    - d(001, 111) = 2
    - d(101, 111) = 1
    - d(111, 111) = 0
    - d(010, 101) = 3
    """
    xor = bitstring1 ^ bitstring2
    d = bitcount(xor)
    # Bits of xor are 1 if-and-only-if the bits of x and y differ.
    # Therefore, xor.bit_count() is the number of bits that differ. 
    return d

@njit
def bitcount(x: int):
    """
    Counts the number of 1s in the binary representation of an integer x.
    """
    count = 0
    while x:
        count += x & 1
        x >>= 1
    return count

def pad_to_shape(arr, target_shape):
    """Pad arr with zeros to match target_shape."""
    pad_width = [(0, max(0, t - s)) for s, t in zip(arr.shape, target_shape)]
    return np.pad(arr, pad_width, mode='constant', constant_values=0)

def pad_and_stack(arrays):
    """Pad all arrays to the largest shape among them, then stack along a new axis (0)."""
    assert all(arr.ndim == 3 for arr in arrays), "All arrays must be 3D"
    max_shape = tuple(max(arr.shape[i] for arr in arrays) for i in range(3))
    padded = [pad_to_shape(arr, max_shape) for arr in arrays]
    return np.stack(padded, axis=0)

def average_distributions(distributions):
    """
    Compute the average of multiple 3D distributions.
    distributions: list of np.ndarray, each shape (c', d, c)
    Returns: np.ndarray of shape (max_c', max_d, max_c)
    """
    stacked = pad_and_stack(distributions)
    return np.mean(stacked, axis=0)

def stddev_distributions(distributions):
    """
    Compute the standard deviation of multiple 3D distributions.
    distributions: list of np.ndarray, each shape (c', d, c)
    Returns: np.ndarray of shape (max_c', max_d, max_c)
    """
    stacked = pad_and_stack(distributions)
    return np.std(stacked, axis=0)

def mean_and_stddev(distributions):
    """
    Compute mean and stddev of multiple 3D distributions.
    Returns: (mean, stddev), both np.ndarray of shape (max_c', max_d, max_c)
    """
    stacked = pad_and_stack(distributions)
    return np.mean(stacked, axis=0), np.std(stacked, axis=0)

def plot_stddev_div_mean_heatmap(distributions, cost):
    """
    Plot a heatmap of the mean divided by the standard deviation of multiple
    3D distributions.
    distributions: list of np.ndarray, each shape (c', d, c)
    cost: int, the cost to fix c' at for the heatmap (since we can only plot 2D)
    """
    mean, stddev = mean_and_stddev(distributions)
    ratio = stddev / mean
    ## Copilot suggestion, don't think it's necessary
    #with np.errstate(divide='ignore', invalid='ignore'):
    #    ratio = np.where(stddev != 0, mean / stddev, 0)
    plt.imshow(ratio, aspect='auto', origin='lower')
    plt.colorbar(label='Dev/Avg')
    plt.xlabel("Cost c")
    plt.ylabel("Hamming Distance d")
    plt.title("Mean divided by Stddev Heatmap")
    plt.title(f'Dev/Avg N from cost {cost}')
    plt.show()