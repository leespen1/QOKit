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

def get_homogeneous_distribution(graphs, max_number_of_edges=0, simulator_name="auto"):
    """
    Given a graph (or graphs), compute the real distribution n(x; d, c) (as a 3D array),
    then compute N(c'; d, c) by averaging n(x; d, c) over each x with cost c'.

    max_number_of_edges can be provided to provide additional padding with
    zeros to make size compatible with homodists of other graphs which have the
    same number of vertices, but more edges.

    "The better that N(c'; d, c) estimates the average of n(x; d, c) over all
    c(x) = c', and the less n(x; d, c) deviates over all x with c(x) = c', the
    better N(c'; d, c) should estimate n(x; d, c) for all x with c(x) = c'" 
    """
    if isinstance(graphs, nx.Graph):
        graph = graphs 
        costs = get_costs(graph)
        num_edges = graph.number_of_edges()
        num_vertices = graph.number_of_nodes()
        real_distribution = get_real_distribution_from_costs(
            costs, num_edges, num_vertices, max_number_of_edges
        )
        homodist =  get_homogeneous_distribution_from_costs(
            costs, real_distribution, max_number_of_edges
        )
    elif isinstance(graphs, list) and all(isinstance(G, nx.Graph) for G in graphs):
        max_number_of_edges = max(max_number_of_edges,
                                  max(G.number_of_edges() for G in graphs))
        num_vertices = graphs[0].number_of_nodes()
        assert all(G.number_of_nodes() == num_vertices for G in graphs), "All graphs must have the same number of vertices."
        max_number_of_costs = 1+max_number_of_edges
        num_distances = 1+num_vertices
        homodist = np.zeros((max_number_of_costs, num_distances, max_number_of_costs))
        for graph in graphs:
            costs = get_costs(graph)
            num_edges = graph.number_of_edges()
            num_vertices = graph.number_of_nodes()
            real_distribution = get_real_distribution_from_costs(
                costs, num_edges, num_vertices, max_number_of_edges
            )
            homodist += get_homogeneous_distribution_from_costs(
                costs, real_distribution, max_number_of_edges
            )
        homodist /= len(graphs) # Take the average
    else:
        raise TypeError("graphs must be an nx.Graph or a list of nx.Graphs.")

    return homodist



def get_real_distribution(graph):
    """
    Given a graph, compute the real distribution n(x; d, c) from the paper.

    n(x; d, c) gives the number of bitstrings that have cost c and are Hamming
    distance c from bitstring x.

    The distribution is implemented as an 3D array:
    - Index 1 chooses the bitstring, e.g. 000, 001, 010, 011, etc
    - Index 2 chooses the Hamming distance, e.g. 0, 1, 2, ..., num_qubits
    - Index 3 chooses the cost, e.g. 0, 1, 2, ... 

    This function is for taking a particular graph, and obtaining the real
    distribution n(x; d, c) from the paper. That is, the number of bitstrings
    with cost c that are Hamming distance d from the bitstring x.
    """

    costs = get_costs(graph)
    num_edges = graph.number_of_edges()
    num_vertices = graph.number_of_nodes()
    return get_real_distribution_from_costs(costs, num_edges, num_vertices)


def get_costs(graph, simulator_name="auto"):
    """
    Given a graph, get the cost associated with each bitstring (for the maxcut
    problem), as a 1D array.
    """
    num_vertices = graph.number_of_nodes()

    # Get the cost associated with each bitstring, as a 1D array
    ising_model = mc.get_maxcut_terms(graph)
    sim = get_simulator(num_vertices, ising_model, simulator_name=simulator_name)
    costs = sim.get_cost_diagonal()

    return costs



@njit
def get_real_distribution_from_costs(costs, num_edges, num_vertices, max_num_edges=0):
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

    `max_num_edges` may be specified to allocate additional zeros (for
    compatibility across multiple graphs)
    """

    num_bitstrings = 2 ** num_vertices
    num_distances = 1 + num_vertices
    num_costs = 1 + num_edges
    cost_axis_size = max(num_costs, 1+max_num_edges)


    assert len(costs.shape) == 1, "costs must be a vector"
    assert len(costs) == num_bitstrings, "length of costs vector must match number of bitstrings"

    n_distribution = np.zeros((num_bitstrings, num_distances, cost_axis_size))

    for x in range(num_bitstrings): # 0:num_bitstrings-1
        for y in range(num_bitstrings):
            d = hamming_distance(x, y)
            cost_y = int(costs[y]) # Convert from float to int
            n_distribution[x, d, cost_y] += 1

    return n_distribution



@njit
def get_homogeneous_distribution_from_costs(costs, real_distribution, max_num_edges=0):
    """
    Given the costs associated with each bitstring, and the real distribution
    n(x; d, c) (as a 3D array), compute N(c'; d, c) by averaging n(x; d, c)
    over each x with cost c'.
    """
    assert len(real_distribution.shape) == 3
    num_bitstrings = real_distribution.shape[0]
    num_distances = real_distribution.shape[1]
    num_costs = real_distribution.shape[2]

    cost_axis_size = max(num_costs, 1+max_num_edges)

    homogeneous_distribution = np.zeros((cost_axis_size, num_distances, cost_axis_size))
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
    arr_shape = np.array(arr.shape)
    assert np.all(arr_shape <= target_shape), "arr shape must be <= target_shape elementwise"
    pad_width = [(0, max(0, t - s)) for s, t in zip(arr.shape, target_shape)]
    return np.pad(arr, pad_width, mode='constant', constant_values=0)


def pad_to_match(a, b):
    """
    Pads the smaller array to match the shape of the larger array.
    Assumes one shape is elementwise <= the other.
    """
    shape_a = np.array(a.shape)
    shape_b = np.array(b.shape)

    if np.array_equal(shape_a, shape_b):
        return a, b  # already the same shape

    # determine which is smaller
    if np.all(shape_a <= shape_b):
        smaller, larger = a, b
    elif np.all(shape_b <= shape_a):
        smaller, larger = b, a
    else:
        raise ValueError("Shapes are not broadcast-compatible for padding.")

    # compute padding widths
    pad_widths = [(0, t - s) for s, t in zip(smaller.shape, larger.shape)]
    padded = np.pad(smaller, pad_widths, constant_values=0)

    # return in original order
    if smaller is a:
        return padded, larger
    else:
        return larger, padded



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



def distributions_mean_and_stddev(distributions):
    """
    Compute mean and stddev of multiple 3D distributions.
    Returns: (mean, stddev), both np.ndarray of shape (max_c', max_d, max_c)
    """
    stacked = pad_and_stack(distributions)
    return np.mean(stacked, axis=0), np.std(stacked, axis=0)



def plot_distribution_lines(distribution, cost_prime):
    """
    Make a line plot where each line is N(c'; d, c) for fixed cost c' and
    Hamming distance d.
    """
    num_costs_prime, num_distances, num_costs = distribution.shape
    assert 0 <= cost_prime < num_costs_prime, "cost_prime out of bounds"
    plt.figure()
    for d in range(num_distances):
        plt.plot(range(num_costs), distribution[cost_prime, d, :], label=f'd={d}')
    plt.xlabel('Cost c')
    plt.ylabel(f'N(c\'={cost_prime}; d, c)')
    plt.title(f'Distribution Lines for cost c\'={cost_prime}')
    plt.legend()
    plt.show()



def plot_distribution_lines_all(distribution, suptitle="Distribution Lines for all c'"):
    """
    Create a grid of subplots, each showing N(c'; d, c) for a fixed cost c'.
    """
    num_costs_prime, num_distances, num_costs = distribution.shape
    ncols = min(4, num_costs_prime)
    nrows = (num_costs_prime + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
    for cost_prime in range(num_costs_prime):
        row, col = divmod(cost_prime, ncols)
        ax = axes[row][col]
        for d in range(num_distances):
            ax.plot(range(num_costs), distribution[cost_prime, d, :], label=f'd={d}')
        ax.set_xlabel('Cost c')
        ax.set_ylabel(f'N(c\'={cost_prime}; d, c)')
        ax.set_title(f'c\'={cost_prime}')
        ax.legend(fontsize='small')
    # Hide unused subplots
    for idx in range(num_costs_prime, nrows * ncols):
        row, col = divmod(idx, ncols)
        fig.delaxes(axes[row][col])
    fig.suptitle(suptitle)
    #plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # leave space for suptitle


def plot_stddev_div_mean_heatmap(distributions, cost_prime):
    """
    Plot a heatmap of the mean divided by the standard deviation of multiple
    3D distributions.
    distributions: list of np.ndarray, each shape (c', d, c)
    cost: int, the cost to fix c' at for the heatmap (since we can only plot 2D)
    """
    mean, stddev = distributions_mean_and_stddev(distributions)
    ratio = stddev / mean
    ## Copilot suggestion, don't think it's necessary
    #with np.errstate(divide='ignore', invalid='ignore'):
    #    ratio = np.where(stddev != 0, mean / stddev, 0)
    plt.imshow(ratio, aspect='auto', origin='lower')
    plt.colorbar(label='Dev/Avg')
    plt.xlabel("Cost c")
    plt.ylabel("Hamming Distance d")
    plt.title("Mean divided by Stddev Heatmap")
    plt.title(f'Dev/Avg N from cost {cost_prime}')
    plt.show()

def distribution_array_to_dict(distribution_array):
    """
    Given a 3D distribution array (e.g. n(x; d, c) or N(c'; d, c)), return a
    the array in dictionary format. In this format the keys of the dictionary
    are integer-tuples of the form (x, d, c) or (c', d, c).

    If the value of the distribution for a certain index (x, d, c) is 0, then
    no key-value pair is created in the dictionary.

    This is done to make it easier to work with distributions for multiple
    graphs, since the 
    """
    assert len(distribution_array.shape) == 3
    # Keys are tuples of 3 integers, values are integers
    distribution_dict = dict()
    for indices in np.ndindex(distribution_array.shape):
        distribution_value = distribution_array[indices]
        if distribution_value != 0:
            distribution_dict[indices] = distribution_value

    return distribution_dict

def distribution_mean_squared_error(proxy, homodist):
    """
    Given a proxy object and a homogeneous distribution array homodist (e.g. one
    computed by averaging over bitstrings), compute the mean-squared error between
    the two distributions.

    Currently only works for python proxies, not julia ones
    """
    num_constraints = homodist.shape[0] - 1

    predicted = get_homogeneous_distribution_from_proxy(proxy, num_constraints)

    ## Normalize both to sum to 1 (or same scale)
    #predicted /= predicted.sum()
    #homodist_norm = homodist / homodist.sum()

    # Compute MSE
    mse = np.mean((predicted - homodist)**2)

    return mse

def get_homogeneous_distribution_from_proxy(proxy, num_constraints=0):
    """
    Given a proxy object, compute the homogeneous distribution array N(c'; d, c)
    from the proxy. Optionally, provide num_constraints to pad the output array
    to have that many costs (in case the proxy has fewer costs than desired).

    Currently only works for python proxies, not julia ones
    """
    num_constraints = max(num_constraints, proxy.num_constraints)

    distribution = np.zeros((num_constraints+1, proxy.num_qubits+1, num_constraints+1))

    # Loop over all cost_1, distance, cost_2
    #note: this is parallelizable
    for cost_1 in range(num_constraints+1):
        for distance in range(proxy.num_qubits+1):
            for cost_2 in range(num_constraints+1):
                distribution[cost_1, distance, cost_2] = proxy.N_cost_distance_distribution(cost_1, distance, cost_2)

    return distribution

def get_pearson_correlation_coefficients(homodist1, homodist2):
    """
    Given two homogeneous distribution arrays, compute the Pearson correlation
    coefficients between each pair of 2D distributions N(c'; :, :), where the
    2D distributions are "raveled" into 1D arrays.
    """
    assert homodist1.ndim == 3 and homodist2.ndim == 3, "Both distributions must be 3D arrays"
    homodist1, homodist2 = pad_to_match(homodist1, homodist2)
    pearson_function = lambda arr1, arr2: np.corrcoef(arr1.flatten(), arr2.flatten())[0, 1]
    return [pearson_function(homodist1[i,:,:], homodist2[i,:,:]) for i in range(homodist1.shape[0])]
