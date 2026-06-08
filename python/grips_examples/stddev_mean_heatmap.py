"""Recreating figure 2 from the Parameter Setting paper."""

from grips import mean_and_stddev, get_homogeneous_distribution
from matplotlib import pyplot as plt
from networkx import erdos_renyi_graph


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
    # with np.errstate(divide='ignore', invalid='ignore'):
    #    ratio = np.where(stddev != 0, mean / stddev, 0)
    plt.imshow(ratio[cost, :, :], aspect="auto", origin="lower")
    plt.colorbar(label="Dev/Avg")
    plt.xlabel("Cost c")
    plt.ylabel("Hamming Distance d")
    plt.title("Mean divided by Stddev Heatmap")
    plt.title(f"Dev/Avg N from cost {cost}")
    plt.show()


N_vertices = 10
edge_probability = 1 / 3
N_graphs = 10
graphs = [erdos_renyi_graph(N_vertices, edge_probability, seed=i) for i in range(N_graphs)]
distributions = [get_homogeneous_distribution(g) for g in graphs]
plot_stddev_div_mean_heatmap(distributions, cost=7)
