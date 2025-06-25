import scipy, numpy as np
from scipy.stats import multivariate_normal, norm

class NormalProxy:
    """
    Class for implementing the "Normal" parameterized proxy for QAOA.

    Required arguments:
    - num_constraints: int
    - num_qubits: int
    - cost_mean: float
    - cov_1: float
    - cov_2: float

    For a graph with M edges and N vertices, we create a NormalProxy like this:
        proxy = NormalProxy(M, N, cost_mean, cov_1, cov_2)

    ---------------------------------------------------------------------------

    The idea seems to be to make N(c'; d, c) a multivariate normal distribution
    using variables d and c.
    (https://en.wikipedia.org/wiki/Multivariate_normal_distribution)

    The mean for the variable c is a parameter we choose (`cost_mean`). 

    We also choose the parameters cov_1 and cov_2, which are used to construct
    a diagonal matrix, representing the variance of two variables. This
    diagonal matrix is scaled and rotated by `cov_mat = P*D*P^-1`.

    TODO: Nadav, could you explain what is happening here and why?
    """
    def __init__(self, num_constraints: int, num_qubits: int, cost_mean: float,
                 cov_1: float, cov_2: float):
        self.num_constraints = num_constraints
        self.num_qubits = num_qubits
        self.cost_mean = cost_mean
        self.cov_1 = cov_1
        self.cov_2 = cov_2

    # P(c') from paper
    def P_cost_distribution(self, cost: int) -> float:
        prob_cost_mean = self.num_qubits / 2
        prob_cost_cov = self.num_qubits / 4
        return norm.pdf(cost, loc = prob_cost_mean, scale = prob_cost_cov)

    # N(c') from paper
    def N_cost_distribution(self, cost: int) -> float:
        scale = 1 << self.num_qubits
        return self.P_cost_distribution(cost) * scale


    # N(c'; d, c) from paper
    def N_cost_distance_distribution(self, cost_1: int, distance: int, cost_2: int) -> float:
        distance_mean = self.num_qubits / 2

        P = np.matrix([[cost_1 - self.cost_mean, distance_mean],
                       [-distance_mean, cost_1 - self.cost_mean]])
        P_inv = scipy.linalg.inv(P)
        cov_mat = P@np.matrix([[self.cov_1, 0], [0, self.cov_2]])@P_inv
        cov_mat[0, 1] = cov_mat[1, 0] # cov_mat must be symmetric and is prone to floating point error
        scale = (1 << self.num_qubits)
        return multivariate_normal([self.cost_mean, distance_mean], cov_mat).pdf([cost_2, distance])*scale

