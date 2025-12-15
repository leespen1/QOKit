'''
Refactoring the code a bit: 
Let's make one py file where we have all the proxies, here. 
This will make incorporating the cost distribution (P(c'))
a bit easier. 

TODO later: Transfer these to the Julia implementation.
'''

from numba import njit, jit, int64, float64
from numba.experimental import jitclass
import numpy as np

##########################################################################

##########################################################################


'''TRIANGLE PROXY CLASS'''
# Removing jitclass to support CostDistribution object
#TODO: Try to reincorporate jitclass in a compatible way
#BUT: Probably the Julia implementation is what will 
#matter anyways. 

# triangle_spec = [
#     ('num_constraints', int64),
#     ('num_qubits', int64),
#     ('h_tweak_sub', float64),
#     ('hc_tweak_add', float64),
#     ('l_tweak_mul', float64),
#     ('r_tweak_mul', float64),
#     ('h_peak', float64),
#     ('center_at_h_peak', float64),
# ]

# @jitclass(triangle_spec)


class TriangleProxy:
    """
    Class for implementing the "Triangle" parameterized proxy for QAOA.

    Required arguments:
    - num_constraints: int
    - num_qubits: int
    Argument with defaults (optional):
    - h_tweak_sub
    - hc_tweak_add
    - l_tweak_mul
    - r_tweak_mul
    - cost_distribution: Optional[CostDistribution]

    For a graph with M edges and N vertices, we create a TriangleProxy like this:
        proxy = TriangleProxy(M, N, h_tweak_sub, hc_tweak_add, l_tweak_mul, r_tweak_mul)

    """
    def __init__(self, num_constraints: int, num_qubits: int,
        h_tweak_sub: float = 0, hc_tweak_add: float = 0, l_tweak_mul: float = 1,
        r_tweak_mul: float = 1, cost_distribution=None):

        self.num_constraints = num_constraints
        self.num_qubits = num_qubits
        self.h_tweak_sub = h_tweak_sub # Shifts the peak of the pyramid down (Default 0)
        self.hc_tweak_add = hc_tweak_add # Moves the cost_2 of the peak to the right (Default 0)
        self.l_tweak_mul = l_tweak_mul # Defines the (inverse of the) slope of the left side of the pyramid (Default 1)
        self.r_tweak_mul = r_tweak_mul # Defines the (inverse of the) slope of the right side of the pyramid (Default 1)
        self.cost_distribution = cost_distribution
        
        # Approximate the peak value of the paper's multinomial distribution (roughly)
        assert num_qubits >= 4, "num_qubits must be at least 4"
        h_peak = (1 << (num_qubits - 4)) - h_tweak_sub
        if h_peak < 0.0:
            print("WARNING: h_peak is negative, setting to 0")
        self.h_peak = max(h_peak, 0.0)  
        self.center_at_h_peak = num_constraints / 2 + hc_tweak_add

    #this is to simplify the opt in Sednai opt 
    def set_params(self, params):
        """
        Sets the base parameters and immediately updates the derived parameters.
        """
        # Set the base parameters from the input array
        self.h_tweak_sub = params[0]
        self.hc_tweak_add = params[1]
        self.l_tweak_mul = params[2]
        self.r_tweak_mul = params[3]

        # This logic was missing and caused the previous issues. 
        h_peak = (1 << (self.num_qubits - 4)) - self.h_tweak_sub
        self.h_peak = max(h_peak, 0.0)
        self.center_at_h_peak = self.num_constraints / 2 + self.hc_tweak_add


    # P(c') from paper
    def P_cost_distribution(self, cost: int) -> float:
        if self.cost_distribution is not None:
            return self.cost_distribution(cost)

        # Compute the normalization constant for the discrete triangle distribution
        # Sum of min(c+1, n+1-c) for c=0 to n equals:
        # - (n+2)^2 / 4 if n is even
        # - (n+1)(n+3) / 4 if n is odd
        if self.num_constraints % 2 == 0:
            normalization = (self.num_constraints + 2) ** 2 / 4
        else:
            normalization = (self.num_constraints + 1) * (self.num_constraints + 3) / 4
        
        return min(cost + 1, self.num_constraints + 1 - cost) / normalization


    # N(c') from paper
    def N_cost_distribution(self, cost: int) -> float:
        scale = 1 << self.num_qubits
        return self.P_cost_distribution(cost) * scale

    # N(c'; d, c) from paper
    def N_cost_distance_distribution(self, cost_1: int, distance: int, cost_2: int) -> float:
        # Want distance to be between 0 and self.num_qubits//2 since further distance corresponds to being near the bitwise complement (which has the same cost)
        reflected_distance = distance
        if distance > self.num_qubits // 2:
            reflected_distance = self.num_qubits - distance

        # Take the peak height at reflected_distance to be on the straight line between (0 or self.num_qubits, 1) and (self.num_qubits/2, h_peak)
        h_at_cost_2 = line_between(reflected_distance, 0, 1, self.num_qubits / 2, self.h_peak)
        # Let the peak height at reflected_distance occur where cost_2 is on the stright line between cost_1 and self.num_constraints/2
        center = line_between(reflected_distance, 0, cost_1, self.num_qubits / 2, self.center_at_h_peak)
        left = center - self.l_tweak_mul * reflected_distance - 1
        right = center + self.r_tweak_mul * reflected_distance + 1

        return triangle_value(cost_2, left, right, h_at_cost_2)

hard_coded_triangle_spec = [
    ('num_constraints', int64),
    ('num_qubits', int64),
]

@jitclass(hard_coded_triangle_spec)
class HardCodedTriangleProxy:
    """
    Class for the original, hard-coded triangle distribution.

    This class is now obselete, because the same result can be achieved using
    the parameterized triangle proxy with the default parameters.
    """
    def __init__(self, num_constraints, num_qubits):
        self.num_constraints = num_constraints
        self.num_qubits = num_qubits


    def P_cost_distribution(self, cost: int) -> float:
        """
        P(c') from paper but dumber
        """
        # Compute the normalization constant for the discrete triangle distribution
        # Sum of min(c+1, n+1-c) for c=0 to n equals:
        # - (n+2)^2 / 4 if n is even
        # - (n+1)(n+3) / 4 if n is odd
        if self.num_constraints % 2 == 0:
            normalization = (self.num_constraints + 2) ** 2 / 4
        else:
            normalization = (self.num_constraints + 1) * (self.num_constraints + 3) / 4
        
        return min(cost + 1, self.num_constraints + 1 - cost) / normalization


    def N_cost_distribution(self, cost: int) -> float:
        """
        N(c') from paper but dumber
        """

        scale = 1 << self.num_qubits
        return self.P_cost_distribution(cost) * scale

    def N_cost_distance_distribution(self, cost_1: int, distance: int, cost_2: int) -> float:
        """
        N(c'; d, c) from paper but instead of a multinomial distribution, we just approximate by a prism whose cross-sections at fixed distances are triangles
        TODO: This only works for prob_edge = 0.5
        """
        # Want distance to be between 0 and num_qubits//2 since further distance corresponds to being near the bitwise complement (which has the same cost)
        reflected_distance = distance
        if distance > self.num_qubits // 2:
            reflected_distance = self.num_qubits - distance

        # Approximate the peak value of the paper's multinomial distribution (roughly)
        h_peak = 1 << (self.num_qubits - 4)
        # Take the peak height at reflected_distance to be on the straight line between (0 or self.num_qubits, 1) and (self.num_qubits/2, h_peak)
        h_at_cost_2 = line_between(reflected_distance, 0, 1, self.num_qubits / 2, h_peak)
        # Let the peak height at reflected_distance occur where cost_2 is on the stright line between cost_1 and num_constraints/2
        center = line_between(reflected_distance, 0, cost_1, self.num_qubits / 2, self.num_constraints / 2)
        left = center - reflected_distance - 1
        right = center + reflected_distance + 1

        return triangle_value(cost_2, left, right, h_at_cost_2)



@njit
def line_between(current_time: float, start_time: float, start_value: float, end_time: float, end_value: float) -> float:
    """
    Gives the y-value at x=current_time on the line between (start_time, start_value) and (end_time, end_value)
    """
    # Goes from 0 to 1 as current_time goes from start_time to end_time
    relative_time = (current_time - start_time) / (end_time - start_time)

    # Goes from start_value to end_value as relative_time goes from 0 to 1
    return (1 - relative_time) * start_value + relative_time * end_value


@njit
def triangle_value(x: int, left: int | float, right: int | float, height: int | float) -> float:
    r"""
                /\height
               /  \
    _ _ _ left/    \right _ _ _

    Given x, returns the corresponding y value on the preceeding curve
    """
    return max(0, min(x - left, right - x) * 2 * height / (right - left))


##########################################################################

##########################################################################


'''NORMAL PROXY CLASS'''
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
    - cost_distribution: Optional[CostDistribution]

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
                 cov_1: float, cov_2: float, cost_distribution=None):
        self.num_constraints = num_constraints
        self.num_qubits = num_qubits
        self.cost_mean = cost_mean
        self.cov_1 = cov_1
        self.cov_2 = cov_2
        self.cost_distribution = cost_distribution

    #This is to simplify in the optimization in sendai_opt.py
    def set_params(self, params):
        self.cost_mean = params[0]
        self.cov_1 = params[1]
        self.cov_2 = params[2]

    # P(c') from paper
    def P_cost_distribution(self, cost: int) -> float:
        if self.cost_distribution is not None:
            return self.cost_distribution(cost)

        prob_cost_mean = self.num_qubits / 2
        prob_cost_cov = self.num_qubits / 4
        
        # Compute PDF value at this cost
        pdf_value = norm.pdf(cost, loc = prob_cost_mean, scale = prob_cost_cov)
        
        # Compute normalization factor: sum of PDF values over all discrete costs [0, num_constraints]
        # Note: costs range from 0 to num_constraints (num_edges), not num_qubits (num_vertices)
        normalization = sum(norm.pdf(c, loc = prob_cost_mean, scale = prob_cost_cov) 
                          for c in range(self.num_constraints + 1))
        
        # Return normalized probability
        return pdf_value / normalization

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




##########################################################################

##########################################################################

'''PAPER PROXY CLASS'''
import numpy as np
import math, typing, time, scipy, os
from scipy.stats import binom, multinomial


class PaperProxy:
    """
    This class implements the QAOA proxy algorithm for MaxCut from:
    https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.6.023171

    Required arguments: 
    - num_constraints: int
    - num_qubits: int
    - prob_edge: float
    - cost_distribution: Optional[CostDistribution]
    """
    def __init__(self, num_constraints, num_qubits, prob_edge=0.5, cost_distribution=None):
        self.num_constraints = num_constraints
        self.num_qubits = num_qubits
        self.prob_edge = prob_edge
        self.cost_distribution = cost_distribution

    # P(c') from paper
    def P_cost_distribution(self, cost: int) -> float:
        if self.cost_distribution is not None:
            return self.cost_distribution(cost)
        return binom.pmf(cost, self.num_constraints, self.prob_edge)


    # N(c') from paper
    def N_cost_distribution(self, cost: int) -> float:
        scale = 1 << self.num_qubits
        return self.P_cost_distribution(cost) * scale

    # P(b, c'-b, c-b | d) from paper
    def prob_common_at_distance_paper(
            self,
            common_constraints: int,
            cost_1: int, 
            distance: int,
            cost_2: int,
        ) -> float:

        prob_same = (math.comb(self.num_qubits - distance, 2) + math.comb(distance, 2)) / math.comb(self.num_qubits, 2)
        prob_neither = prob_same / 2
        prob_both = prob_neither
        prob_one = (1 - prob_neither - prob_both) / 2
        return multinomial.pmf(
            [common_constraints, cost_1 - common_constraints, cost_2 - common_constraints, self.num_constraints + common_constraints - (cost_1 + cost_2)],
            self.num_constraints,
            [prob_both, prob_one, prob_one, prob_neither],
        )

    # N(c'; d, c) from paper
    def N_cost_distance_distribution(self, cost_1: int, distance: int, cost_2: int) -> float:
        my_sum = 0
        for common_constraints in range(max(0, cost_1 + cost_2 - self.num_constraints), min(cost_1, cost_2) + 1):
            my_sum += self.prob_common_at_distance_paper(
                common_constraints, cost_1, distance, cost_2
            )

        p_cost = self.P_cost_distribution(cost_1)
        if p_cost <= 1e-9:
            return 0.0
        return (math.comb(self.num_qubits, distance) / p_cost) * my_sum



