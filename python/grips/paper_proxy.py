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
    """
    def __init__(self, num_constraints, num_qubits, prob_edge=0.5):
        self.num_constraints = num_constraints
        self.num_qubits = num_qubits
        self.prob_edge = prob_edge

    # P(c') from paper
    def P_cost_distribution(self, cost: int) -> float:
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
        return (math.comb(self.num_qubits, distance) / p_cost) * my_sum
