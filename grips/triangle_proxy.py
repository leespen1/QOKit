from numba import njit, jit, int64, float64
from numba.experimental import jitclass
import numpy as np

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
def triangle_value(x: int, left: int, right: int, height: float) -> float:
    r"""
                /\height
               /  \
    _ _ _ left/    \right _ _ _

    Given x, returns the corresponding y value on the preceeding curve
    """
    return max(0, min(x - left, right - x) * 2 * height / (right - left))




hard_coded_triangle_spec = [
    ('num_constraints', int64),
    ('num_qubits', int64),
]


@jitclass(hard_coded_triangle_spec)
class HardCodedTriangleProxy:
    """
    Class for the original, hard-coded triangle distribution.
    """
    def __init__(self, num_constraints, num_qubits):
        self.num_constraints = num_constraints
        self.num_qubits = num_qubits


    def P_cost_distribution(self, cost: int) -> float:
        """
        P(c') from paper but dumber
        """
        return 4 / ((self.num_constraints + 1) ** 2) * min(cost + 1, self.num_constraints + 1 - cost)


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


triangle_spec = [
    ('num_constraints', int64),
    ('num_qubits', int64),
    ('h_tweak_sub', float64),
    ('hc_tweak_add', float64),
    ('l_tweak_mul', float64),
    ('r_tweak_mul', float64),
]

@jitclass(triangle_spec)
class TriangleProxy:
    def __init__(self, num_constraints: int, num_qubits: int,
        h_tweak_sub: float, hc_tweak_add: float, l_tweak_mul: float,
        r_tweak_mul: float):

        self.num_constraints = num_constraints
        self.num_qubits = num_qubits
        self.h_tweak_sub = h_tweak_sub # Shifts the peak of the pyramid down (Default 0)
        self.hc_tweak_add = hc_tweak_add # Moves the cost_2 of the peak to the right (Default 0)
        self.l_tweak_mul = l_tweak_mul # Defines the (inverse of the) slope of the left side of the pyramid (Default 1)
        self.r_tweak_mul = r_tweak_mul # Defines the (inverse of the) slope of the right side of the pyramid (Default 1)


    def P_cost_distribution(self, cost: int) -> float:
        """
        P(c') from paper but dumber
        """
        return 4 / ((self.num_constraints + 1) ** 2) * min(cost + 1, self.num_constraints + 1 - cost)


    def N_cost_distribution(self, cost: int) -> float:
        """
        N(c') from paper but dumber
        """

        scale = 1 << self.num_qubits
        return self.P_cost_distribution(cost) * scale

    def N_cost_distance_distribution(self, cost_1: int, distance: int, cost_2: int) -> float:
        # Want distance to be between 0 and self.num_qubits//2 since further distance corresponds to being near the bitwise complement (which has the same cost)
        reflected_distance = distance
        if distance > self.num_qubits // 2:
            reflected_distance = self.num_qubits - distance

        # Approximate the peak value of the paper's multinomial distribution (roughly)
        h_peak = (1 << (self.num_qubits - 4)) - self.h_tweak_sub
        center_at_h_peak = self.num_constraints / 2 + self.hc_tweak_add
        # Take the peak height at reflected_distance to be on the straight line between (0 or self.num_qubits, 1) and (self.num_qubits/2, h_peak)
        h_at_cost_2 = line_between(reflected_distance, 0, 1, self.num_qubits / 2, h_peak)
        # Let the peak height at reflected_distance occur where cost_2 is on the stright line between cost_1 and self.num_constraints/2
        center = line_between(reflected_distance, 0, cost_1, self.num_qubits / 2, center_at_h_peak)
        left = center - self.l_tweak_mul * reflected_distance - 1
        right = center + self.r_tweak_mul * reflected_distance + 1

        return triangle_value(cost_2, left, right, h_at_cost_2)

