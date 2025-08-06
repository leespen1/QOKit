################################################################################
#
# Test that the function "hamming_distance" is implemented correctly.
#
################################################################################
from grips.real_distribution import hamming_distance

assert hamming_distance(0b000, 0b000) == 0
assert hamming_distance(0b000, 0b001) == 1
assert hamming_distance(0b000, 0b010) == 1
assert hamming_distance(0b000, 0b011) == 2
assert hamming_distance(0b000, 0b100) == 1
assert hamming_distance(0b000, 0b101) == 2
assert hamming_distance(0b000, 0b110) == 2
assert hamming_distance(0b000, 0b111) == 3

assert hamming_distance(0b101, 0b000) == 2
assert hamming_distance(0b101, 0b001) == 1
assert hamming_distance(0b101, 0b010) == 3
assert hamming_distance(0b101, 0b011) == 2
assert hamming_distance(0b101, 0b100) == 1
assert hamming_distance(0b101, 0b101) == 0
assert hamming_distance(0b101, 0b110) == 2
assert hamming_distance(0b101, 0b111) == 1
print("All hamming distances were correct.")


