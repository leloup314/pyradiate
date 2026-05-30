import math


def decay_constant(half_life: float) -> float:
    return math.log(2) / half_life
