import math


def decay_constant(half_life: float) -> float:
    return math.log(2) / half_life


def decay_law(initial: float, half_life: float, time: float) -> float:
    return initial * math.exp(-decay_constant(half_life) * time)
