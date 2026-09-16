"""Small dependency-free utilities for paired experimental evaluation."""

from math import isfinite, sqrt
from typing import Sequence


def _pairs(predicted: Sequence[float], measured: Sequence[float]) -> tuple[list[float], list[float]]:
    if len(predicted) != len(measured) or len(predicted) < 2:
        raise ValueError("paired vectors must have equal length of at least two")
    x, y = [float(v) for v in predicted], [float(v) for v in measured]
    if not all(isfinite(v) for v in x + y):
        raise ValueError("paired vectors must be finite")
    return x, y


def pearson(predicted: Sequence[float], measured: Sequence[float]) -> float:
    x, y = _pairs(predicted, measured)
    mx, my = sum(x) / len(x), sum(y) / len(y)
    dx, dy = [v - mx for v in x], [v - my for v in y]
    den = sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    if den == 0:
        raise ValueError("correlation is undefined for a constant vector")
    return sum(a * b for a, b in zip(dx, dy)) / den


def _ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + j - 1) / 2 + 1
        for index in order[i:j]:
            ranks[index] = rank
        i = j
    return ranks


def spearman(predicted: Sequence[float], measured: Sequence[float]) -> float:
    x, y = _pairs(predicted, measured)
    return pearson(_ranks(x), _ranks(y))
