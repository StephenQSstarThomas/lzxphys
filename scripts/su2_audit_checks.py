"""Reproduce the handoff's numerical baselines; not a closed-form solver."""

import itertools
import math

import numpy as np
from numpy.polynomial.legendre import leggauss


PAULI = np.array(
    [[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]],
    dtype=complex,
)
IDENTITY = np.eye(2, dtype=complex)


def exp_vector(vector):
    vector = np.asarray(vector, dtype=float)
    angle = np.linalg.norm(vector)
    return np.cos(angle) * IDENTITY + 1j * np.sinc(angle / np.pi) * np.einsum(
        "a,aij->ij", vector, PAULI
    )


def geometry(matrices):
    product = IDENTITY.copy()
    vertices = [[1.0, 0.0, 0.0, 0.0]]
    for matrix in matrices:
        product = product @ matrix
        vertices.append(
            [np.trace(product).real / 2]
            + [(np.trace(s @ product) / (2j)).real for s in PAULI]
        )
    q = np.array(vertices).T
    return np.linalg.det(q), q.T @ q


def volume_integral(matrices, order):
    determinant, gram = geometry(matrices)
    nodes, weights = leggauss(order)
    nodes, weights = (nodes + 1) / 2, weights / 2
    a, b, c = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    barycentric = np.array(
        [(1 - a) * (1 - b) * (1 - c), a, (1 - a) * b, (1 - a) * (1 - b) * c]
    )
    squared_norm = np.einsum("iabc,ij,jabc->abc", barycentric, gram, barycentric)
    return determinant * np.einsum(
        "i,j,k,ijk->",
        weights,
        weights,
        weights,
        (1 - a) ** 2 * (1 - b) / squared_norm**2,
    )


def compositions(total, count):
    if count == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, count - 1):
                yield (first,) + rest


def volume_series(matrices, cutoff):
    determinant, gram = geometry(matrices)
    edges = list(itertools.combinations(range(4), 2))
    terms = []
    for total in range(cutoff + 1):
        for powers in compositions(total, 6):
            degrees = [0] * 4
            factor = 1.0
            for (i, j), power in zip(edges, powers):
                degrees[i] += power
                degrees[j] += power
                factor *= (2 * (1 - gram[i, j])) ** power / math.factorial(power)
            terms.append(
                math.factorial(total + 1)
                * math.prod(math.factorial(degree) for degree in degrees)
                / math.factorial(2 * total + 3)
                * factor
            )
    return determinant * math.fsum(terms)


def main():
    axes = np.eye(3)
    matrices = [exp_vector(axis * angle) for axis, angle in zip(axes, [.7, .8, .9])]
    determinant, _ = geometry(matrices)
    values = [volume_integral(matrices, order) for order in (32, 64)]
    series = volume_series(matrices, 12)
    print(f"numpy={np.__version__}; float64; deterministic, no random samples")
    print(f"axis determinant={determinant:.17g}")
    print(f"axis volumes (orders 32,64)={values}")
    print(f"axis series M<=12={series:.17g}; integral gap={values[1]-series:.9g}")
    assert abs(determinant - .25220893803677946) < 1e-14
    assert max(abs(value - .07082150248008584) for value in values) < 1e-13
    assert abs(series - .07082148767853955) < 1e-13

    alpha, beta = 1.1, 2.0
    family = [
        exp_vector([alpha, 0, 0]),
        1j * (np.cos(alpha) * PAULI[1] + np.sin(alpha) * PAULI[2]),
        exp_vector([beta, 0, 0]),
    ]
    values = [volume_integral(family, order) for order in (32, 64)]
    print(f"family volumes (orders 32,64)={values}; expected={alpha*beta/2}")
    assert max(abs(value - alpha * beta / 2) for value in values) < 1e-12

    trace = sum(
        np.linalg.det(axes[list(p)])
        * np.trace((1j * PAULI[p[0]]) @ (1j * PAULI[p[1]]) @ (1j * PAULI[p[2]]))
        for p in itertools.permutations(range(3))
    )
    print(f"alternating trace of i*sigma generators={trace}; expected=12")
    assert abs(trace - 12) < 1e-14
    for epsilon in (.1, .03, .01):
        volume = volume_integral([exp_vector(epsilon * axis) for axis in axes], 24)
        print(f"epsilon={epsilon}; V/(epsilon^3/6)={volume/(epsilon**3/6):.15g}")
    print("Baseline checks passed. General closed form and global branches are not tested.")


if __name__ == "__main__":
    main()
