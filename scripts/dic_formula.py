"""Closed Dic_n (binary dihedral) volume table from the 0922 calculation.

The table uses the convention g_i = V**s_i U**r_i,
U=exp(i*pi*sigma_1/n), V=i*sigma_2, and the source-note phase
omega_A=exp(-i*k*Vol/pi). Elements supplied as DicElement are U^r V^epsilon;
for epsilon=1 the table exponent is -r modulo 2n.
"""
import mpmath as mp
from fractions import Fraction

from scripts.ade_subgroups import DicElement


def residue(value, modulus):
    return int(value) % int(modulus)


def N_n(a, b, n):
    """The carry floor((a+b)/(2n)) for residues in [0,2n)."""
    return (int(a) + int(b)) // (2 * int(n))


def table_coordinates(element):
    if not isinstance(element, DicElement):
        raise TypeError("DicElement inputs required")
    return (residue((-1 if element.epsilon else 1) * element.r, 2 * element.n), element.epsilon)


def dic_oriented_volume(n, g1, g2, g3):
    """Return the signed source-note volume for all Dic_n triples."""
    if not all(isinstance(g, DicElement) and g.n == n for g in (g1, g2, g3)):
        raise ValueError("all elements must belong to the same Dic_n")
    (r1, s1), (r2, s2), (r3, s3) = map(table_coordinates, (g1, g2, g3))
    mod = 2 * n
    bracket = lambda x: residue(x, mod)
    carry = lambda a, b: N_n(bracket(a), bracket(b), n)
    if (s1, s2, s3) == (0, 0, 0):
        return -mp.pi**2 / n * r1 * carry(r2, r3)
    if (s1, s2, s3) == (0, 0, 1):
        return mp.pi**2 / n * (bracket(r3-r1-r2) - n) * carry(r1, r2)
    if (s1, s2, s3) == (0, 1, 0):
        return mp.pi**2 / (2*n**2) * r1 * r3
    if (s1, s2, s3) == (0, 1, 1):
        return -mp.pi**2 / n * (bracket(r2-r1) - n) * carry(r1, n-r2+r3)
    if (s1, s2, s3) == (1, 0, 0):
        return -mp.pi**2 / n * r1 * carry(r2, r3)
    if (s1, s2, s3) == (1, 0, 1):
        return mp.pi**2 / (2*n**2) * bracket(n-r1-r2+r3) * r2
    if (s1, s2, s3) == (1, 1, 0):
        return mp.pi**2 / n * (r1-n) * carry(n-r1+r2, r3)
    return -mp.pi**2 / (2*n**2) * bracket(n-r1+r2) * bracket(n-r2+r3)


def dic_volume_units(n, g1, g2, g3):
    """Exact rational u with Vol=pi^2*u, for symbolic closedness checks."""
    if not all(isinstance(g, DicElement) and g.n == n for g in (g1, g2, g3)):
        raise ValueError("all elements must belong to the same Dic_n")
    (r1, s1), (r2, s2), (r3, s3) = map(table_coordinates, (g1, g2, g3))
    mod = 2 * n
    bracket = lambda x: residue(x, mod)
    carry = lambda a, b: N_n(bracket(a), bracket(b), n)
    if (s1, s2, s3) == (0, 0, 0):
        return Fraction(-r1 * carry(r2, r3), n)
    if (s1, s2, s3) == (0, 0, 1):
        return Fraction((bracket(r3-r1-r2)-n) * carry(r1, r2), n)
    if (s1, s2, s3) == (0, 1, 0):
        return Fraction(r1*r3, 2*n*n)
    if (s1, s2, s3) == (0, 1, 1):
        return Fraction(-(bracket(r2-r1)-n) * carry(r1, n-r2+r3), n)
    if (s1, s2, s3) == (1, 0, 0):
        return Fraction(-r1 * carry(r2, r3), n)
    if (s1, s2, s3) == (1, 0, 1):
        return Fraction(bracket(n-r1-r2+r3) * r2, 2*n*n)
    if (s1, s2, s3) == (1, 1, 0):
        return Fraction((r1-n) * carry(n-r1+r2, r3), n)
    return Fraction(-bracket(n-r1+r2) * bracket(n-r2+r3), 2*n*n)


def dic_phase(n, g1, g2, g3, k=1):
    return mp.exp(-1j * int(k) * dic_oriented_volume(n, g1, g2, g3) / mp.pi)
