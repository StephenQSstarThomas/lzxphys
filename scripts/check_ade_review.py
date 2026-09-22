"""Exact ADE closure checks and Q8 exhaustive closedness check.

The phase convention here follows the 0922 source note (its Eq. (7.1)):
omega_A is the inverse of the previously published positive-volume helper.
The inherited E/F/T proof establishes closedness for all Dic_n; Q8 is an
exhaustive implementation check on the exact rational subgroup.
"""
import itertools
import json
from pathlib import Path

import mpmath as mp

from scripts.ade_subgroups import (
    binary_icosahedral, binary_octahedral, binary_tetrahedral,
    dic_elements, dic_quaternion_exact_q8, identity_dic,
)
from src.su2_omega import omega_quaternions, omega_quaternions_source, quaternion_multiply


def source_phase(*args, k=1, dps=70):
    """A Eq. (7.1) phase, using the existing positive-volume representative's inverse."""
    return omega_quaternions_source(*args, k=k, dps=dps)


def q8_elements():
    return dic_elements(2)


def q8_check():
    with mp.workdps(85):
        elements = q8_elements()
        quaternions = {element: dic_quaternion_exact_q8(element) for element in elements}
        phase_cache = {
            (a, b, c): source_phase(quaternions[a], quaternions[b], quaternions[c], dps=70)
            for a, b, c in itertools.product(elements, repeat=3)
        }
        phase = lambda a, b, c: phase_cache[(a, b, c)]
        max_error = mp.mpf(0)
        worst = None
        for a, b, c, d in itertools.product(elements, repeat=4):
            lhs = phase(b, c, d) * phase(a, b * c, d) * phase(a, b, c)
            rhs = phase(a * b, c, d) * phase(a, b, c * d)
            error = abs(lhs - rhs)
            if error > max_error:
                max_error = error
                worst = (a, b, c, d, lhs, rhs)
        generator = next(element for element in elements if element.r == 1 and element.epsilon == 0)
        power = identity_dic(2)
        product = mp.mpc(1)
        for _ in range(4):
            product *= phase(generator, power, generator)
            power = power * generator
        return {
            "group_order": len(elements),
            "all_quadruples": len(elements) ** 4,
            "pentagon_max_error": float(max_error),
            "worst_pentagon": [str(value) for value in worst[:4]] if worst else None,
            "worst_lhs": mp.nstr(worst[4], 40) if worst else None,
            "worst_rhs": mp.nstr(worst[5], 40) if worst else None,
            "C4_cycle_invariant": mp.nstr(product, 50),
            "C4_cycle_real": float(mp.re(product)),
            "C4_cycle_imag": float(mp.im(product)),
            "C4_expected_source_convention": mp.nstr(1j, 50),
        }


def main():
    closure = {
        "2T": len(binary_tetrahedral()),
        "2O": len(binary_octahedral()),
        "2I": len(binary_icosahedral()),
    }
    # exact_closure is intentionally imported lazily: it uses SymPy's exact algebraic field.
    from scripts.ade_subgroups import exact_closure
    import sympy as sp
    exact_closure_counts = {
        "2T": exact_closure(binary_tetrahedral()),
        "2O": exact_closure(binary_octahedral(), (sp.Symbol("s2"), 2)),
        "2I": exact_closure(binary_icosahedral(), (sp.Symbol("s5"), 5)),
    }
    dic_closure = {}
    for n in range(2, 9):
        elements = dic_elements(n)
        failures = []
        for left in elements:
            for right in elements:
                product = left * right
                if product not in elements:
                    failures.append((left, right, product))
        dic_closure[str(n)] = {"order": len(elements), "failures": len(failures)}
    result = {
        "convention": "source A Eq. (7.1): omega_A = inverse of positive-volume repository helper",
        "Dic_n": dic_closure,
        "binary_polyhedral_counts": closure,
        "binary_polyhedral_exact_closure_counts": exact_closure_counts,
        "Q8": q8_check(),
        "class_identification": {
            "theorem": "Epa-Ganter arXiv:1605.09192v1 Theorem 1.1",
            "restriction_order": "|Gamma| for every finite Gamma <= SU(2)",
            "ADE_periods": {
                "C_n": "n",
                "Dic_n": "4n",
                "2T": 24,
                "2O": 48,
                "2I": 120,
            },
        },
    }
    output = Path(__file__).resolve().parents[1] / "docs" / "review-0922" / "SU2_ADE_validation.json"
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    assert all(item["failures"] == 0 for item in dic_closure.values())
    assert closure == exact_closure_counts == {"2T": 24, "2O": 48, "2I": 120}
    assert result["Q8"]["pentagon_max_error"] < 1e-60
    assert abs(result["Q8"]["C4_cycle_real"]) < 1e-15
    assert abs(result["Q8"]["C4_cycle_imag"] - 1) < 1e-15


if __name__ == "__main__":
    main()
