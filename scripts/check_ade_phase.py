"""Representative exact-algebraic E-type phase and pentagon checks."""
import json
from pathlib import Path

import mpmath as mp

from scripts.ade_phase import ade_group, phase_from_exact_quaternions
from scripts.ade_subgroups import sympy_quaternion_multiply


CASES = {
    "2T": ((1, 5, 9, 13), (2, 7, 14, 20)),
    "2O": ((1, 5, 9, 17), (3, 11, 21, 31)),
    "2I": ((1, 5, 17, 41), (3, 13, 29, 67)),
}


def product(points, i, j):
    return sympy_quaternion_multiply(points[i], points[j])


def check_group(name, index_sets):
    points = ade_group(name)
    errors = []
    phases = []
    with mp.workdps(70):
        for a_i, b_i, c_i, d_i in index_sets:
            a, b, c, d = (points[i] for i in (a_i, b_i, c_i, d_i))
            lhs = (phase_from_exact_quaternions(b, c, d, dps=65)
                   * phase_from_exact_quaternions(a, product(points, b_i, c_i), d, dps=65)
                   * phase_from_exact_quaternions(a, b, c, dps=65))
            rhs = (phase_from_exact_quaternions(product(points, a_i, b_i), c, d, dps=65)
                   * phase_from_exact_quaternions(a, b, product(points, c_i, d_i), dps=65))
            errors.append(float(abs(lhs-rhs)))
            phases.append(mp.nstr(phase_from_exact_quaternions(a, b, c, dps=65), 30))
    return {"order": len(points), "cases": len(index_sets),
            "max_pentagon_error": max(errors), "sample_phases": phases}


def main():
    result = {name: check_group(name, cases) for name, cases in CASES.items()}
    path = Path(__file__).resolve().parents[1] / "docs" / "review-0922" / "SU2_ADE_phase_validation.json"
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    assert all(item["max_pentagon_error"] < 1e-55 for item in result.values())


if __name__ == "__main__":
    main()
