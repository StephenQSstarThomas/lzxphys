"""Exact rational closedness check for the eight Dic_n volume branches."""
import itertools
import json
from pathlib import Path

from scripts.ade_subgroups import dic_elements
from scripts.dic_formula import dic_volume_units


def main():
    result = {}
    for n in range(2, 6):
        elements = dic_elements(n)
        max_twice_fraction = 0
        count = 0
        for a, b, c, d in itertools.product(elements, repeat=4):
            delta = (dic_volume_units(n, b, c, d)
                     - dic_volume_units(n, a*b, c, d)
                     + dic_volume_units(n, a, b*c, d)
                     - dic_volume_units(n, a, b, c*d)
                     + dic_volume_units(n, a, b, c))
            # exp(-i*pi*delta)=1 iff delta is an even integer.
            if delta.denominator != 1 or delta.numerator % 2:
                raise AssertionError((n, a, b, c, d, delta))
            max_twice_fraction = max(max_twice_fraction, abs(delta.numerator // 2))
            count += 1
        result[str(n)] = {"order": len(elements), "quadruples": count,
                          "max_abs_integer_period": max_twice_fraction}
    path = Path(__file__).resolve().parents[1] / "docs" / "review-0922" / "SU2_Dic_table_validation.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
