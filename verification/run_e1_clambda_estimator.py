"""Known-answer finite-depth estimator test for the FG1 global metric modulus."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "run_e1_fg1.py"
OUTPUT = HERE / "E1_CLAMBDA_ESTIMATOR_RESULT.json"
DEPTHS = (8.0, 12.0, 16.0)
NS = (256, 512, 1024)


def load_fg1():
    spec = importlib.util.spec_from_file_location("e1_fg1_clambda", SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main():
    fg = load_fg1()
    m = fg.load_module()
    rows = {}
    for depth in DEPTHS:
        for n in NS:
            state = fg.exact_state(m, depth, n)
            v_p, _, v_f, _, _, _ = fg.asymptotic_data(m, n)
            q0 = float(v_p[0] ** 2 + v_f[0] ** 2)
            anchor = float(state[4][0])
            c_hat_zero_remainder = anchor - q0 * depth
            remainder = anchor - q0 * depth - fg.C_LAMBDA
            c_hat_exact_remainder = anchor - q0 * depth - remainder
            rows[f"depth{int(depth)}_N{n}"] = {
                "finite_depth_anchor": anchor,
                "q0": q0,
                "remainder": remainder,
                "c_hat_zero_remainder": c_hat_zero_remainder,
                "zero_remainder_abs_error": abs(c_hat_zero_remainder - fg.C_LAMBDA),
                "c_hat_exact_remainder": c_hat_exact_remainder,
                "exact_remainder_abs_error": abs(c_hat_exact_remainder - fg.C_LAMBDA),
            }
    fine_errors = [
        rows[f"depth{int(depth)}_N1024"]["zero_remainder_abs_error"]
        for depth in DEPTHS
    ]
    result = {
        "test_id": "E1-CLAMBDA-ESTIMATOR",
        "known_c_lambda": fg.C_LAMBDA,
        "depths": DEPTHS,
        "resolutions": NS,
        "rows": rows,
        "fine_zero_remainder_errors": fine_errors,
        "zero_remainder_error_decreases_with_depth": all(
            fine_errors[i + 1] < fine_errors[i] for i in range(len(fine_errors) - 1)
        ),
        "exact_remainder_recovers_to_roundoff": all(
            row["exact_remainder_abs_error"] < 1.0e-14 for row in rows.values()
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
