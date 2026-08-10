"""Bounded FG1 extension: third resolution and metric-route error budget.

This does not overwrite the frozen E1_FG1_RESULT.json campaign.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "run_e1_fg1.py"
OUTPUT = HERE / "E1_FG1_EXTENSION_RESULT.json"
DEPTH = 16.0
NS = (256, 512, 1024)


def load_fg1():
    spec = importlib.util.spec_from_file_location("e1_fg1_extension_source", SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def window_data(fg, m, n):
    dx = m.dtheta_of(n)
    exact_end = fg.exact_state(m, DEPTH, n)
    earlier = fg.exact_state(m, DEPTH - 2.0, n)
    v_p = (exact_end[0] - earlier[0]) / 2.0
    v_f = (exact_end[2] - earlier[2]) / 2.0
    p_inf = exact_end[0] - v_p * DEPTH
    f_inf = exact_end[2] - v_f * DEPTH
    lam_inf_theta = 2.0 * (
        v_p * m.d1(p_inf, dx) + v_f * m.d1(f_inf, dx)
    )
    lam_inf = m.periodic_antiderivative(lam_inf_theta, dx)
    lam_inf += fg.C_LAMBDA
    q = v_p**2 + v_f**2
    asymptotic_metric = q * DEPTH + lam_inf

    finite_theta = 2.0 * (
        v_p * m.d1(exact_end[0], dx) + v_f * m.d1(exact_end[2], dx)
    )
    finite_metric = m.periodic_antiderivative(finite_theta, dx)
    finite_metric += exact_end[4][0]
    return exact_end, [exact_end[0], v_p, exact_end[2], v_f, asymptotic_metric], [
        exact_end[0], v_p, exact_end[2], v_f, finite_metric
    ]


def orders(rows, key):
    values = [rows[str(n)][key] for n in NS]
    return [float(np.log2(values[i] / values[i + 1])) for i in range(2)]


def main():
    fg = load_fg1()
    m = fg.load_module()
    rows = {}
    for n in NS:
        dx = m.dtheta_of(n)
        truth_observe = fg.exact_state(m, fg.OBS_TAU, n)
        exact_end, window_asym, window_finite = window_data(fg, m, n)
        exact_initial, _, _, _ = fg.construct_arm(m, DEPTH, n, "FUCHSIAN")

        arm_rows = {}
        for name, initial in (
            ("exact_fuchsian", exact_initial),
            ("window_asymptotic_metric", window_asym),
            ("window_finite_metric", window_finite),
        ):
            evolved, _, _ = m.evolve(initial, DEPTH, fg.OBS_TAU, n)
            comparison = fg.compare(m, evolved, truth_observe, fg.OBS_TAU, n)
            comparison["initial_lambda_error"] = float(
                np.max(np.abs(initial[4] - exact_end[4]))
            )
            comparison["initial_momentum_residual"] = m.momentum_residual(
                initial, DEPTH, dx
            )
            arm_rows[name] = comparison
        rows[str(n)] = arm_rows

    convergence = {}
    for arm in rows["256"]:
        convergence[arm] = {}
        for key in ("P", "Pt", "phi", "phit", "lambda", "shape_D", "full_D",
                    "momentum_residual"):
            arm_values = {str(n): rows[str(n)][arm][key] for n in NS}
            convergence[arm][key] = {
                "values": arm_values,
                "orders": [
                    float(np.log2(arm_values[str(NS[i])] / arm_values[str(NS[i + 1])]))
                    for i in range(2)
                ],
            }

    result = {
        "test_id": "E1-FG1-EXTENSION",
        "preserves_frozen_fg1": True,
        "depth": DEPTH,
        "resolutions": NS,
        "multimode_total_field_constraint_with_cross_terms": True,
        "rows": rows,
        "convergence": convergence,
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
