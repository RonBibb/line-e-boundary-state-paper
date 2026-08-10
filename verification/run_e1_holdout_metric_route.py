"""Bounded production holdout test using an asymptotic metric route.

The test preserves the finite-depth anchor at theta=0 and changes only the
spatial construction of lambda. It does not claim independent extraction of
the limiting c_lambda modulus.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ACADEMIC = HERE.parents[1]
SOURCE = ACADEMIC / "bb-cbe2i" / "run_cbe2i.py"
OUTPUT = HERE / "E1_HOLDOUT_METRIC_ROUTE_RESULT.json"
DEPTHS = (16.0, 24.0)
NS = (256, 512, 1024)
SURFACES = (2.0, 4.0, 8.0, 12.0)


def load_model():
    spec = importlib.util.spec_from_file_location("cbe2i_metric_route", SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def compare(m, a, b, tau, dx):
    ia, ib = m.invariants(a, tau, dx), m.invariants(b, tau, dx)
    keys = ("Omega", "Sigma2", "s_P", "clock_rate")
    return float(np.max(np.sqrt(sum((ia[k] - ib[k]) ** 2 for k in keys))))


def main():
    m = load_model()
    rows = {}
    for depth in DEPTHS:
        for n in NS:
            dx = m.dtheta_of(n)
            stored_times = tuple(sorted(set(SURFACES + (depth - 4.0, depth - 2.0, depth))))
            end, _, col = m.evolve(m.production_data(n), 0.0, depth, n,
                                   store_taus=stored_times)
            col[depth] = [f.copy() for f in end]
            v_p = (col[depth - 2.0][0] - col[depth - 4.0][0]) / 2.0
            v_f = (col[depth - 2.0][2] - col[depth - 4.0][2]) / 2.0
            p_inf = col[depth - 2.0][0] - v_p * (depth - 2.0)
            f_inf = col[depth - 2.0][2] - v_f * (depth - 2.0)
            p_end = v_p * depth + p_inf
            f_end = v_f * depth + f_inf
            q = v_p**2 + v_f**2

            lam_inf_theta = 2.0 * (
                v_p * m.d1(p_inf, dx) + v_f * m.d1(f_inf, dx)
            )
            lam_inf = m.periodic_antiderivative(lam_inf_theta, dx)
            c_hat_no_remainder = float(col[depth][4][0] - q[0] * depth)
            lam_asym = q * depth + lam_inf + c_hat_no_remainder
            y0 = [p_end, v_p, f_end, v_f, lam_asym]
            _, _, rec = m.evolve(y0, depth, 0.0, n, store_taus=SURFACES)
            surface_d = {
                str(t): compare(m, col[t], rec[t], t, dx) for t in SURFACES
            }
            rows[f"depth{int(depth)}_N{n}"] = {
                "D_by_surface": surface_d,
                "D_max": max(surface_d.values()),
                "c_hat_no_remainder": c_hat_no_remainder,
                "initial_lambda_error": float(np.max(np.abs(lam_asym - col[depth][4]))),
                "initial_momentum_residual": m.momentum_residual(y0, depth, dx),
            }

    result = {
        "test_id": "E1-HOLDOUT-METRIC-ROUTE",
        "interpretation": (
            "Changes only the spatial metric route while preserving the finite-depth "
            "theta=0 anchor; c_hat_no_remainder is not an independent c_lambda estimate."
        ),
        "rows": rows,
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
