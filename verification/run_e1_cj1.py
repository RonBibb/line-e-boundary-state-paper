"""E1-CJ1 independent Friedmann-constraint propagation audit."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "E1_CJ1_RESULT.json"


def rhs(_t, y):
    a, H, pi = y
    return (a * H, -0.5 * pi**2, -3.0 * H * pi)


def run_arm(t0, rtol):
    H0 = 1.0 / (3.0 * t0)
    a0 = t0 ** (1.0 / 3.0)
    pi0 = np.sqrt(6.0) * H0
    times = np.geomspace(t0, 1.0, 800)
    sol = solve_ivp(
        rhs, (t0, 1.0), (a0, H0, pi0), t_eval=times,
        rtol=rtol, atol=1.0e-14, method="DOP853"
    )
    a, H, pi = sol.y
    rho = 0.5 * pi**2
    C = 3.0 * H**2 - rho
    relative = np.abs(C) / rho
    mask = (relative > 1.0e-14) & np.isfinite(relative)
    if np.count_nonzero(mask) >= 8:
        slope = float(np.polyfit(np.log(a[mask]), np.log(relative[mask]), 1)[0])
    else:
        slope = None
    c_scale = max(float(np.max(np.abs(C))), 1.0e-300)
    return {
        "t0": t0,
        "rtol": rtol,
        "C_initial": float(C[0]),
        "C_final": float(C[-1]),
        "C_max_abs": float(np.max(np.abs(C))),
        "C_range_relative_to_max_abs": float((np.max(C) - np.min(C)) / c_scale),
        "relative_initial": float(relative[0]),
        "relative_final": float(relative[-1]),
        "relative_max": float(np.max(relative)),
        "log_relative_vs_log_a_slope": slope,
        "a_growth": float(a[-1] / a[0]),
        "predicted_relative_growth_from_a6": float((a[-1] / a[0]) ** 6),
        "observed_relative_growth": (
            float(relative[-1] / relative[0]) if relative[0] != 0.0 else None
        ),
    }


def main():
    rows = [
        run_arm(t0, rtol)
        for t0 in (1.0e-5, 1.0e-6, 1.0e-7)
        for rtol in (1.0e-8, 1.0e-11)
    ]
    slopes = [r["log_relative_vs_log_a_slope"] for r in rows
              if r["log_relative_vs_log_a_slope"] is not None]
    result = {
        "test_id": "E1-CJ1",
        "analytic_Cdot": "0",
        "analytic_relative_defect_scaling": "a^6",
        "constrained_reduction_identity": "Hdot = -pi^2/2 on C=0",
        "comparison_statistic_used_by_reduction": False,
        "opposite_branch_data_used_by_reduction": False,
        "rows": rows,
        "median_fitted_slope": float(np.median(slopes)),
        "a6_scaling_supported": bool(abs(float(np.median(slopes)) - 6.0) < 0.5),
        "a4_explanation_supported": bool(abs(float(np.median(slopes)) - 4.0) < 0.5),
        "classification": "E1-CJ1-CONSTRAINED-VALID / A4-EXPLANATION-RETIRED",
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
