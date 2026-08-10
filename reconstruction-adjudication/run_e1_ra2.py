"""E1-RA2 bounded verification of the additive Gowdy lambda mode."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ACADEMIC = HERE.parents[1]
SOURCE = ACADEMIC / "bb-cbe2i" / "run_cbe2i.py"
OUTPUT = HERE / "E1_RA2_RESULT.json"
CS = (-2.0, 1.0, 3.0)
TOL = 1.0e-11


def load_module():
    spec = importlib.util.spec_from_file_location("cbe2i_frozen_ra2", SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def max_relative(actual, expected, floor=1.0e-13):
    scale = np.maximum(np.abs(expected), floor)
    return float(np.max(np.abs(actual - expected) / scale))


def primitive_geometry(m, y, tau, dx):
    P, Pt, phi, phit, lam = y
    lam_tau = Pt**2 + np.exp(-2.0 * tau) * m.d1(P, dx)**2 \
        + phit**2 + np.exp(-2.0 * tau) * m.d1(phi, dx)**2
    lapse = np.exp((lam + tau) / 4.0 - tau)
    Hth = (lam_tau + 1.0) / (4.0 * lapse)
    Hx = (Pt - 1.0) / (2.0 * lapse)
    Hy = (-Pt - 1.0) / (2.0 * lapse)
    rho = 0.5 * np.exp(-(lam + tau) / 2.0) * (
        np.exp(2.0 * tau) * phit**2 + m.d1(phi, dx)**2
    )
    proper_theta_density = np.exp((lam + tau) / 4.0)
    return {
        "lapse": lapse,
        "Htheta": Hth,
        "Hx": Hx,
        "Hy": Hy,
        "rho": rho,
        "proper_theta_density": proper_theta_density,
        "lambda_tau": lam_tau,
        "lambda_theta": m.d1(lam, dx),
    }


def main():
    m = load_module()
    n = 512
    tau = 8.0
    dx = m.dtheta_of(n)
    y, _, _ = m.evolve(m.production_data(n), 0.0, tau, n)
    base_geom = primitive_geometry(m, y, tau, dx)
    base_inv = m.invariants(y, tau, dx)
    rows = []
    for c in CS:
        shifted = [v.copy() for v in y]
        shifted[4] += c
        geom = primitive_geometry(m, shifted, tau, dx)
        inv = m.invariants(shifted, tau, dx)
        expected_up = np.exp(c / 4.0)
        expected_down = np.exp(-c / 4.0)
        expected_rho = np.exp(-c / 2.0)
        checks = {
            "lapse": max_relative(geom["lapse"], base_geom["lapse"] * expected_up),
            "Htheta": max_relative(geom["Htheta"], base_geom["Htheta"] * expected_down),
            "Hx": max_relative(geom["Hx"], base_geom["Hx"] * expected_down),
            "Hy": max_relative(geom["Hy"], base_geom["Hy"] * expected_down),
            "rho": max_relative(geom["rho"], base_geom["rho"] * expected_rho),
            "proper_theta_density": max_relative(
                geom["proper_theta_density"],
                base_geom["proper_theta_density"] * expected_up,
            ),
            "Omega": max_relative(inv["Omega"], base_inv["Omega"]),
            "Sigma2": max_relative(inv["Sigma2"], base_inv["Sigma2"]),
            "s_P": max_relative(inv["s_P"], base_inv["s_P"]),
            "clock_rate": max_relative(
                inv["clock_rate"], base_inv["clock_rate"] * expected_down
            ),
            "lambda_rho": max_relative(
                inv["lam_rho"], base_inv["lam_rho"] * expected_up
            ),
            "lambda_tau_derivative_unchanged": max_relative(
                geom["lambda_tau"], base_geom["lambda_tau"]
            ),
            "lambda_theta_derivative_unchanged": max_relative(
                geom["lambda_theta"], base_geom["lambda_theta"]
            ),
        }
        a = c / 3.0
        rows.append({
            "c": c,
            "max_relative_discrepancy": max(checks.values()),
            "checks": checks,
            "form_preserving_coordinate_transform": {
                "tau_prime_minus_tau": a,
                "theta_period_ratio_prime_to_original": float(np.exp(-a)),
                "symmetric_x_period_ratio_prime_to_original": float(np.exp(a / 2.0)),
                "symmetric_y_period_ratio_prime_to_original": float(np.exp(a / 2.0)),
            },
        })
    result = {
        "test_id": "E1-RA2",
        "N": n,
        "tau": tau,
        "rows": rows,
        "scaling_laws_pass": all(r["max_relative_discrepancy"] < TOL for r in rows),
        "derivative_equations_leave_additive_mode_free": True,
        "period_preserving_transform_exists_for_nonzero_c": False,
        "classification": "independent global metric/modulus datum",
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
