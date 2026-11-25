"""E1-RA1: bounded reconstruction-adjudication test.

Imports the frozen CBE-2I implementation read-only and evaluates extraction
and metric-anchor ablations specified in RECONSTRUCTION_ADJUDICATION_SPEC_2026-08-09.md.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUTS = HERE.parent / "reproduction" / "inputs"
SOURCE = INPUTS / "bb-cbe2i" / "run_cbe2i.py"
OUTPUT = HERE / "E1_RA1_RESULT.json"

RESOLUTIONS = (256, 512, 1024)
DEPTHS = (16.0, 24.0)
SURFACES = (2.0, 4.0, 8.0, 12.0)
ARMS = {
    "V0-A0": ("endpoint", "collapse"),
    "V0-A1": ("endpoint", "zero"),
    "V1-A0": ("holdout", "collapse"),
    "V1-A1": ("holdout", "zero"),
}


def load_frozen_module():
    spec = importlib.util.spec_from_file_location("cbe2i_frozen", SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def specified_momentum_residual(m, y, dx):
    P, Pt, phi, phit, lam = y
    target = 2.0 * (Pt * m.d1(P, dx) + phit * m.d1(phi, dx))
    got = m.d1(lam, dx)
    numerator = np.max(np.abs(got - target))
    denominator = np.max(np.abs(got) + 1.0e-30)
    return float(numerator / denominator)


def extract(col_store, depth, kind):
    if kind == "endpoint":
        t1, t0 = depth, depth - 2.0
    elif kind == "holdout":
        t1, t0 = depth - 2.0, depth - 4.0
    else:
        raise ValueError(kind)
    dt = t1 - t0
    vP = (col_store[t1][0] - col_store[t0][0]) / dt
    vF = (col_store[t1][2] - col_store[t0][2]) / dt
    Pinf = col_store[t1][0] - vP * t1
    Finf = col_store[t1][2] - vF * t1
    return vP, vF, Pinf, Finf


def construct_endpoint(m, col_store, depth, dx, extraction, anchor):
    vP, vF, Pinf, Finf = extract(col_store, depth, extraction)
    P = vP * depth + Pinf
    phi = vF * depth + Finf
    lam_theta = 2.0 * (vP * m.d1(P, dx) + vF * m.d1(phi, dx))
    lam = m.periodic_antiderivative(lam_theta, dx)
    if anchor == "collapse":
        anchor_value = float(col_store[depth][4][0])
    elif anchor == "zero":
        anchor_value = 0.0
    else:
        raise ValueError(anchor)
    lam += anchor_value
    endpoint = [P, vP.copy(), phi, vF.copy(), lam]
    truth = col_store[depth]
    errors = {
        "P": float(np.max(np.abs(P - truth[0]))),
        "Pt": float(np.max(np.abs(vP - truth[1]))),
        "phi": float(np.max(np.abs(phi - truth[2]))),
        "phit": float(np.max(np.abs(vF - truth[3]))),
        "lambda": float(np.max(np.abs(lam - truth[4]))),
        "anchor_imported": anchor == "collapse",
        "anchor_value": anchor_value,
    }
    return endpoint, errors


def compare_surface(m, collapse, expansion, tau, dx):
    a = m.invariants(collapse, tau, dx)
    b = m.invariants(expansion, tau, dx)
    keys = ("Omega", "Sigma2", "s_P", "clock_rate")
    component_sup = {key: float(np.max(np.abs(a[key] - b[key]))) for key in keys}
    common_normal_D = float(np.max(np.sqrt(sum((a[key] - b[key]) ** 2 for key in keys))))

    # Expansion-future normal is opposite to the increasing-tau normal.
    expansion_future_clock_rate = -b["clock_rate"]
    parity_transformed_clock_rate = -expansion_future_clock_rate
    parity_components = [a[key] - b[key] for key in ("Omega", "Sigma2", "s_P")]
    parity_components.append(a["clock_rate"] - parity_transformed_clock_rate)
    explicit_parity_D = float(np.max(np.sqrt(sum(c**2 for c in parity_components))))
    raw_future_clock_gap = float(np.max(np.abs(a["clock_rate"] - expansion_future_clock_rate)))
    clock_gap = float(np.max(np.abs(a["lam_rho"] - b["lam_rho"]) / a["lam_rho"]))
    return {
        "D_common_normal": common_normal_D,
        "D_explicit_parity": explicit_parity_D,
        "parity_vs_common_abs": abs(explicit_parity_D - common_normal_D),
        "raw_future_clock_gap_without_T": raw_future_clock_gap,
        "lambda_rho_relative_gap": clock_gap,
        "component_sup": component_sup,
    }


def run_case(m, n, depth):
    dx = m.dtheta_of(n)
    store_times = tuple(sorted(set(SURFACES + (depth - 4.0, depth - 2.0, depth))))
    end, _, col_store = m.evolve(m.production_data(n), 0.0, depth, n, store_taus=store_times)
    col_store[depth] = [f.copy() for f in end]
    result = {
        "collapse_constraint_implemented": m.momentum_residual(col_store[depth], depth, dx),
        "collapse_constraint_specified": specified_momentum_residual(m, col_store[depth], dx),
        "arms": {},
    }
    for arm, (extraction, anchor) in ARMS.items():
        endpoint, endpoint_errors = construct_endpoint(
            m, col_store, depth, dx, extraction, anchor
        )
        _, _, exp_store = m.evolve(endpoint, depth, 0.0, n, store_taus=SURFACES)
        surfaces = {
            str(tau): compare_surface(m, col_store[tau], exp_store[tau], tau, dx)
            for tau in SURFACES
        }
        result["arms"][arm] = {
            "endpoint_errors": endpoint_errors,
            "reconstruction_constraint_implemented": m.momentum_residual(endpoint, depth, dx),
            "reconstruction_constraint_specified": specified_momentum_residual(m, endpoint, dx),
            "surfaces": surfaces,
            "D_max": max(row["D_common_normal"] for row in surfaces.values()),
            "parity_closure_max": max(row["parity_vs_common_abs"] for row in surfaces.values()),
        }
    return result


def summarize(cases, original):
    summary = {arm: {} for arm in ARMS}
    for depth in DEPTHS:
        for n in RESOLUTIONS:
            case = cases[f"depth{int(depth)}_N{n}"]
            for arm in ARMS:
                summary[arm][f"depth{int(depth)}_N{n}"] = case["arms"][arm]["D_max"]

    baseline_checks = {}
    for depth in DEPTHS:
        for n in RESOLUTIONS:
            key = f"depth{int(depth)}_N{n}"
            original_max = max(v["D"] for v in original["campaign"][key]["surfaces"].values())
            new_max = summary["V0-A0"][key]
            tolerance = max(0.05 * original_max, 1.0e-8)
            baseline_checks[key] = {
                "original": original_max,
                "new": new_max,
                "abs_difference": abs(new_max - original_max),
                "tolerance": tolerance,
                "passes": abs(new_max - original_max) <= tolerance,
            }

    def fine_both_depths_below(arm, threshold=1.0e-3):
        return all(summary[arm][f"depth{int(d)}_N1024"] < threshold for d in DEPTHS)

    def decreases_with_resolution(arm):
        for depth in DEPTHS:
            vals = [summary[arm][f"depth{int(depth)}_N{n}"] for n in RESOLUTIONS]
            if not (vals[1] <= vals[0] and vals[2] <= vals[1]):
                return False
        return True

    parity_pass = all(
        cases[key]["arms"][arm]["parity_closure_max"] <= 1.0e-12
        for key in cases for arm in ARMS
    )
    specified_constraint_pass = all(
        cases[key]["collapse_constraint_specified"] < 1.0e-3
        for key in cases if key.endswith("N1024")
    )
    return {
        "D_max_by_arm": summary,
        "baseline_reproduction": {
            "cases": baseline_checks,
            "passes": all(v["passes"] for v in baseline_checks.values()),
        },
        "anchor_ablation_V0_A1_below_gate": fine_both_depths_below("V0-A1"),
        "holdout_V1_A0": {
            "below_gate": fine_both_depths_below("V1-A0"),
            "decreases_with_resolution": decreases_with_resolution("V1-A0"),
        },
        "strict_V1_A1": {
            "below_gate": fine_both_depths_below("V1-A1"),
            "decreases_with_resolution": decreases_with_resolution("V1-A1"),
        },
        "specified_constraint_gate_passes_fine": specified_constraint_pass,
        "explicit_parity_closes": parity_pass,
    }


def main():
    m = load_frozen_module()
    original = json.loads((INPUTS / "bb-cbe2i" / "CBE2I.json").read_text())
    cases = {}
    for depth in DEPTHS:
        for n in RESOLUTIONS:
            key = f"depth{int(depth)}_N{n}"
            print(f"running {key}", flush=True)
            cases[key] = run_case(m, n, depth)
    result = {
        "test_id": "E1-RA1",
        "source": str(SOURCE),
        "spec": "RECONSTRUCTION_ADJUDICATION_SPEC_2026-08-09.md",
        "cases": cases,
        "summary": summarize(cases, original),
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
