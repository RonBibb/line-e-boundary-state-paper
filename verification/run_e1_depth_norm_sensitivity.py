"""Depth and norm-sensitivity audit for the Paper-E1 production holdout.

The audit preserves the manuscript's predictive field extraction and
asymptotic metric route.  It adds shallower extraction depths and evaluates
all cases on the common observation surfaces tau=2,4.  No production truth is
assumed; this is a cross-depth and scoring-sensitivity test.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ACADEMIC = HERE.parents[1]
SOURCE = ACADEMIC / "bb-cbe2i" / "run_cbe2i.py"
OUTPUT = HERE / "E1_DEPTH_NORM_SENSITIVITY_RESULT.json"

DEPTHS = (8.0, 12.0, 16.0, 24.0)
NS = (256, 512, 1024)
SURFACES = (2.0, 4.0)
KEYS = ("Omega", "Sigma2", "s_P", "clock_rate")
WEIGHTS = {
    "identity_l2": (1.0, 1.0, 1.0, 1.0),
    "shape_emphasis_l2": (2.0, 2.0, 2.0, 0.5),
    "clock_emphasis_l2": (0.5, 0.5, 0.5, 2.0),
}


def load_model():
    spec = importlib.util.spec_from_file_location("cbe2i_depth_norm", SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def compare(m, collapse, reconstruction, tau, dx):
    a = m.invariants(collapse, tau, dx)
    b = m.invariants(reconstruction, tau, dx)
    deltas = {key: a[key] - b[key] for key in KEYS}
    component_sup = {
        key: float(np.max(np.abs(delta))) for key, delta in deltas.items()
    }
    norms = {}
    for name, weights in WEIGHTS.items():
        pointwise = np.sqrt(sum((w * deltas[key]) ** 2
                                for key, w in zip(KEYS, weights)))
        norms[name] = float(np.max(pointwise))
    norms["component_linf"] = max(component_sup.values())
    density_gap = float(np.max(
        np.abs(a["lam_rho"] - b["lam_rho"]) / a["lam_rho"]
    ))
    return {
        "component_sup": component_sup,
        "norms": norms,
        "density_clock_relative_gap": density_gap,
    }


def run_case(m, depth, n):
    dx = m.dtheta_of(n)
    store = tuple(sorted(set(SURFACES + (depth - 4.0, depth - 2.0, depth))))
    end, _, collapse = m.evolve(
        m.production_data(n), 0.0, depth, n, store_taus=store
    )
    collapse[depth] = [f.copy() for f in end]

    t0, t1 = depth - 4.0, depth - 2.0
    v_p = (collapse[t1][0] - collapse[t0][0]) / (t1 - t0)
    v_f = (collapse[t1][2] - collapse[t0][2]) / (t1 - t0)
    p_inf = collapse[t1][0] - v_p * t1
    f_inf = collapse[t1][2] - v_f * t1
    p_end = v_p * depth + p_inf
    f_end = v_f * depth + f_inf
    q = v_p**2 + v_f**2
    lam_inf_theta = 2.0 * (
        v_p * m.d1(p_inf, dx) + v_f * m.d1(f_inf, dx)
    )
    lam_inf = m.periodic_antiderivative(lam_inf_theta, dx)
    anchor = float(collapse[depth][4][0])
    c_hat = anchor - q[0] * depth
    lam_end = q * depth + lam_inf + c_hat
    reconstructed_endpoint = [p_end, v_p, f_end, v_f, lam_end]

    _, _, reconstruction = m.evolve(
        reconstructed_endpoint, depth, 0.0, n, store_taus=SURFACES
    )
    surfaces = {
        str(tau): compare(m, collapse[tau], reconstruction[tau], tau, dx)
        for tau in SURFACES
    }
    endpoint_field_errors = {
        "P": float(np.max(np.abs(p_end - collapse[depth][0]))),
        "Pt": float(np.max(np.abs(v_p - collapse[depth][1]))),
        "phi": float(np.max(np.abs(f_end - collapse[depth][2]))),
        "phit": float(np.max(np.abs(v_f - collapse[depth][3]))),
    }
    maxima = {
        "component_sup": {
            key: max(s["component_sup"][key] for s in surfaces.values())
            for key in KEYS
        },
        "norms": {
            name: max(s["norms"][name] for s in surfaces.values())
            for name in (*WEIGHTS, "component_linf")
        },
        "density_clock_relative_gap": max(
            s["density_clock_relative_gap"] for s in surfaces.values()
        ),
    }
    return {
        "finite_depth_anchor_imported": True,
        "finite_depth_anchor": anchor,
        "c_hat_zero_remainder": c_hat,
        "endpoint_field_errors": endpoint_field_errors,
        "surfaces": surfaces,
        "maxima": maxima,
    }


def main():
    m = load_model()
    rows = {}
    for depth in DEPTHS:
        for n in NS:
            key = f"depth{int(depth)}_N{n}"
            print(f"running {key}", flush=True)
            rows[key] = run_case(m, depth, n)

    summary = {}
    for depth in DEPTHS:
        depth_key = str(int(depth))
        summary[depth_key] = {}
        for norm in (*WEIGHTS, "component_linf"):
            vals = [rows[f"depth{int(depth)}_N{n}"]["maxima"]["norms"][norm]
                    for n in NS]
            summary[depth_key][norm] = {
                "values": vals,
                "difference_ratio": ((vals[0] - vals[1]) / (vals[1] - vals[2])
                                     if vals[1] != vals[2] else None),
            }
        summary[depth_key]["fine_component_sup"] = rows[
            f"depth{int(depth)}_N1024"
        ]["maxima"]["component_sup"]
        summary[depth_key]["fine_endpoint_field_errors"] = rows[
            f"depth{int(depth)}_N1024"
        ]["endpoint_field_errors"]

    result = {
        "test_id": "E1-DEPTH-NORM-SENSITIVITY",
        "source": str(SOURCE),
        "depths": DEPTHS,
        "resolutions": NS,
        "common_observation_surfaces": SURFACES,
        "components": KEYS,
        "weights": WEIGHTS,
        "interpretation": (
            "Observable-relative production audit with imported finite-depth "
            "metric anchor; tests shallower extraction and norm sensitivity, "
            "not source-independent recovery of the limiting modulus."
        ),
        "rows": rows,
        "summary": summary,
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
