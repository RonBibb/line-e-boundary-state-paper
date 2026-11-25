"""E1-RA3 open-family local rank and conditioning audit."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUTS = HERE.parent / "reproduction" / "inputs"
SOURCE = INPUTS / "bb-cbe2i" / "run_cbe2i.py"
OUTPUT = HERE / "E1_RA3_RESULT.json"

N_VALUES = (256, 512)
TAU_INITIAL = 16.0
TAU_OBSERVE = 2.0
DELTA = 1.0e-4
PARAMS = ("a0", "a2", "b1", "c0", "c2", "d1", "d0", "c_lambda")
BASES = {
    "A": np.array((0.40, 0.30, 0.60, 0.80, -0.20, 0.50, 0.10, 0.50)),
    "B": np.array((0.55, -0.18, 0.35, 0.68, 0.14, -0.42, -0.20, -0.75)),
}
SHAPE_ACTIVE = (0, 1, 2, 3, 4, 5)
FULL_ACTIVE = (0, 1, 2, 3, 4, 5, 7)
NULL_SHAPE = (6, 7)
NULL_FULL = (6,)


def load_module():
    spec = importlib.util.spec_from_file_location("cbe2i_frozen_ra3", SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def boundary_state(m, p, n):
    a0, a2, b1, c0, c2, d1, d0, c_lambda = p
    dx = m.dtheta_of(n)
    th = np.arange(n) * dx
    vP = a0 + a2 * np.cos(2.0 * th)
    Pinf = b1 * np.cos(th)
    vF = c0 + c2 * np.cos(2.0 * th)
    Finf = d0 + d1 * np.sin(th)
    P = vP * TAU_INITIAL + Pinf
    phi = vF * TAU_INITIAL + Finf
    lam_theta = 2.0 * (vP * m.d1(P, dx) + vF * m.d1(phi, dx))
    lam = m.periodic_antiderivative(lam_theta, dx)
    q0 = float(vP[0] ** 2 + vF[0] ** 2)
    desired_anchor = q0 * TAU_INITIAL + c_lambda
    lam += desired_anchor
    compatibility = float(abs(np.sum(lam_theta) * dx))
    return [P, vP, phi, vF, lam], compatibility


def output_vector(m, p, n, layer):
    y0, compatibility = boundary_state(m, p, n)
    y, _, _ = m.evolve(y0, TAU_INITIAL, TAU_OBSERVE, n)
    inv = m.invariants(y, TAU_OBSERVE, m.dtheta_of(n))
    stride = n // 64
    keys = ["Omega", "Sigma2", "s_P"]
    if layer == "full":
        keys += ["clock_rate"]
    pieces = [inv[key][::stride] for key in keys]
    if layer == "full":
        pieces.append(np.log(inv["lam_rho"][::stride]))
    return np.concatenate(pieces), compatibility


def normalized_jacobian(m, p, n, layer, delta):
    base, compatibility = output_vector(m, p, n, layer)
    block_count = 3 if layer == "shape" else 5
    scales = []
    for block in np.array_split(base, block_count):
        scales.extend([max(float(np.max(np.abs(block))), 1.0e-3)] * len(block))
    scales = np.asarray(scales)
    columns = []
    compatibilities = [compatibility]
    for j in range(len(PARAMS)):
        plus = p.copy()
        minus = p.copy()
        plus[j] += delta
        minus[j] -= delta
        yp, cp = output_vector(m, plus, n, layer)
        ym, cm = output_vector(m, minus, n, layer)
        compatibilities.extend((cp, cm))
        derivative = (yp - ym) / (2.0 * delta)
        pscale = max(abs(float(p[j])), 0.1)
        columns.append(derivative * pscale / scales)
    J = np.column_stack(columns)
    return J, max(compatibilities)


def analyze(J, active, null):
    column_norms = np.linalg.norm(J, axis=0)
    Ja = J[:, active]
    singular = np.linalg.svd(Ja, compute_uv=False)
    rank = int(np.sum(singular > 1.0e-7 * singular[0]))
    condition = float(singular[0] / singular[-1])
    return {
        "column_norms": {PARAMS[i]: float(column_norms[i]) for i in range(len(PARAMS))},
        "active_all_above_gate": bool(all(column_norms[i] > 1.0e-5 for i in active)),
        "null_all_below_gate": bool(all(column_norms[i] < 1.0e-8 for i in null)),
        "singular_values": [float(v) for v in singular],
        "rank": rank,
        "condition_number": condition,
    }


def factor_agreement(a, b):
    ratio = max(a, b) / min(a, b)
    return float(ratio), bool(ratio <= 2.0)


def main():
    m = load_module()
    cases = {}
    jacobians = {}
    for base_name, p in BASES.items():
        for n in N_VALUES:
            key = f"base{base_name}_N{n}"
            cases[key] = {}
            jacobians[key] = {}
            for layer, active, null in (
                ("shape", SHAPE_ACTIVE, NULL_SHAPE),
                ("full", FULL_ACTIVE, NULL_FULL),
            ):
                J, compatibility = normalized_jacobian(m, p, n, layer, DELTA)
                jacobians[key][layer] = J
                cases[key][layer] = analyze(J, active, null)
                cases[key][layer]["max_periodicity_defect"] = compatibility

    J_half, _ = normalized_jacobian(m, BASES["A"], 256, "full", DELTA / 2.0)
    J_ref = jacobians["baseA_N256"]["full"]
    active_ref = J_ref[:, FULL_ACTIVE]
    active_half = J_half[:, FULL_ACTIVE]
    linearity = float(np.linalg.norm(active_half - active_ref) / np.linalg.norm(active_ref))

    resolution_checks = {}
    for base_name in BASES:
        resolution_checks[base_name] = {}
        for layer in ("shape", "full"):
            lo = cases[f"base{base_name}_N256"][layer]
            hi = cases[f"base{base_name}_N512"][layer]
            sv_ratio, sv_pass = factor_agreement(lo["singular_values"][-1], hi["singular_values"][-1])
            cond_ratio, cond_pass = factor_agreement(lo["condition_number"], hi["condition_number"])
            resolution_checks[base_name][layer] = {
                "smallest_singular_value_factor": sv_ratio,
                "smallest_singular_value_pass": sv_pass,
                "condition_number_factor": cond_ratio,
                "condition_number_pass": cond_pass,
            }

    expected_ranks = all(
        row[layer]["rank"] == (6 if layer == "shape" else 7)
        for row in cases.values() for layer in ("shape", "full")
    )
    direction_gates = all(
        row[layer]["active_all_above_gate"] and row[layer]["null_all_below_gate"]
        for row in cases.values() for layer in ("shape", "full")
    )
    resolution_pass = all(
        check[gate]
        for base in resolution_checks.values()
        for check in base.values()
        for gate in ("smallest_singular_value_pass", "condition_number_pass")
    )
    result = {
        "test_id": "E1-RA3",
        "cases": cases,
        "linearity_relative_difference": linearity,
        "linearity_pass": linearity <= 0.02,
        "resolution_checks": resolution_checks,
        "expected_ranks_pass": expected_ranks,
        "direction_gates_pass": direction_gates,
        "resolution_stability_pass": resolution_pass,
        "classification": (
            "E1-RA3-LOCAL-INJECTIVITY-A" if
            expected_ranks and direction_gates and resolution_pass and linearity <= 0.02
            else "E1-RA3-NONCLOSURE"
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
