"""E1-FG1 exact Bessel/Fuchsian and geometric-data benchmark."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import j0, j1, y0, y1


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
INPUTS = PAPER / "reproduction" / "inputs"
SOURCE = INPUTS / "bb-cbe2i" / "run_cbe2i.py"
OUTPUT = HERE / "E1_FG1_RESULT.json"
DEPTHS = (8.0, 12.0, 16.0)
NS = (256, 512)
OBS_TAU = 2.0
C_LAMBDA = 0.4
EULER_GAMMA = float(np.euler_gamma)

P_SPEC = {
    "constant": 0.2, "velocity": 0.35,
    "modes": ((2, "cos", 0.30, 0.12), (3, "cos", -0.15, 0.05)),
}
F_SPEC = {
    "constant": -0.1, "velocity": 0.65,
    "modes": ((1, "sin", 0.25, -0.08), (4, "sin", 0.10, 0.04)),
}


def load_module():
    spec = importlib.util.spec_from_file_location("cbe2i_frozen_fg1", SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def trig(kind, k, th):
    return np.cos(k * th) if kind == "cos" else np.sin(k * th)


def trig_d(kind, k, th):
    return -k * np.sin(k * th) if kind == "cos" else k * np.cos(k * th)


def exact_field(spec, tau, th):
    value = spec["constant"] + spec["velocity"] * tau + np.zeros_like(th)
    velocity = spec["velocity"] + np.zeros_like(th)
    theta_derivative = np.zeros_like(th)
    zfac = np.exp(-tau)
    for k, kind, A, B in spec["modes"]:
        z = k * zfac
        radial = A * j0(z) + B * y0(z)
        radial_tau = z * (A * j1(z) + B * y1(z))
        value += radial * trig(kind, k, th)
        velocity += radial_tau * trig(kind, k, th)
        theta_derivative += radial * trig_d(kind, k, th)
    return value, velocity, theta_derivative


def asymptotic_field(spec, th):
    v = spec["velocity"] + np.zeros_like(th)
    finf = spec["constant"] + np.zeros_like(th)
    vth = np.zeros_like(th)
    finfth = np.zeros_like(th)
    for k, kind, A, B in spec["modes"]:
        vk = -2.0 * B / np.pi
        ck = A + (2.0 * B / np.pi) * (np.log(k / 2.0) + EULER_GAMMA)
        v += vk * trig(kind, k, th)
        finf += ck * trig(kind, k, th)
        vth += vk * trig_d(kind, k, th)
        finfth += ck * trig_d(kind, k, th)
    return v, finf, vth, finfth


def exact_lambda(m, tau, th, dx):
    P, Pt, Pth = exact_field(P_SPEC, tau, th)
    F, Ft, Fth = exact_field(F_SPEC, tau, th)
    vP, _, _, _ = asymptotic_field(P_SPEC, th)
    vF, _, _, _ = asymptotic_field(F_SPEC, th)
    q0 = float(vP[0] ** 2 + vF[0] ** 2)

    def remainder(s):
        zero = np.array([0.0])
        _, ps, pths = exact_field(P_SPEC, s, zero)
        _, fs, fths = exact_field(F_SPEC, s, zero)
        return float(ps[0] ** 2 + np.exp(-2.0 * s) * pths[0] ** 2
                     + fs[0] ** 2 + np.exp(-2.0 * s) * fths[0] ** 2 - q0)

    correction = quad(remainder, tau, 40.0, epsabs=1.0e-13, epsrel=1.0e-13)[0]
    anchor = q0 * tau + C_LAMBDA - correction
    lam_theta = 2.0 * (Pt * Pth + Ft * Fth)
    lam = m.periodic_antiderivative(lam_theta, dx)
    lam += anchor
    return lam


def exact_state(m, tau, n):
    dx = m.dtheta_of(n)
    th = np.arange(n) * dx
    P, Pt, _ = exact_field(P_SPEC, tau, th)
    F, Ft, _ = exact_field(F_SPEC, tau, th)
    lam = exact_lambda(m, tau, th, dx)
    return [P, Pt, F, Ft, lam]


def asymptotic_data(m, n):
    dx = m.dtheta_of(n)
    th = np.arange(n) * dx
    vP, Pinf, _, _ = asymptotic_field(P_SPEC, th)
    vF, Finf, _, _ = asymptotic_field(F_SPEC, th)
    lam_inf_theta = 2.0 * (vP * m.d1(Pinf, dx) + vF * m.d1(Finf, dx))
    lam_inf = m.periodic_antiderivative(lam_inf_theta, dx)
    lam_inf += C_LAMBDA
    compatibility = float(abs(np.sum(lam_inf_theta) * dx))
    return vP, Pinf, vF, Finf, lam_inf, compatibility


def construct_arm(m, depth, n, arm):
    dx = m.dtheta_of(n)
    exact_end = exact_state(m, depth, n)
    vP_exact, Pinf_exact, vF_exact, Finf_exact, lam_inf, compatibility = asymptotic_data(m, n)
    if arm == "FUCHSIAN":
        vP, Pinf, vF, Finf = vP_exact, Pinf_exact, vF_exact, Finf_exact
        P = vP * depth + Pinf
        F = vF * depth + Finf
        q = vP**2 + vF**2
        lam = q * depth + lam_inf
    elif arm == "WINDOW":
        earlier = exact_state(m, depth - 2.0, n)
        vP = (exact_end[0] - earlier[0]) / 2.0
        vF = (exact_end[2] - earlier[2]) / 2.0
        Pinf = exact_end[0] - vP * depth
        Finf = exact_end[2] - vF * depth
        P, F = exact_end[0].copy(), exact_end[2].copy()
        lam_theta = 2.0 * (vP * m.d1(P, dx) + vF * m.d1(F, dx))
        lam = m.periodic_antiderivative(lam_theta, dx)
        lam += exact_end[4][0]
    else:
        raise ValueError(arm)
    data_errors = {
        "vP": float(np.max(np.abs(vP - vP_exact))),
        "Pinf": float(np.max(np.abs(Pinf - Pinf_exact))),
        "vphi": float(np.max(np.abs(vF - vF_exact))),
        "phiinf": float(np.max(np.abs(Finf - Finf_exact))),
    }
    return [P, vP, F, vF, lam], exact_end, data_errors, compatibility


def compare(m, y, truth, tau, n):
    dx = m.dtheta_of(n)
    inv = m.invariants(y, tau, dx)
    exact_inv = m.invariants(truth, tau, dx)
    shape_keys = ("Omega", "Sigma2", "s_P")
    full_keys = shape_keys + ("clock_rate",)
    shape_D = float(np.max(np.sqrt(sum((inv[k] - exact_inv[k])**2 for k in shape_keys))))
    full_parts = [inv[k] - exact_inv[k] for k in full_keys]
    full_parts.append(np.log(inv["lam_rho"] / exact_inv["lam_rho"]))
    full_D = float(np.max(np.sqrt(sum(v**2 for v in full_parts))))
    return {
        "P": float(np.max(np.abs(y[0] - truth[0]))),
        "Pt": float(np.max(np.abs(y[1] - truth[1]))),
        "phi": float(np.max(np.abs(y[2] - truth[2]))),
        "phit": float(np.max(np.abs(y[3] - truth[3]))),
        "lambda": float(np.max(np.abs(y[4] - truth[4]))),
        "shape_D": shape_D,
        "full_D": full_D,
        "momentum_residual": m.momentum_residual(y, tau, dx),
    }


def main():
    m = load_module()
    cases = {}
    for depth in DEPTHS:
        for n in NS:
            truth = exact_state(m, OBS_TAU, n)
            for arm in ("FUCHSIAN", "WINDOW"):
                key = f"{arm}_depth{int(depth)}_N{n}"
                y0, _, data_errors, compatibility = construct_arm(m, depth, n, arm)
                y, _, _ = m.evolve(y0, depth, OBS_TAU, n)
                cases[key] = {
                    "data_errors": data_errors,
                    "compatibility": compatibility,
                    "comparison": compare(m, y, truth, OBS_TAU, n),
                }

    gates = {}
    for arm in ("FUCHSIAN", "WINDOW"):
        fine_data = [
            max(cases[f"{arm}_depth{int(d)}_N512"]["data_errors"].values())
            for d in DEPTHS
        ]
        gates[f"{arm}_data_depth_decreases"] = bool(
            fine_data[1] < fine_data[0] and fine_data[2] < fine_data[1]
        )
        coarse = cases[f"{arm}_depth16_N256"]["comparison"]
        fine = cases[f"{arm}_depth16_N512"]["comparison"]
        gates[f"{arm}_resolution_improves"] = bool(
            fine["full_D"] < coarse["full_D"] and fine["lambda"] < coarse["lambda"]
        )
        gates[f"{arm}_fine_full_below_gate"] = fine["full_D"] < 1.0e-3
        gates[f"{arm}_constraint_below_gate"] = fine["momentum_residual"] < 1.0e-3
    gates["compatibility"] = all(c["compatibility"] < 1.0e-12 for c in cases.values())
    result = {
        "test_id": "E1-FG1",
        "cases": cases,
        "gates": gates,
        "all_pass": all(gates.values()),
        "classification": "E1-FG1-FUCHSIAN-GEOMETRIC-A" if all(gates.values())
        else "E1-FG1-NONCLOSURE",
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
