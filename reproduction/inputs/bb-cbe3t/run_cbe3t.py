"""BB-CBE-3T: exact Lim spike known-answer control battery (Line E).

Executes Discussion/BB_CBE3T_EXACT_SPIKE_CONTROL_SPEC_2026-08-04.md.
Exact solutions transcribed from Lim arXiv:0710.0628v2 eqs. (17), (22),
(28)-(30), (33)-(35); conventions verified identical to run_cbe3 E1/E2.
All time derivatives of exact fields via 8th-order FD on closed forms.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "bb-cbe3"))
from run_cbe3 import (d1, evolve, hyperbolic_distance,  # noqa: E402
                      spike_components)

OUTPUT = ROOT / "CBE3T.json"

W = 0.55
Q0 = 0.5
TAU_F = 6.0
SCORE_HALF_WIDTH = np.pi - 1.2
CONTROL_N = (1024, 2048, 4096)
T4_TAUS = (2.0, 4.0)
T3_TAUS = (1.0, 2.0, 4.0)


def grid(n):
    dx = 2.0 * np.pi / n
    return np.arange(-n // 2, n // 2) * dx, dx


def scored_mask(x):
    return np.abs(x) <= SCORE_HALF_WIDTH


# ------------------------------------------------------- exact solutions ---

def true_spike(x, tau, w=W, q0=Q0):
    """Ernst (Lim eq. 7) of the rotated Taub solution.

    T0 falsified the PDF-extracted composite constants; these forms are
    re-derived from Lim's verified transformations (6)-(7): the Q
    coefficient is -q0*w (not -4*q0*w), from integrating
    Qhat_x = 4*q0*w*(1 - w*tanh(w tau))*x with the Qhat_tau cross-check
    Qhat_tau = 2*q0*w*e^{-2tau}(1 - k^2 x^2)."""
    k = w * np.exp(tau) / np.cosh(w * tau)
    P = (2.0 * tau - np.log(np.cosh(w * tau))
         - np.log(k**2 * x**2 + 1.0) - np.log(2.0 * q0))
    Q = -q0 * w * (np.exp(-2.0 * tau)
                   + 2.0 * (w * np.tanh(w * tau) - 1.0) * x**2)
    lam_lim = (4.0 * np.log(np.cosh(w * tau))
               + 2.0 * np.log(k**2 * x**2 + 1.0) - (w**2 + 4.0) * tau)
    return P, Q, -lam_lim


def false_spike(x, tau, w=W, q0=Q0):
    """Rotation (Lim eq. 6) of the Taub solution; +ln(k^2 x^2 + 1) sign
    re-derived (the inversion places 1/(k^2 x^2 + 1) inside e^{-P})."""
    k = w * np.exp(tau) / np.cosh(w * tau)
    P = (-tau + np.log(np.cosh(w * tau))
         + np.log(k**2 * x**2 + 1.0) + np.log(2.0 * q0))
    Q = (-(1.0 / (2.0 * q0)) * w * x * np.exp(2.0 * tau)
         / np.cosh(w * tau)**2 / (k**2 * x**2 + 1.0))
    lam_lim = (2.0 * np.log(np.cosh(w * tau)) - (w**2 + 1.0) * tau)
    return P, Q, -lam_lim


def kasner(x, tau, w=W):
    return w * tau * np.ones_like(x), np.zeros_like(x), -(w**2) * tau * np.ones_like(x) * (-1.0)


def rotated_kasner(x, tau, w=W, q0=Q0):
    # Lim eq. (22) with P0 + ln Q0 = 0
    P = -np.log(1.0 / np.cosh(w * tau) / (2.0 * q0))
    Q = -(1.0 / (2.0 * q0)) * (1.0 + np.tanh(w * tau))
    return P * np.ones_like(x), Q * np.ones_like(x), (w**2) * tau * np.ones_like(x)


def tau_derivative(fn, x, tau, h=1.0e-3):
    """8th-order central FD in tau on closed forms."""
    c = np.array([1/280, -4/105, 1/5, -4/5, 0.0, 4/5, -1/5, 4/105, -1/280])
    P = np.zeros_like(x); Q = np.zeros_like(x); L = np.zeros_like(x)
    for i, ci in enumerate(c):
        if ci == 0.0:
            continue
        Pi, Qi, Li = fn(x, tau + (i - 4) * h)
        P += ci * Pi; Q += ci * Qi; L += ci * Li
    return P / h, Q / h, L / h


def exact_state(fn, x, tau):
    P, Q, lam = fn(x, tau)
    Pt, Qt, _ = tau_derivative(fn, x, tau)
    return [P, Pt, Q, Qt, lam]


# ------------------------------------------------------------------- T0 ----

def pde_residual(fn, x, tau, h):
    """Residual of E1/E2 on closed forms via 4th-order FD in tau and x."""
    dx = x[1] - x[0]

    def f_at(dt_steps):
        return fn(x, tau + dt_steps * h)

    c2 = [(-1, 16), (-2, -1)]  # 4th-order second derivative stencil parts
    P0, Q0_, _ = f_at(0)
    Ptt = (-f_at(2)[0] + 16*f_at(1)[0] - 30*P0 + 16*f_at(-1)[0] - f_at(-2)[0]) / (12*h**2)
    Qtt = (-f_at(2)[1] + 16*f_at(1)[1] - 30*Q0_ + 16*f_at(-1)[1] - f_at(-2)[1]) / (12*h**2)
    Pt, Qt, _ = tau_derivative(fn, x, tau)

    def dsp(f):
        return (np.roll(f, -2) - 8*np.roll(f, -1) + 8*np.roll(f, 1)
                - np.roll(f, 2)) / (-12.0 * dx)

    def dsp2(f):
        return (-np.roll(f, -2) + 16*np.roll(f, -1) - 30*f
                + 16*np.roll(f, 1) - np.roll(f, 2)) / (12.0 * dx**2)

    e2 = np.exp(-2.0 * tau)
    e2P = np.exp(2.0 * P0)
    r1 = Ptt - e2 * dsp2(P0) - e2P * (Qt**2 - e2 * dsp(Q0_)**2)
    r2 = Qtt - e2 * dsp2(Q0_) + 2.0 * (Pt * Qt - e2 * dsp(P0) * dsp(Q0_))
    inner = np.abs(x) < SCORE_HALF_WIDTH
    return float(np.max(np.abs(r1[inner]))), float(np.max(np.abs(r2[inner])))


def control_T0():
    x, _ = grid(4096)
    out = {}
    for name, fn in (("true_spike", true_spike), ("false_spike", false_spike)):
        rows = []
        for h in (2.0e-3, 1.0e-3, 5.0e-4):
            r1, r2 = pde_residual(fn, x, 3.0, h)
            rows.append({"h": h, "residual_P": r1, "residual_Q": r2})
        out[name] = rows
    passes = all(min(r["residual_P"] for r in rows) < 1e-8
                 and min(r["residual_Q"] for r in rows) < 1e-8
                 for rows in out.values())
    return {"residuals": out, "passes": bool(passes)}


# ------------------------------------------------------------------- T1 ----

def smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t**3 * (10.0 - 15.0 * t + 6.0 * t**2)


def wrap_taper(x, inner=np.pi - 0.4, outer=np.pi - 0.1):
    """C^2 window: 1 for |x|<=inner, 0 for |x|>=outer. Applied to the odd
    (non-periodic) Q of the false spike so the wrap jump cannot generate
    NaNs that spread at stencil speed; the taper's PHYSICAL influence stays
    outside the scored region by causality."""
    return 1.0 - smoothstep((np.abs(x) - inner) / (outer - inner))


def control_T1(fn, label, taper_q=False, half_width=SCORE_HALF_WIDTH):
    rows = []
    for n in CONTROL_N:
        x, dx = grid(n)
        y0 = exact_state(fn, x, 0.0)
        if taper_q:
            chi = wrap_taper(x)
            y0[2] = y0[2] * chi
            y0[3] = y0[3] * chi
        y, _, _ = evolve(y0, 0.0, TAU_F, n)
        Pex, Qex, _ = fn(x, TAU_F)
        dH = hyperbolic_distance(y[0], y[2], Pex, Qex)
        mask = np.abs(x) <= half_width
        rows.append({"N": n,
                     "sup_dH2_error": float(np.max(dH[mask]))})
    orders = [float(np.log2(rows[i]["sup_dH2_error"]
                            / rows[i + 1]["sup_dH2_error"])) for i in range(2)]
    passes = all(abs(o - 2.0) <= 0.3 for o in orders)
    return {"label": label, "rows": rows, "orders": orders,
            "passes": bool(passes)}


# ------------------------------------------------------------------- T2 ----

def detector_on(fn, n=4096):
    x, dx = grid(n)
    P_f, Q_f, _ = fn(x, TAU_F)
    P_w, _, _ = fn(x, TAU_F - 2.0)
    v = (P_f - P_w) / 2.0
    comps = spike_components(v, Q_f, dx)
    keep = []
    for c in comps:
        # map cells to x, keep only scored region
        xc = x[c["center"]]
        if abs(xc) <= SCORE_HALF_WIDTH:
            c["x_center"] = float(xc)
            keep.append(c)
    return keep


def control_T2():
    radius = (1.0 / W) * np.exp(-TAU_F) * np.cosh(W * TAU_F)
    results = {}
    ts = detector_on(true_spike)
    results["true_spike"] = ts
    ok_true = (len(ts) >= 1 and any(c["kind"] == "TRUE"
                                    and abs(c["x_center"]) < radius
                                    for c in ts)
               and all(c["kind"] == "TRUE" for c in ts))
    fs = detector_on(false_spike)
    results["false_spike"] = fs
    ok_false = len(fs) >= 1 and all(c["kind"] == "TRANSIENT" for c in fs)
    kn = detector_on(kasner)
    rk = detector_on(rotated_kasner)
    results["kasner"], results["rotated_kasner"] = kn, rk
    ok_none = len(kn) == 0 and len(rk) == 0
    return {"matrix": {k: [(c["kind"], round(c["x_center"], 4))
                           for c in v] for k, v in results.items()},
            "exact_radius_at_tau_f": float(radius),
            "passes": bool(ok_true and ok_false and ok_none)}


# ------------------------------------------------------------------- T3 ----

def vtd_reconstruct_and_compare(n=4096):
    x, dx = grid(n)
    P_f, Q_f, lam_f = true_spike(x, TAU_F)
    Pt_f, Qt_f, _ = tau_derivative(true_spike, x, TAU_F)
    P_w, _, _ = true_spike(x, TAU_F - 2.0)
    v = (P_f - P_w) / 2.0
    Pinf = P_f - v * TAU_F
    piQ = np.exp(2.0 * P_f) * Qt_f
    v_safe = np.where(np.abs(v) < 0.05, 0.05, v)
    Qinf = Q_f + np.exp(-2.0 * P_f) * piQ / (2.0 * v_safe)
    P0 = v * TAU_F + Pinf
    Q0_ = Qinf - np.exp(-2.0 * P0) * piQ / (2.0 * v_safe)
    Qt0 = np.exp(-2.0 * P0) * piQ
    lam_th = 2.0 * (v * d1(P0, dx) + np.exp(2 * P0) * Qt0 * d1(Q0_, dx))
    from run_cbe3 import periodic_antiderivative
    lam0, _ = periodic_antiderivative(lam_th, dx)
    lam0 += lam_f[len(x) // 4]      # anchor inside scored region
    _, _, stored = evolve([P0, v.copy(), Q0_, Qt0, lam0], TAU_F, 0.0, n,
                          store_taus=T3_TAUS)
    radius = (1.0 / W) * np.exp(-TAU_F) * np.cosh(W * TAU_F)
    rows = {}
    for tc in T3_TAUS:
        Pex, Qex, _ = true_spike(x, tc)
        dH = hyperbolic_distance(stored[tc][0], stored[tc][2], Pex, Qex)
        center = np.abs(x) < 3.0 * radius
        flank = (np.abs(x) > 10.0 * radius) & scored_mask(x)
        rows[str(tc)] = {"D_center_max": float(np.max(dH[center])),
                         "D_flank_max": float(np.max(dH[flank]))}
    return rows


# ------------------------------------------------------------------- T4 ----

def control_T4(n=4096):
    x, dx = grid(n)
    P_f, Q_f, lam_f = true_spike(x, TAU_F)
    P_w, _, _ = true_spike(x, TAU_F - 2.0)
    v = (P_f - P_w) / 2.0
    comps = spike_components(v, Q_f, dx)
    comps = [c for c in comps if abs(x[c["center"]]) <= SCORE_HALF_WIDTH
             and c["kind"] == "TRUE"]
    if comps:
        # paired TRUE components are the two edges of one spike: the
        # physical center is their midpoint (T0 battery round-1 lesson)
        x0 = float(np.mean([x[c["center"]] for c in comps]))
        center_cell = int(np.argmin(np.abs(x - x0)))
    else:
        center_cell = int(np.argmax(v))
        x0 = float(x[center_cell])
    v_center = float(v[center_cell])
    # spec Amendment A2: exact window-matched estimator for w
    from scipy.optimize import brentq

    def window_mismatch(w):
        return (2.0 - (np.log(np.cosh(w * TAU_F))
                       - np.log(np.cosh(w * (TAU_F - 2.0)))) / 2.0
                - v_center)
    w_hat = float(brentq(window_mismatch, 1.0e-3, 1.999))
    # ln(2 Q0_hat) = 2 tau_f - ln cosh(w_hat tau_f) - P(center): exact
    q0_hat = 0.5 * np.exp(2.0 * TAU_F
                          - np.log(np.cosh(w_hat * TAU_F))
                          - float(P_f[center_cell]))
    # A1-corrected family: Q(center) = -q0 w e^{-2 tau} + q2
    q2_hat = (float(Q_f[center_cell])
              + q0_hat * w_hat * np.exp(-2.0 * TAU_F))

    def family(xa, tau):
        P, Q, lam = true_spike(xa - x0, tau, w=w_hat, q0=q0_hat)
        return P, Q + q2_hat, lam

    y_f = exact_state(family, x, TAU_F)
    _, _, stored = evolve(y_f, TAU_F, 0.0, n, store_taus=T4_TAUS)
    rows = {}
    for tc in T4_TAUS:
        Pex, Qex, _ = true_spike(x, tc)
        dH = hyperbolic_distance(stored[tc][0], stored[tc][2], Pex, Qex)
        rows[str(tc)] = float(np.max(dH[scored_mask(x)]))
    passes = all(val <= 1.0e-2 for val in rows.values())
    return {"fitted": {"w_hat": w_hat, "q0_hat": q0_hat, "q2_hat": q2_hat,
                       "x0_hat": x0},
            "sup_D_scored": rows, "passes": bool(passes)}


def main():
    t0 = control_T0()
    result = {"program": "BB-CBE-3T", "line": "E",
              "spec": "BB_CBE3T_EXACT_SPIKE_CONTROL_SPEC_2026-08-04.md",
              "T0_solutionhood": t0}
    if not t0["passes"]:
        result["outcome"] = "T0 failed — transcription not certified; stop"
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"outcome": result["outcome"],
                          "T0": t0["residuals"]}, indent=2))
        return
    t1a = control_T1(true_spike, "true spike")
    t1b = control_T1(false_spike, "false spike", taper_q=True,
                     half_width=np.pi - 1.5)
    t2 = control_T2()
    result.update({"T1_true": t1a, "T1_false": t1b, "T2_confusion": t2})
    t3 = vtd_reconstruct_and_compare()
    t4 = control_T4()
    result.update({"T3_information_loss": t3, "T4_spike_aware": t4})
    detector_validated = t1a["passes"] and t1b["passes"] and t2["passes"]
    if detector_validated and t4["passes"]:
        outcome = ("CBE-3T-PASS — detector and spike taxonomy validated on "
                   "the exact family; four-function VTD loss quantified "
                   "(T3); spike-aware reconstruction restores equivalence "
                   "for the exact family (T4). CBE-3 production finding "
                   "upgraded to validated-class spike mismatch against "
                   "four-function VTD; enlarged-data question for "
                   "production spikes remains open (CBE-3b).")
    elif detector_validated:
        outcome = ("CBE-3T-PARTIAL — detector validated; spike-aware "
                   "reconstruction failed its gate; retain "
                   "CBE-D_architecture for interpretation.")
    else:
        outcome = ("CBE-3T-FAIL — detector or evolution failed the exact "
                   "control; CBE-3 wording retracted to 'localized feature "
                   "mismatch'; CBE-D retained.")
    result["outcome"] = outcome
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "T0": t0["passes"], "T1_true": t1a["passes"],
        "T1_false": t1b["passes"], "T2": t2["passes"],
        "T2_matrix": t2["matrix"], "T3": t3,
        "T4": {"passes": t4["passes"], "sup_D": t4["sup_D_scored"],
               "w_hat": t4["fitted"]["w_hat"]},
        "outcome": outcome}, indent=2))


if __name__ == "__main__":
    main()
