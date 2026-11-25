"""BB-CBE-3: unpolarized Gowdy T^3 spike-sector equivalence (Line E).

Executes Discussion/BB_CBE3_UNPOLARIZED_SPIKE_EQUIVALENCE_SPEC_2026-08-04.md.
Controls gate the campaign; results reported per point class with the
adverse-preservation rule enforced in the classifier.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.special import j0, j1, y0 as bessel_y0, y1 as bessel_y1

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "CBE3.json"

DT_FACTOR = 0.2
COMPARISON_TAUS = (2.0, 4.0, 8.0)
CAMPAIGN_N = (1024, 2048, 4096)     # spec Amendment A1
DEPTHS = (12.0, 18.0)
D_TOL = 1.0e-3
D_MISMATCH = 1.0e-2
CONSTRAINT_GATE = 1.0e-3
SPIKE_PAD = 8
V_UNSETTLED = 0.05


def dtheta_of(n):
    return 2.0 * np.pi / n


def d1(f, dx):
    return (np.roll(f, -1) - np.roll(f, 1)) / (2.0 * dx)


def d2(f, dx):
    return (np.roll(f, -1) - 2.0 * f + np.roll(f, 1)) / dx**2


def periodic_antiderivative(f, dx):
    n = len(f)
    k = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    F = np.fft.fft(f)
    mean = float(np.real(F[0]) / n)
    F[0] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        A = np.where(k != 0.0, F / (1j * k), 0.0)
    out = np.real(np.fft.ifft(A))
    return out - out[0], mean


def rhs(tau, y, dx):
    P, Pt, Q, Qt, lam = y
    e2 = np.exp(-2.0 * tau)
    e2P = np.exp(2.0 * P)
    Pth, Qth = d1(P, dx), d1(Q, dx)
    return [
        Pt,
        e2 * d2(P, dx) + e2P * (Qt**2 - e2 * Qth**2),
        Qt,
        e2 * d2(Q, dx) - 2.0 * (Pt * Qt - e2 * Pth * Qth),
        Pt**2 + e2 * Pth**2 + e2P * (Qt**2 + e2 * Qth**2),
    ]


def rk4(y, tau, dt, dx):
    def add(a, b, s):
        return [ai + s * bi for ai, bi in zip(a, b)]
    k1 = rhs(tau, y, dx)
    k2 = rhs(tau + dt / 2, add(y, k1, dt / 2), dx)
    k3 = rhs(tau + dt / 2, add(y, k2, dt / 2), dx)
    k4 = rhs(tau + dt, add(y, k3, dt), dx)
    return [yi + dt / 6 * (a + 2 * b + 2 * c + d)
            for yi, a, b, c, d in zip(y, k1, k2, k3, k4)]


def evolve(y0, tau0, tau1, n, store_taus=(), dt_cap=None):
    dx = dtheta_of(n)
    step_scale = min(DT_FACTOR * dx, dt_cap) if dt_cap else DT_FACTOR * dx
    backward = tau1 < tau0
    targets = sorted((t for t in store_taus
                      if (tau0 < t <= tau1) or (tau1 <= t < tau0)),
                     reverse=backward)
    y, tau, stored = [f.copy() for f in y0], tau0, {}
    for t_next in targets + [tau1]:
        seg = t_next - tau
        if seg == 0.0:
            stored[t_next] = [f.copy() for f in y]
            continue
        steps = max(1, int(np.ceil(abs(seg) / step_scale)))
        dt = seg / steps
        for _ in range(steps):
            y = rk4(y, tau, dt, dx)
            tau += dt
        tau = t_next
        if t_next != tau1 or t_next in targets:
            stored[t_next] = [f.copy() for f in y]
    return y, tau, stored


def momentum_residual(y, tau, dx):
    P, Pt, Q, Qt, lam = y
    target = 2.0 * (Pt * d1(P, dx) + np.exp(2 * P) * Qt * d1(Q, dx))
    got = d1(lam, dx)
    return float(np.max(np.abs(got - target))
                 / (np.max(np.abs(target)) + 1e-30))


def hyperbolic_distance(P1, Q1, P2, Q2):
    y1_, y2_ = np.exp(-P1), np.exp(-P2)
    arg = 1.0 + ((Q1 - Q2)**2 + (y1_ - y2_)**2) / (2.0 * y1_ * y2_)
    return np.arccosh(np.maximum(arg, 1.0))


def comparison_D(ycol, yexp, tau, dx):
    Pc, Ptc, Qc, Qtc, lamc = ycol
    Pe, Pte, Qe, Qte, lame = yexp
    e2 = np.exp(-2.0 * tau)
    dH = hyperbolic_distance(Pc, Qc, Pe, Qe)
    v_c = np.sqrt(Ptc**2 + np.exp(2 * Pc) * Qtc**2)
    v_e = np.sqrt(Pte**2 + np.exp(2 * Pe) * Qte**2)
    lt_c = Ptc**2 + e2 * d1(Pc, dx)**2 + np.exp(2 * Pc) * (Qtc**2 + e2 * d1(Qc, dx)**2)
    lt_e = Pte**2 + e2 * d1(Pe, dx)**2 + np.exp(2 * Pe) * (Qte**2 + e2 * d1(Qe, dx)**2)
    return np.sqrt(dH**2 + (v_c - v_e)**2 + ((lt_c - lt_e) / (1 + np.abs(lt_c)))**2)


# ------------------------------------------------- SL(2,R) exact controls --

def polarized_seed(theta, tau, a_amp, b_amp, mode):
    """Exact polarized solution: each angular mode n pairs with Bessel
    argument n*e^-tau. P~ = a*J0(mode e^-tau)cos(mode th)
    + b*Y0(1*e^-tau)cos(th); the Y0 n=1 mode gives sign-varying asymptotic
    velocity ~ -(2/pi) b cos(th) when b != 0."""
    za = mode * np.exp(-tau)
    zb = np.exp(-tau)
    Pt_ = a_amp * j0(za) * np.cos(mode * theta) \
        + b_amp * bessel_y0(zb) * np.cos(theta)
    # d/dtau: dz/dtau = -z; dJ0/dz = -J1, dY0/dz = -Y1
    dPt = a_amp * za * j1(za) * np.cos(mode * theta) \
        + b_amp * zb * bessel_y1(zb) * np.cos(theta)
    return Pt_, dPt


def transform_c(Ptil, Ptil_tau, Ptil_th, c):
    """zeta -> zeta/(c zeta + 1) applied to polarized zeta = i e^{-Ptil}."""
    yv = np.exp(-Ptil)
    den = 1.0 + c**2 * yv**2
    Q = c * yv**2 / den
    eP = den / yv                      # e^{P'} = (1 + c^2 y^2)/y
    P = np.log(eP)
    # derivatives via y_t = -y * Ptil_t
    def dfield(dPtil):
        dy = -yv * dPtil
        dden = 2.0 * c**2 * yv * dy
        dQ = (2.0 * c * yv * dy * den - c * yv**2 * dden) / den**2
        dP = dden / den - dy / yv
        return dP, dQ
    dP_t, dQ_t = dfield(Ptil_tau)
    dP_th, dQ_th = dfield(Ptil_th)
    return P, dP_t, dP_th, Q, dQ_t, dQ_th


def exact_transformed_state(n, tau, a_amp, b_amp, mode, c):
    th = np.arange(n) * dtheta_of(n)
    Ptil, Ptil_tau = polarized_seed(th, tau, a_amp, b_amp, mode)
    Ptil_th = -mode * a_amp * j0(mode * np.exp(-tau)) * np.sin(mode * th) \
        - b_amp * bessel_y0(np.exp(-tau)) * np.sin(th)
    P, Pt, Pth, Q, Qt, Qth = transform_c(Ptil, Ptil_tau, Ptil_th, c)
    return [P, Pt, Q, Qt, np.zeros(n)], th


def control_exact_unpolarized(a_amp, b_amp, label):
    rows = []
    for n in (128, 256, 512):
        y0_, th = exact_transformed_state(n, 0.0, a_amp, b_amp, 3, 1.0)
        y, _, _ = evolve(y0_, 0.0, 8.0, n)
        yex, _ = exact_transformed_state(n, 8.0, a_amp, b_amp, 3, 1.0)
        err = float(np.max(hyperbolic_distance(y[0], y[2], yex[0], yex[2])))
        rows.append({"N": n, "sup_dH2_error": err})
    orders = [float(np.log2(rows[i]["sup_dH2_error"]
                            / rows[i + 1]["sup_dH2_error"])) for i in range(2)]
    passes = all(abs(o - 2.0) <= 0.3 for o in orders)
    return {"label": label, "rows": rows, "orders": orders,
            "passes": bool(passes)}


def spike_components(v, Q, dx):
    """Spec Amendment A2: two-channel detector on O(1)-normalized fields."""
    channels = []
    for field in (v, Q):
        g = np.abs(d1(field, dx))
        thr = 10.0 * (np.median(g) + 0.01 * np.max(g))
        channels.append((g, thr))
    mask = np.zeros(len(v), dtype=bool)
    score = np.zeros(len(v))
    for g, thr in channels:
        mask |= g > thr
        score += g / thr
    comps, idx = [], np.flatnonzero(mask)
    if len(idx):
        splits = np.flatnonzero(np.diff(idx) > 1)
        comps = np.split(idx, splits + 1)
    out = []
    for comp in comps:
        center = comp[int(np.argmax(score[comp]))]
        v_center = float(v[center])
        out.append({"cells": [int(comp[0]), int(comp[-1])],
                    "width_cells": int(len(comp)),
                    "center": int(center), "v_center": v_center,
                    "v_max": float(np.max(v[comp])),
                    "kind": "TRUE" if v_center > 1.0 else "TRANSIENT"})
    return out


def control_false_spike():
    """Sign-varying-velocity seed: exact transformed solution develops
    false-spike structure at the velocity zero crossings."""
    conv = control_exact_unpolarized(0.4, 0.8, "false-spike seed")
    n = 512
    tau_f = 8.0
    y0_, th = exact_transformed_state(n, 0.0, 0.4, 0.8, 3, 1.0)
    y, _, stored = evolve(y0_, 0.0, tau_f, n, store_taus=(tau_f - 2.0,))
    v = (y[0] - stored[tau_f - 2.0][0]) / 2.0
    comps = spike_components(v, y[2], dtheta_of(n))
    detector_ok = len(comps) >= 1 and all(c["kind"] == "TRANSIENT"
                                          for c in comps)
    return {"convergence": conv,
            "detected_features": comps,
            "detector_classifies_false": bool(detector_ok),
            "passes": bool(conv["passes"] and detector_ok)}


def control_homogeneous():
    n = 64
    y0_ = [np.zeros(n), 0.8 * np.ones(n), 0.3 * np.ones(n),
           0.4 * np.ones(n), np.zeros(n)]
    y, _, _ = evolve(y0_, 0.0, 10.0, n, dt_cap=1.0e-3)
    piQ0 = np.exp(2 * y0_[0][0]) * y0_[3][0]
    piQ1 = float(np.exp(2 * y[0][0]) * y[2 + 1][0])
    drift = abs(piQ1 - piQ0) / abs(piQ0)
    yb, _, _ = evolve(y, 10.0, 0.0, n, dt_cap=1.0e-3)
    back = float(max(np.max(np.abs(yb[i] - y0_[i])) for i in range(4)))
    return {"piQ_drift": float(drift), "backward_roundtrip_error": back,
            "passes": bool(drift < 1e-10 and back < 1e-7)}


def production_data(n):
    th = np.arange(n) * dtheta_of(n)
    return [np.zeros(n), 2.0 + 1.5 * np.cos(th), np.cos(th),
            np.zeros(n), np.zeros(n)]


def control_constraint():
    rows = []
    for n in CAMPAIGN_N:
        y, _, _ = evolve(production_data(n), 0.0, 8.0, n)
        rows.append({"N": n, "residual": momentum_residual(y, 8.0,
                                                           dtheta_of(n))})
    orders = [float(np.log2(rows[i]["residual"] / rows[i + 1]["residual"]))
              for i in range(2)]
    passes = (rows[-1]["residual"] < CONSTRAINT_GATE
              and all(o > 1.5 for o in orders))
    return {"rows": rows, "orders": orders, "passes": bool(passes)}


# ---------------------------------------------------------------- campaign --

def region_masks(n, comps, v):
    spike_mask = np.zeros(n, dtype=bool)
    for c in comps:
        lo, hi = c["cells"]
        for j in range(lo - SPIKE_PAD, hi + SPIKE_PAD + 1):
            spike_mask[j % n] = True
    unsettled = v < V_UNSETTLED
    for j in np.flatnonzero(unsettled).tolist():
        for k in range(j - SPIKE_PAD, j + SPIKE_PAD + 1):
            unsettled[k % n] = True
    ordinary = ~(spike_mask | unsettled)
    return {"ordinary": ordinary, "spike": spike_mask, "unsettled": unsettled}


def run_campaign():
    out = {}
    for depth in DEPTHS:
        for n in CAMPAIGN_N:
            dx = dtheta_of(n)
            store = tuple(COMPARISON_TAUS) + (depth - 2.0,)
            y_end, _, col = evolve(production_data(n), 0.0, depth, n,
                                   store_taus=store)
            col[depth] = y_end
            P_f, Pt_f, Q_f, Qt_f, lam_f = col[depth]
            v = (P_f - col[depth - 2.0][0]) / 2.0
            Pinf = P_f - v * depth
            with np.errstate(over="ignore"):
                piQ = np.exp(2.0 * P_f) * Qt_f
            v_safe = np.where(np.abs(v) < V_UNSETTLED,
                              np.sign(v) * V_UNSETTLED + (v == 0), v)
            Qinf = Q_f + np.exp(-2.0 * P_f) * piQ / (2.0 * v_safe)
            comps = spike_components(v, Q_f, dx)
            masks = region_masks(n, comps, v)
            P0 = v * depth + Pinf
            Q0 = Qinf - np.exp(-2.0 * P0) * piQ / (2.0 * v_safe)
            Qt0 = np.exp(-2.0 * P0) * piQ
            lam_theta = 2.0 * (v * d1(P0, dx)
                               + np.exp(2 * P0) * Qt0 * d1(Q0, dx))
            lam0, mean_res = periodic_antiderivative(lam_theta, dx)
            lam0 += lam_f[0]
            y_exp = [P0, v.copy(), Q0, Qt0, lam0]
            _, _, exp = evolve(y_exp, depth, 0.0, n,
                               store_taus=COMPARISON_TAUS)
            per_surface = {}
            for tc in COMPARISON_TAUS:
                D = comparison_D(col[tc], exp[tc], tc, dx)
                per_surface[str(tc)] = {
                    reg: float(np.max(D[m])) if np.any(m) else None
                    for reg, m in masks.items()}
            out[f"depth{int(depth)}_N{n}"] = {
                "surfaces": per_surface,
                "spike_components": comps,
                "region_cells": {k: int(np.sum(m)) for k, m in masks.items()},
                "lam_mean_residual": float(mean_res),
                "momentum_residual_collapse_end":
                    momentum_residual(col[depth], depth, dx),
            }
    return out


def classify(campaign):
    fine = campaign[f"depth18_N{CAMPAIGN_N[-1]}"]
    def sup(cfg, reg):
        vals = [cfg["surfaces"][s][reg] for s in cfg["surfaces"]
                if cfg["surfaces"][s][reg] is not None]
        return max(vals) if vals else None
    d_ord = sup(fine, "ordinary")
    d_spk = sup(fine, "spike")
    spikes = fine["spike_components"]
    true_spikes = [c for c in spikes if c["kind"] == "TRUE"]
    widths_ok = all(c["width_cells"] >= 4 for c in true_spikes)
    ord_trend = all(
        sup(campaign[f"depth18_N{n2}"], "ordinary")
        <= sup(campaign[f"depth18_N{n1}"], "ordinary")
        for n1, n2 in zip(CAMPAIGN_N[:-1], CAMPAIGN_N[1:]))
    per_class = {
        "ordinary": ("pass" if d_ord is not None and d_ord <= D_TOL
                     and ord_trend else
                     "mismatch" if d_ord is not None and d_ord > D_MISMATCH
                     and ord_trend else "intermediate"),
        "true_spike": ("absent" if not true_spikes else
                       "unresolved" if not widths_ok else
                       "mismatch" if d_spk is not None and d_spk > D_MISMATCH
                       else "pass" if d_spk is not None and d_spk <= D_TOL
                       else "intermediate"),
        "transient_spike": "present" if any(c["kind"] == "TRANSIENT"
                                            for c in spikes) else "absent",
    }
    if per_class["ordinary"] == "pass" and per_class["true_spike"] in (
            "mismatch", "unresolved"):
        agg = ("CBE-3-AE — equivalence almost everywhere; "
               f"true-spike set {per_class['true_spike']}")
    elif per_class["ordinary"] == "pass" and per_class["true_spike"] in (
            "pass", "absent"):
        agg = "CBE-3-EQ-pointwise (within detector scope)"
    elif per_class["ordinary"] == "mismatch":
        agg = "CBE-C(ordinary) — inhomogeneous mismatch away from spikes"
    else:
        agg = "intermediate/CBE-D — see per-class table"
    return {"per_class": per_class, "aggregate": agg,
            "D_ordinary_fine": d_ord, "D_spike_fine": d_spk,
            "true_spike_widths": [c["width_cells"] for c in true_spikes]}


def main():
    controls = {
        "exact_unpolarized": control_exact_unpolarized(0.9, 0.0,
                                                       "sign-fixed seed"),
        "false_spike": control_false_spike(),
        "homogeneous": control_homogeneous(),
        "constraint": control_constraint(),
    }
    gates = all(c["passes"] for c in controls.values())
    result = {"program": "BB-CBE-3", "line": "E",
              "spec": "BB_CBE3_UNPOLARIZED_SPIKE_EQUIVALENCE_SPEC_2026-08-04.md",
              "controls": controls, "controls_pass": bool(gates)}
    if gates:
        campaign = run_campaign()
        result["campaign"] = campaign
        result["classification"] = classify(campaign)
    else:
        result["classification"] = {"aggregate":
                                    "controls failed — campaign not started"}
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"controls_pass": gates,
                      "classification": result["classification"]}, indent=2))


if __name__ == "__main__":
    main()
