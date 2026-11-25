"""BB-CBE-3b: almost-everywhere + centered-spike test (Line E).

Executes the frozen specification
Discussion/BB_CBE3B_ALMOST_EVERYWHERE_AND_CENTERED_SPIKE_SPEC_2026-08-04.md
verbatim: launch controls C0-C3, one stored-history production run per
resolution, resolution-independent gradient-energy footprint with the
zero-measure gate, ordinary Fuchsian gates O1-O4 (forward prediction only,
no backward integration), the centered Lim-family hypothesis test S1-S3
(training tau=12,14,16; tau=18 untouched validation), the Section-8
convergence-to-limit predicate, and the Section-9 outcome mapping.
"""

from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, differential_evolution, least_squares, minimize

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "bb-cbe3"))
sys.path.insert(0, str(ROOT.parent / "bb-cbe3t"))
from run_cbe3 import (comparison_D, d1, evolve, hyperbolic_distance,  # noqa: E402
                      momentum_residual, spike_components)
import run_cbe3 as cbe3  # noqa: E402
import run_cbe3t as cbe3t  # noqa: E402
import run_cbe3t_r1 as r1  # noqa: E402

OUTPUT = ROOT / "CBE3B_RESULTS.json"
RESOLUTIONS = (1024, 2048, 4096)
STORE_TAUS = (8.0, 10.0, 12.0, 14.0, 16.0, 18.0)
FOOT_TAUS = (12.0, 14.0, 16.0, 18.0)
PARITY_GATE = 1.0e-9
CONSTRAINT_GATE = 1.0e-3
SEED = 20260804
DX_FINE = 2.0 * np.pi / 4096


def grid(n):
    dx = 2.0 * np.pi / n
    return (np.arange(n) - n // 2) * dx, dx


def production_data(n):
    x, _ = grid(n)
    return [np.zeros(n), 2.0 + 1.5 * np.cos(x), np.cos(x),
            np.zeros(n), np.zeros(n)]


def parity_asym(f):
    n = len(f)
    j = np.arange(n)
    mirror = (n - (j - n // 2) + n // 2) % n   # x -> -x on centered grid
    return float(np.max(np.abs(f - f[mirror])) / (np.max(np.abs(f)) + 1e-30))


def richardson(x1024, x2048, x4096):
    e1 = (4.0 * x2048 - x1024) / 3.0
    e2 = (4.0 * x4096 - x2048) / 3.0
    closes = abs(e2 - e1) <= max(0.2 * abs(e2), 5.0e-4)
    if not closes:
        return {"e1": e1, "e2": e2, "closes": False, "class": "intermediate"}
    cls = ("zero" if e2 <= 1.0e-3 else
           "nonzero" if e2 >= 1.0e-2 else "intermediate")
    return {"e1": e1, "e2": e2, "closes": True, "class": cls}


# ---------------------------------------------------------------- campaign --

def run_campaign():
    out = {}
    for n in RESOLUTIONS:
        t0 = time.monotonic()
        y0 = production_data(n)
        _, dx = grid(n)
        y_end, _, stored = evolve(y0, 0.0, 18.0, n, store_taus=STORE_TAUS)
        stored[18.0] = y_end
        ledger = {}
        for tau in STORE_TAUS:
            y = stored[tau]
            finite = all(np.all(np.isfinite(f)) for f in y)
            par = max(parity_asym(y[0]), parity_asym(y[2]),
                      parity_asym(y[1]), parity_asym(y[4]))
            mom = momentum_residual(y, tau, dx)
            ledger[str(tau)] = {"finite": finite, "parity": par,
                                "momentum_residual": mom}
            if not finite:
                return out, f"3b-D — nonfinite field at N={n}, tau={tau}"
            if par > PARITY_GATE:
                return out, f"3b-D-SYMMETRY — parity {par:.2e} at N={n}, tau={tau}"
        np.savez_compressed(ROOT / f"cbe3b_slices_N{n}.npz",
                            **{f"tau{tau}_{i}": stored[tau][i]
                               for tau in STORE_TAUS for i in range(5)})
        out[n] = {"stored": stored, "ledger": ledger,
                  "wall_seconds": time.monotonic() - t0}
    fine_mom = out[4096]["ledger"]["18.0"]["momentum_residual"]
    if fine_mom > CONSTRAINT_GATE:
        return out, f"3b-D — fine constraint residual {fine_mom:.2e}"
    orders = [np.log2(out[1024]["ledger"]["18.0"]["momentum_residual"]
                      / out[2048]["ledger"]["18.0"]["momentum_residual"]),
              np.log2(out[2048]["ledger"]["18.0"]["momentum_residual"]
                      / out[4096]["ledger"]["18.0"]["momentum_residual"])]
    if not all(o > 1.5 for o in orders):
        return out, f"3b-D — constraint order lost: {orders}"
    return out, None


# ---------------------------------------------------------------- controls --

def family_with_offset(w, q0, q2):
    def fn(x, tau):
        P, Q, lam = cbe3t.true_spike(x, tau, w=w, q0=q0)
        return P, Q + q2, lam
    return fn


def recover_parameters(x, dx, P_f, Q_f, P_w, tau_f=6.0):
    v = (P_f - P_w) / 2.0
    comps = spike_components(v, Q_f, dx)
    comps = [c for c in comps if c["kind"] == "TRUE"
             and abs(x[c["center"]]) <= np.pi - 1.2]
    if len(comps) != 2:
        return None, comps
    x0 = float(np.mean([x[c["center"]] for c in comps]))
    cc = int(np.argmin(np.abs(x - x0)))
    v_c = float(v[cc])

    def wm(w):
        return (2.0 - (np.log(np.cosh(w * tau_f))
                       - np.log(np.cosh(w * (tau_f - 2.0)))) / 2.0 - v_c)
    w_hat = float(brentq(wm, 1.0e-3, 1.999))
    q0_hat = 0.5 * np.exp(2.0 * tau_f - np.log(np.cosh(w_hat * tau_f))
                          - float(P_f[cc]))
    q2_hat = float(Q_f[cc]) + q0_hat * w_hat * np.exp(-2.0 * tau_f)
    return {"x0": x0, "w": w_hat, "q0": q0_hat, "q2": q2_hat}, comps


def control_C0():
    conv = cbe3.control_exact_unpolarized(0.9, 0.0, "regression")
    analytic = r1.analytic_residuals()
    amended_ok = analytic.get("passes", False)
    return {"transformed_convergence": conv, "analytic": analytic,
            "passes": bool(conv["passes"] and amended_ok)}


def control_C1():
    cases = [(0.35, 0.30, 0.17), (0.85, 0.80, -0.11)]
    results = []
    n = 4096
    x, dx = grid(n)
    for (w, q0, q2) in cases:
        fn = family_with_offset(w, q0, q2)
        P_f, Q_f, _ = fn(x, 6.0)
        P_w, _, _ = fn(x, 4.0)
        rec, comps = recover_parameters(x, dx, P_f, Q_f, P_w)
        row = {"case": [w, q0, q2], "recovered": rec,
               "n_true_comps": len(comps)}
        ok = rec is not None and (
            abs(rec["x0"]) <= 2 * dx
            and abs(rec["w"] - w) <= 5.0e-4
            and abs(rec["q0"] - q0) / q0 <= 1.0e-3
            and abs(rec["q2"] - q2) <= 1.0e-4)
        if ok:
            fit_fn = family_with_offset(rec["w"], rec["q0"], rec["q2"])
            y_f = cbe3t.exact_state(fit_fn, x, 6.0)
            _, _, back = evolve(y_f, 6.0, 0.0, n, store_taus=(2.0, 4.0))
            ds = {}
            for tc in (2.0, 4.0):
                Pex, Qex, _ = fn(x, tc)
                dH = hyperbolic_distance(back[tc][0], back[tc][2], Pex, Qex)
                ds[str(tc)] = float(np.max(dH[np.abs(x) <= np.pi - 1.2]))
            row["reconstruction_D"] = ds
            ok = all(v <= 1.0e-2 for v in ds.values())
        row["passes"] = bool(ok)
        results.append(row)
    return {"cases": results, "passes": all(r["passes"] for r in results)}


def footprint_radius(y, tau, dx, x):
    P, Pt, Q, Qt, lam = y
    G = np.exp(-2.0 * tau) * (d1(P, dx)**2 + np.exp(2.0 * P) * d1(Q, dx)**2)
    Gp = np.maximum(G - np.median(G), 0.0)
    total = float(np.sum(Gp) * dx)
    if total < 1.0e-14:
        return 0.0
    order = np.argsort(np.abs(x))
    csum = np.cumsum(Gp[order]) * dx
    idx = int(np.searchsorted(csum, 0.999 * total))
    idx = min(idx, len(x) - 1)
    return float(np.abs(x[order][idx]))


def zero_measure_gate(radii_by_n):
    """radii_by_n[n] = {tau: r_F}. Fit r = A e^{-g tau} + r_inf, all >= 0."""
    def fit(taus, rs):
        def resid(p):
            A, g, ri = p
            return A * np.exp(-g * np.asarray(taus)) + ri - np.asarray(rs)
        sol = least_squares(resid, x0=[max(rs[0], 1e-3), 0.3, 1e-4],
                            bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        return sol.x
    report, ok = {}, True
    for n, rd in radii_by_n.items():
        taus = sorted(rd)
        rs = [rd[t] for t in taus]
        A, g, ri = fit(taus, rs)
        loo = []
        for drop in range(len(taus)):
            t2 = [t for i, t in enumerate(taus) if i != drop]
            r2 = [r for i, r in enumerate(rs) if i != drop]
            loo.append(float(fit(t2, r2)[1]))
        report[str(n)] = {"A": float(A), "gamma": float(g),
                          "r_inf": float(ri), "loo_gammas": loo}
        ok = ok and g > 0 and all(l > 0 for l in loo)
    ok = ok and report["4096"]["r_inf"] <= 2.0 * DX_FINE
    rich_ok = True
    for t in FOOT_TAUS:
        r1_ = (4 * radii_by_n[2048][t] - radii_by_n[1024][t]) / 3
        r2_ = (4 * radii_by_n[4096][t] - radii_by_n[2048][t]) / 3
        if abs(r2_ - r1_) > max(0.2 * abs(r2_), 2.0 * DX_FINE):
            rich_ok = False
        report.setdefault("richardson", {})[str(t)] = [float(r1_), float(r2_)]
    return {"report": report, "passes": bool(ok and rich_ok)}


def control_C2():
    x, dx = grid(4096)
    fn_t = cbe3t.true_spike
    P_f, Q_f, _ = fn_t(x, 6.0)
    P_w, _, _ = fn_t(x, 4.0)
    ct = spike_components((P_f - P_w) / 2.0, Q_f, dx)
    ct = [c for c in ct if abs(x[c["center"]]) <= np.pi - 1.2]
    fn_f = cbe3t.false_spike
    P2, Q2_, _ = fn_f(x, 6.0)
    P2w, _, _ = fn_f(x, 4.0)
    cf = spike_components((P2 - P2w) / 2.0, Q2_, dx)
    cf = [c for c in cf if abs(x[c["center"]]) <= np.pi - 1.2]
    ck = spike_components(0.55 * np.ones_like(x), np.zeros_like(x), dx)
    Pr, Qr, _ = cbe3t.rotated_kasner(x, 6.0)
    Prw, Qrw, _ = cbe3t.rotated_kasner(x, 4.0)
    cr = spike_components((Pr - Prw) / 2.0, Qr, dx)
    taxonomy_ok = (len(ct) == 2 and all(c["kind"] == "TRUE" for c in ct)
                   and len(cf) >= 1 and all(c["kind"] == "TRANSIENT"
                                            for c in cf)
                   and len(ck) == 0 and len(cr) == 0)
    # synthetic broad nonshrinking feature must FAIL the zero-measure gate:
    # bump amplitude e^{tau} cancels the e^{-2tau} weight in the gradient
    # density, so the footprint radius is tau-independent (~ bump width 1.0)
    synth = {}
    for n in RESOLUTIONS:
        xs, dxs = grid(n)
        synth[n] = {}
        for t in FOOT_TAUS:
            bump = np.exp(-(xs / 1.0)**2)
            # P = 0.5 tau (so e^{2P} = e^{tau}), Q = e^{tau/2} bump:
            # G = e^{-2tau} e^{2P} Q_theta^2 = bump'^2, tau-independent
            y = [0.5 * t * np.ones_like(xs), 0.5 * np.ones_like(xs),
                 np.exp(t / 2.0) * bump, np.zeros_like(xs),
                 np.zeros_like(xs)]
            synth[n][t] = footprint_radius(y, t, dxs, xs)
    gate = zero_measure_gate(synth)
    return {"taxonomy_ok": bool(taxonomy_ok),
            "synthetic_gate_failed_as_required": bool(not gate["passes"]),
            "synthetic_radii_fine": synth[4096],
            "passes": bool(taxonomy_ok and not gate["passes"])}


def control_C3():
    seqs = {
        "second_order_to_zero": (4.0e-3, 1.0e-3, 2.5e-4, "zero"),
        "second_order_to_0.02": (0.024, 0.021, 0.02025, "nonzero"),
        "noisy_from_below": (0.0185, 0.0197, 0.0199, "nonzero"),
    }
    rows, ok = {}, True
    for name, (a, b, c, want) in seqs.items():
        got = richardson(a, b, c)["class"]
        rows[name] = {"got": got, "want": want}
        ok = ok and got == want
    return {"rows": rows, "passes": bool(ok)}


# ---------------------------------------------------- footprints and masks --

def compute_footprints(campaign):
    radii = {}
    for n in RESOLUTIONS:
        x, dx = grid(n)
        radii[n] = {t: footprint_radius(campaign[n]["stored"][t], t, dx, x)
                    for t in FOOT_TAUS}
    return radii


def window_velocity(stored, a, b):
    return (stored[b][0] - stored[a][0]) / (b - a)


def ordinary_mask_fine(campaign, radii):
    x, dx = grid(4096)
    R = max(radii[4096][t] for t in FOOT_TAUS) + abs(np.exp(-8.0)
                                                     - np.exp(-18.0))
    mask = np.abs(x) > R
    v = window_velocity(campaign[4096]["stored"], 16.0, 18.0)
    mask &= v >= 0.05
    return mask, float(R)


def project_mask(mask_fine, n):
    f = 4096 // n
    m = mask_fine.reshape(n, f)
    return np.all(m, axis=1)


# ------------------------------------------------------------ O1-O4 gates --

def v_hyp_of(y):
    P, Pt, Q, Qt, lam = y
    return np.sqrt(Pt**2 + np.exp(2.0 * P) * Qt**2)


def gate_O1(campaign, masks):
    rows, excl_measure = {}, {}
    for n in RESOLUTIONS:
        m = masks[n]
        for t in (14.0, 16.0, 18.0):
            vh = v_hyp_of(campaign[n]["stored"][t])[m]
            bad = np.sum((vh < 0.05) | (vh > 0.95))
            rows[f"N{n}_tau{t}"] = int(bad)
            excl_measure.setdefault(n, {})[t] = float(
                bad * (2 * np.pi / n))
    all_zero = all(v == 0 for v in rows.values())
    return {"violations": rows, "passes": bool(all_zero),
            "note": "zero violating cells at every resolution and time"
            if all_zero else "nonzero exclusions; measure gate required"}


def spatial_forcing(y, tau, dx):
    P, Pt, Q, Qt, lam = y
    e2 = np.exp(-2.0 * tau)
    Pth, Qth = d1(P, dx), d1(Q, dx)
    Pthth = (np.roll(P, -1) - 2 * P + np.roll(P, 1)) / dx**2
    Qthth = (np.roll(Q, -1) - 2 * Q + np.roll(Q, 1)) / dx**2
    SP = e2 * (Pthth - np.exp(2 * P) * Qth**2)
    SQ = e2 * (Qthth + 2 * Pth * Qth)
    Ptt = e2 * Pthth + np.exp(2 * P) * (Qt**2 - e2 * Qth**2)
    Qtt = e2 * Qthth - 2 * (Pt * Qt - e2 * Pth * Qth)
    num = np.sqrt(SP**2 + np.exp(2 * P) * SQ**2)
    den = 1.0 + np.sqrt(Ptt**2 + np.exp(2 * P) * Qtt**2)
    return num / den


def gate_O2(campaign, masks):
    sups = {}
    for n in RESOLUTIONS:
        _, dx = grid(n)
        m = masks[n]
        sups[n] = {t: float(np.max(spatial_forcing(
            campaign[n]["stored"][t], t, dx)[m])) for t in (14.0, 16.0, 18.0)}
    slopes = {n: float(np.polyfit([14, 16, 18],
                                  np.log([sups[n][t] for t in (14.0, 16.0,
                                                               18.0)]), 1)[0])
              for n in RESOLUTIONS}
    rich = richardson(sups[1024][18.0], sups[2048][18.0], sups[4096][18.0])
    passes = all(s < 0 for s in slopes.values()) and rich["class"] == "zero"
    return {"sup_epsS": {str(k): v for k, v in sups.items()},
            "decay_slopes": {str(k): v for k, v in slopes.items()},
            "richardson_18": rich, "passes": bool(passes)}


def extract_window(stored, a, b):
    P_b, Pt_b, Q_b, Qt_b, _ = stored[b]
    v = (stored[b][0] - stored[a][0]) / (b - a)
    Pinf = P_b - v * b
    piQ = np.exp(2.0 * P_b) * Qt_b
    v_safe = np.where(np.abs(v) < 0.05, 0.05, v)
    Qinf = Q_b + np.exp(-2.0 * P_b) * piQ / (2.0 * v_safe)
    return {"v": v, "Pinf": Pinf, "Qinf": Qinf, "piQ": piQ}


def gate_O3(campaign, masks):
    changes = {}
    for n in RESOLUTIONS:
        st = campaign[n]["stored"]
        m = masks[n]
        wins = [extract_window(st, a, b)
                for a, b in ((12.0, 14.0), (14.0, 16.0), (16.0, 18.0))]
        deltas = []
        for w1, w2 in zip(wins[:-1], wins[1:]):
            d = max(float(np.max(np.abs(w2[k][m] - w1[k][m])
                                 / (1.0 + np.abs(w2[k][m]))))
                    for k in ("v", "Pinf", "Qinf", "piQ"))
            deltas.append(d)
        changes[n] = deltas
    fine_ok = changes[4096][-1] <= 1.0e-3
    rich = richardson(changes[1024][-1], changes[2048][-1],
                      changes[4096][-1])
    # robustness: hyperbolic-velocity variant
    hyp_changes = {}
    for n in RESOLUTIONS:
        st = campaign[n]["stored"]
        m = masks[n]
        vh = [np.sqrt(((st[b][0] - st[a][0]) / 2.0)**2
                      + np.exp(2 * st[b][0]) * st[b][3]**2)
              for a, b in ((14.0, 16.0), (16.0, 18.0))]
        hyp_changes[n] = float(np.max(np.abs(vh[1][m] - vh[0][m])
                                      / (1.0 + np.abs(vh[1][m]))))
    hyp_ok = hyp_changes[4096] <= 1.0e-3
    passes = fine_ok and rich["class"] == "zero" and hyp_ok
    return {"window_changes": {str(k): v for k, v in changes.items()},
            "richardson": rich,
            "hyperbolic_variant": {str(k): v for k, v in hyp_changes.items()},
            "passes": bool(passes)}


def gate_O4(campaign, masks):
    per = {}
    for n in RESOLUTIONS:
        st = campaign[n]["stored"]
        _, dx = grid(n)
        m = masks[n]
        w = extract_window(st, 12.0, 14.0)
        per[n] = {}
        for tp in (16.0, 18.0):
            P_pred = w["v"] * tp + w["Pinf"]
            Qt_pred = np.exp(-2.0 * P_pred) * w["piQ"]
            v_safe = np.where(np.abs(w["v"]) < 0.05, 0.05, w["v"])
            Q_pred = w["Qinf"] - np.exp(-2.0 * P_pred) * w["piQ"] / (2 * v_safe)
            pred = [P_pred, w["v"].copy(), Q_pred, Qt_pred,
                    np.zeros_like(P_pred)]
            D = comparison_D(st[tp], pred, tp, dx)
            per[n][tp] = {"sup": float(np.max(D[m])),
                          "l2": float(np.sqrt(np.mean(D[m]**2)))}
    rich_sup = richardson(per[1024][18.0]["sup"], per[2048][18.0]["sup"],
                          per[4096][18.0]["sup"])
    rich_l2 = richardson(per[1024][18.0]["l2"], per[2048][18.0]["l2"],
                         per[4096][18.0]["l2"])
    depth_ok = per[4096][18.0]["sup"] <= per[4096][16.0]["sup"]
    passes = (rich_sup["class"] == "zero" and rich_l2["class"] == "zero"
              and depth_ok)
    converged_fail = (rich_sup["closes"] and rich_sup["e2"] >= 1.0e-2)
    return {"per": {str(k): {str(t): v for t, v in d.items()}
                    for k, d in per.items()},
            "richardson_sup": rich_sup, "richardson_l2": rich_l2,
            "depth_trend_ok": bool(depth_ok), "passes": bool(passes),
            "converged_fail": bool(converged_fail)}


# ---------------------------------------------------------- S1-S3 Lim fit --

def fit_objective_factory(campaign, n, radii):
    st = campaign[n]["stored"]
    x, dx = grid(n)

    def objective(p):
        w, lnq0, q2 = p
        fn = family_with_offset(w, np.exp(lnq0), q2)
        J = 0.0
        for t in (12.0, 14.0, 16.0):
            rF = radii[n][t]
            m = np.abs(x) <= max(rF, 2 * dx)
            P, Pt, Q, Qt, lam = st[t]
            Pf, Qf, _ = fn(x[m], t)
            Ptf, Qtf, _ = cbe3t.tau_derivative(fn, x[m], t)
            dh2 = hyperbolic_distance(P[m], Q[m], Pf, Qf)**2
            vh_p = np.sqrt(Pt[m]**2 + np.exp(2 * P[m]) * Qt[m]**2)
            vh_f = np.sqrt(Ptf**2 + np.exp(2 * Pf) * Qtf**2)
            J += float(np.mean(dh2 + 0.25 * (vh_p - vh_f)**2))
        return J
    return objective


def run_fit(campaign, radii):
    bounds = [(1.0e-3, 0.999), (-10.0, 10.0), (-10.0, 10.0)]
    fits = {}
    for n in RESOLUTIONS:
        obj = fit_objective_factory(campaign, n, radii)
        de = differential_evolution(obj, bounds, seed=SEED, popsize=20,
                                    maxiter=500, polish=False, tol=1e-12)
        ref = minimize(obj, de.x, method="L-BFGS-B", bounds=bounds)
        fits[n] = {"w": float(ref.x[0]), "q0": float(np.exp(ref.x[1])),
                   "q2": float(ref.x[2]), "J": float(ref.fun),
                   "de_J": float(de.fun), "de_nit": int(de.nit)}
    return fits


def validate_fit(campaign, radii, fits):
    n = 4096
    st = campaign[n]["stored"]
    x, dx = grid(n)
    f = fits[n]
    fn = family_with_offset(f["w"], f["q0"], f["q2"])
    rF = radii[n][18.0]
    mF = np.abs(x) <= max(rF, 2 * dx)
    Pf, Qf, _ = fn(x, 18.0)
    Ptf, Qtf, _ = cbe3t.tau_derivative(fn, x, 18.0)
    fam_state = [Pf, Ptf, Qf, Qtf, np.zeros_like(Pf)]
    D = comparison_D(st[18.0], fam_state, 18.0, dx)
    sup_D = float(np.max(D[mF]))
    P, Pt, Q, Qt, lam = st[18.0]
    dh2 = hyperbolic_distance(P[mF], Q[mF], Pf[mF], Qf[mF])**2
    vh_p = np.sqrt(Pt[mF]**2 + np.exp(2 * P[mF]) * Qt[mF]**2)
    vh_f = np.sqrt(Ptf[mF]**2 + np.exp(2 * Pf[mF]) * Qtf[mF]**2)
    profile = float(np.sqrt(np.mean(dh2 + 0.25 * (vh_p - vh_f)**2)))
    v_prod = window_velocity(st, 16.0, 18.0)
    comps_prod = [c for c in spike_components(v_prod, Q, dx)
                  if abs(x[c["center"]]) <= np.pi - 1.2]
    Pw_f, _, _ = fn(x, 16.0)
    v_fam = (Pf - Pw_f) / 2.0
    comps_fam = [c for c in spike_components(v_fam, Qf, dx)
                 if abs(x[c["center"]]) <= np.pi - 1.2]
    same_structure = (len(comps_prod) == len(comps_fam)
                      and all(a["kind"] == b["kind"] for a, b in
                              zip(comps_prod, comps_fam)))
    edge_tol = max(2 * dx, 0.1 * rF)
    edges_ok = same_structure and all(
        abs(x[a["center"]] - x[b["center"]]) <= edge_tol
        for a, b in zip(comps_prod, comps_fam))
    cc = n // 2
    vc_ok = abs(vh_p_center(st, cc) - float(np.sqrt(
        Ptf[cc]**2 + np.exp(2 * Pf[cc]) * Qtf[cc]**2))) <= 1.0e-2
    widths = {}
    widths_ok = True
    for t in (14.0, 16.0, 18.0):
        vp = window_velocity(st, t - 2.0, t)
        cp = [c for c in spike_components(vp, st[t][2], dx)
              if abs(x[c["center"]]) <= np.pi - 1.2]
        Pft, Qft, _ = fn(x, t)
        Pfw, _, _ = fn(x, t - 2.0)
        cf = [c for c in spike_components((Pft - Pfw) / 2.0, Qft, dx)
              if abs(x[c["center"]]) <= np.pi - 1.2]
        wp = sum(c["width_cells"] for c in cp) * dx
        wf = sum(c["width_cells"] for c in cf) * dx
        widths[str(t)] = [wp, wf]
        if wp == 0 or abs(wp - wf) > 0.1 * wp:
            widths_ok = False
    s2_pass = (sup_D <= 1.0e-2 and profile <= 1.0e-3 and edges_ok
               and vc_ok and widths_ok)
    ws = [fits[m]["w"] for m in RESOLUTIONS]
    q0s = [fits[m]["q0"] for m in RESOLUTIONS]
    q2s = [fits[m]["q2"] for m in RESOLUTIONS]
    def rich_pair(vals):
        a = (4 * vals[1] - vals[0]) / 3
        b = (4 * vals[2] - vals[1]) / 3
        return a, b, abs(b - a) <= 0.05 * max(abs(b), 1e-30)
    s3_rows = {k: rich_pair(v) for k, v in
               (("w", ws), ("q0", q0s), ("q2", q2s))}
    s3_params_ok = all(r[2] for r in s3_rows.values())
    sup_by_n = {}
    for m in RESOLUTIONS:
        xm, dxm = grid(m)
        stm = campaign[m]["stored"]
        fm = fits[m]
        fnm = family_with_offset(fm["w"], fm["q0"], fm["q2"])
        Pfm, Qfm, _ = fnm(xm, 18.0)
        Ptfm, Qtfm, _ = cbe3t.tau_derivative(fnm, xm, 18.0)
        Dm = comparison_D(stm[18.0], [Pfm, Ptfm, Qfm, Qtfm,
                                      np.zeros_like(Pfm)], 18.0, dxm)
        mFm = np.abs(xm) <= max(radii[m][18.0], 2 * dxm)
        sup_by_n[m] = float(np.max(Dm[mFm]))
    rich_val = richardson(sup_by_n[1024], sup_by_n[2048], sup_by_n[4096])
    s3_resid_reject = rich_val["closes"] and rich_val["e2"] >= 1.0e-2
    if s2_pass and s3_params_ok and rich_val["class"] in ("zero",
                                                          "intermediate"):
        reading = "S-LIM-A"
    elif s3_resid_reject or (not widths_ok and rich_val["closes"]) \
            or (not edges_ok and rich_val["closes"]):
        reading = "S-LIM-B"
    else:
        reading = "S-LIM-D"
    return {"sup_D_18": sup_D, "profile_error_18": profile,
            "same_structure": bool(same_structure),
            "edges_ok": bool(edges_ok), "v_center_ok": bool(vc_ok),
            "widths": widths, "widths_ok": bool(widths_ok),
            "S3_params": {k: [float(r[0]), float(r[1]), bool(r[2])]
                          for k, r in s3_rows.items()},
            "validation_sup_by_N": {str(k): v for k, v in sup_by_n.items()},
            "validation_richardson": rich_val,
            "reading": reading}


def vh_p_center(st, cc):
    P, Pt, Q, Qt, lam = st[18.0]
    return float(np.sqrt(Pt[cc]**2 + np.exp(2 * P[cc]) * Qt[cc]**2))


# ---------------------------------------------------------------- outcome --

def main():
    t_start = time.monotonic()
    controls = {"C0": control_C0(), "C1": control_C1(),
                "C2": control_C2(), "C3": control_C3()}
    result = {"program": "BB-CBE-3b", "line": "E",
              "spec": "BB_CBE3B_ALMOST_EVERYWHERE_AND_CENTERED_SPIKE_SPEC_2026-08-04.md",
              "controls": {k: v for k, v in controls.items()}}
    if not all(c["passes"] for c in controls.values()):
        result["outcome"] = "3b-D-CONTROL"
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True,
                                     default=str) + "\n")
        print(json.dumps({"outcome": "3b-D-CONTROL",
                          "controls": {k: v["passes"]
                                       for k, v in controls.items()}},
                         indent=2))
        return
    campaign, stop = run_campaign()
    result["ledgers"] = {str(n): campaign[n]["ledger"] for n in campaign}
    result["wall_seconds"] = {str(n): campaign[n]["wall_seconds"]
                              for n in campaign}
    if stop:
        result["outcome"] = stop
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True,
                                     default=str) + "\n")
        print(json.dumps({"outcome": stop}, indent=2))
        return
    radii = compute_footprints(campaign)
    result["footprints"] = {str(n): {str(t): r for t, r in d.items()}
                            for n, d in radii.items()}
    zm = zero_measure_gate(radii)
    result["zero_measure_gate"] = zm
    mask_fine, R = ordinary_mask_fine(campaign, radii)
    masks = {n: project_mask(mask_fine, n) for n in RESOLUTIONS}
    result["ordinary_mask"] = {"R": R,
                               "fraction": float(np.mean(mask_fine))}
    o1 = gate_O1(campaign, masks)
    o2 = gate_O2(campaign, masks)
    o3 = gate_O3(campaign, masks)
    o4 = gate_O4(campaign, masks)
    result.update({"O1": o1, "O2": o2, "O3": o3, "O4": o4})
    fits = run_fit(campaign, radii)
    result["S1_fits"] = {str(k): v for k, v in fits.items()}
    val = validate_fit(campaign, radii, fits)
    result["S2_S3"] = val
    ord_prereq = o1["passes"] and o2["passes"] and o3["passes"]
    ord_pass = ord_prereq and o4["passes"]
    spike_reading = val["reading"]
    if ord_pass and zm["passes"] and spike_reading == "S-LIM-A":
        outcome = "3b-EQ-POINTWISE"
    elif ord_pass and zm["passes"]:
        outcome = "3b-AE"
    elif ord_prereq and o4["converged_fail"]:
        outcome = "3b-C-ORDINARY"
    elif spike_reading == "S-LIM-A":
        outcome = "3b-S-LIM-ONLY"
    else:
        outcome = "3b-D"
    result["outcome"] = outcome
    result["peak_rss_mb"] = resource.getrusage(
        resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)
    result["total_wall_seconds"] = time.monotonic() - t_start
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True,
                                 default=str) + "\n")
    print(json.dumps({
        "outcome": outcome, "zero_measure": zm["passes"],
        "O": [o1["passes"], o2["passes"], o3["passes"], o4["passes"]],
        "spike_reading": spike_reading,
        "fits": {str(k): {kk: fits[k][kk] for kk in ("w", "q0", "q2", "J")}
                 for k in fits},
        "validation": {"sup_D_18": val["sup_D_18"],
                       "profile": val["profile_error_18"]}}, indent=2))


if __name__ == "__main__":
    main()
