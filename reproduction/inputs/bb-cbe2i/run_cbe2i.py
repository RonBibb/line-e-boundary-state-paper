"""BB-CBE-2I: polarized Gowdy T^3 + massless scalar, inhomogeneous equivalence.

Line E (BB-CBE). Executes exactly the frozen specification
Discussion/BB_CBE2I_INHOMOGENEOUS_EQUIVALENCE_SPEC_2026-08-04.md:
controls (manufactured Bessel solution, homogeneous limit, constraint
convergence) gate the bounded three-resolution campaign.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.special import j0, j1

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "CBE2I.json"

DT_FACTOR = 0.2
STORE_DTAU = 0.5
COMPARISON_TAUS = (2.0, 4.0, 8.0, 12.0)
CAMPAIGN_N = (256, 512, 1024)
DEPTHS = (16.0, 24.0)
D_TOL = 1.0e-3
D_MISMATCH = 1.0e-2
CONSTRAINT_GATE = 1.0e-3


def dtheta_of(n):
    return 2.0 * np.pi / n


def d1(f, dx):
    return (np.roll(f, -1) - np.roll(f, 1)) / (2.0 * dx)


def d2(f, dx):
    return (np.roll(f, -1) - 2.0 * f + np.roll(f, 1)) / dx**2


def periodic_antiderivative(f, dx):
    """Spectral antiderivative of a zero-mean periodic sample, anchored at 0."""
    n = len(f)
    k = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    F = np.fft.fft(f)
    F[0] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        A = np.where(k != 0.0, F / (1j * k), 0.0)
    out = np.real(np.fft.ifft(A))
    return out - out[0]


def rhs(tau, y, dx):
    P, Pt, phi, phit, lam = y
    e2 = np.exp(-2.0 * tau)
    return [Pt, e2 * d2(P, dx), phit, e2 * d2(phi, dx),
            Pt**2 + e2 * d1(P, dx)**2 + phit**2 + e2 * d1(phi, dx)**2]


def rk4(y, tau, dt, dx):
    def add(a, b, s):
        return [ai + s * bi for ai, bi in zip(a, b)]
    k1 = rhs(tau, y, dx)
    k2 = rhs(tau + dt / 2, add(y, k1, dt / 2), dx)
    k3 = rhs(tau + dt / 2, add(y, k2, dt / 2), dx)
    k4 = rhs(tau + dt, add(y, k3, dt), dx)
    return [yi + dt / 6 * (a + 2 * b + 2 * c + d)
            for yi, a, b, c, d in zip(y, k1, k2, k3, k4)]


def evolve(y0, tau0, tau1, n, store_taus=()):
    """Segment-wise integration landing EXACTLY on every stored surface."""
    dx = dtheta_of(n)
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
        steps = max(1, int(np.ceil(abs(seg) / (DT_FACTOR * dx))))
        dt = seg / steps
        for _ in range(steps):
            y = rk4(y, tau, dt, dx)
            tau += dt
        tau = t_next                      # suppress roundoff accumulation
        if t_next != tau1 or t_next in targets:
            stored[t_next] = [f.copy() for f in y]
    return y, tau, stored


def invariants(y, tau, dx):
    P, Pt, phi, phit, lam = y
    lam_tau = Pt**2 + np.exp(-2 * tau) * d1(P, dx)**2 \
        + phit**2 + np.exp(-2 * tau) * d1(phi, dx)**2
    N = np.exp((lam + tau) / 4.0 - tau)
    H_th = (lam_tau + 1.0) / (4.0 * N)
    H_x = (Pt - 1.0) / (2.0 * N)
    H_y = (-Pt - 1.0) / (2.0 * N)
    Theta = H_th + H_x + H_y
    H = Theta / 3.0
    sigma2 = (H_th - H)**2 + (H_x - H)**2 + (H_y - H)**2
    rho = 0.5 * np.exp(-(lam + tau) / 2.0) * (np.exp(2 * tau) * phit**2
                                              + d1(phi, dx)**2)
    lam_A = np.exp(-tau / 2.0)
    return {
        "Omega": rho / (3.0 * H**2),
        "Sigma2": sigma2 / (6.0 * H**2),
        "s_P": (H_x - H_y) / Theta,
        "clock_rate": Theta * lam_A,
        "lam_rho": rho**-0.5,
    }


def momentum_residual(y, tau, dx):
    P, Pt, phi, phit, lam = y
    target = 2.0 * (Pt * d1(P, dx) + phit * d1(phi, dx))
    got = d1(lam, dx)
    return float(np.max(np.abs(got - target))
                 / (np.max(np.abs(target)) + 1e-30))


def production_data(n):
    th = np.arange(n) * dtheta_of(n)
    P = 0.6 * np.cos(th)
    Pt = 0.4 + 0.3 * np.cos(2 * th)
    phi = 0.5 * np.sin(th)
    phit = 0.8 - 0.2 * np.cos(2 * th)
    integrand = 2.0 * (Pt * d1(P, dtheta_of(n)) + phit * d1(phi, dtheta_of(n)))
    lam = periodic_antiderivative(integrand, dtheta_of(n))
    return [P, Pt, phi, phit, lam]


# ---------------------------------------------------------------- controls --

def control_manufactured():
    rows = []
    for n in (128, 256, 512):
        th = np.arange(n) * dtheta_of(n)
        z0 = 3.0
        y0 = [j0(z0) * np.cos(3 * th), z0 * j1(z0) * np.cos(3 * th),
              np.zeros(n), np.zeros(n), np.zeros(n)]
        y, tau, _ = evolve(y0, 0.0, 8.0, n)
        exact = j0(3.0 * np.exp(-8.0)) * np.cos(3 * th)
        rows.append({"N": n, "sup_error": float(np.max(np.abs(y[0] - exact)))})
    orders = [float(np.log2(rows[i]["sup_error"] / rows[i + 1]["sup_error"]))
              for i in range(2)]
    passes = all(abs(o - 2.0) <= 0.3 for o in orders)
    return {"rows": rows, "orders": orders, "passes": bool(passes)}


def control_homogeneous():
    n = 64
    y0 = [np.zeros(n), 0.5 * np.ones(n), np.zeros(n), 0.7 * np.ones(n),
          np.zeros(n)]
    y, tau, stored = evolve(y0, 0.0, 12.0, n, store_taus=(4.0, 8.0, 12.0))
    dx = dtheta_of(n)
    vals = [invariants(stored[t], t, dx) for t in (4.0, 8.0, 12.0)]
    s2 = [float(v["Sigma2"][0]) for v in vals]
    om = [float(v["Omega"][0]) for v in vals]
    drift = max(max(s2) - min(s2), max(om) - min(om))
    # expansion branch: exact VTD data reproduces the same line trivially;
    # backward evolution must return the initial state to roundoff
    yb, _, _ = evolve(y, 12.0, 0.0, n)
    back = float(max(np.max(np.abs(yb[i] - y0[i])) for i in range(4)))
    passes = drift < 1e-10 and back < 1e-8
    return {"Sigma2_series": s2, "Omega_series": om,
            "invariant_drift": float(drift),
            "backward_roundtrip_error": back, "passes": bool(passes)}


def control_constraint():
    rows = []
    for n in (256, 512, 1024):
        y, tau, _ = evolve(production_data(n), 0.0, 8.0, n)
        rows.append({"N": n,
                     "residual": momentum_residual(y, 8.0, dtheta_of(n))})
    orders = [float(np.log2(rows[i]["residual"] / rows[i + 1]["residual"]))
              for i in range(2)]
    passes = (rows[-1]["residual"] < CONSTRAINT_GATE
              and all(o > 1.5 for o in orders))
    return {"rows": rows, "orders": orders, "passes": bool(passes)}


# ---------------------------------------------------------------- campaign --

def spike_scan(v, n):
    grad = np.max(np.abs(d1(v, dtheta_of(n))))
    return {"max_grad_v": float(grad), "flagged": bool(grad > n / 8.0)}


def run_campaign():
    out = {}
    for depth in DEPTHS:
        for n in CAMPAIGN_N:
            dx = dtheta_of(n)
            store = tuple(COMPARISON_TAUS) + (depth - 2.0, depth)
            y_end, _, col_store = evolve(production_data(n), 0.0, depth, n,
                                         store_taus=store)
            col_store[depth] = [f.copy() for f in y_end]
            vP = (col_store[depth][0] - col_store[depth - 2.0][0]) / 2.0
            vF = (col_store[depth][2] - col_store[depth - 2.0][2]) / 2.0
            Pinf = col_store[depth][0] - vP * depth
            Finf = col_store[depth][2] - vF * depth
            spikes = spike_scan(vP, n)
            e2 = np.exp(-2.0 * depth)
            lam_t = vP**2 + e2 * d1(vP * depth + Pinf, dx)**2 \
                + vF**2 + e2 * d1(vF * depth + Finf, dx)**2
            lam_theta = 2.0 * (vP * d1(vP * depth + Pinf, dx)
                               + vF * d1(vF * depth + Finf, dx))
            lam0 = periodic_antiderivative(lam_theta, dx)
            lam0 += col_store[depth][4][0]      # anchor at theta = 0
            y_exp = [vP * depth + Pinf, vP.copy(), vF * depth + Finf,
                     vF.copy(), lam0]
            _, _, exp_store = evolve(y_exp, depth, 0.0, n,
                                     store_taus=COMPARISON_TAUS)
            per_surface = {}
            for tc in COMPARISON_TAUS:
                a = invariants(col_store[tc], tc, dx)
                b = invariants(exp_store[tc], tc, dx)
                comps = [a[k] - b[k] for k in
                         ("Omega", "Sigma2", "s_P", "clock_rate")]
                D = float(np.max(np.sqrt(sum(c**2 for c in comps))))
                clock_gap = float(np.max(np.abs(a["lam_rho"] - b["lam_rho"])
                                         / a["lam_rho"]))
                per_surface[str(tc)] = {"D": D, "lam_rho_relgap": clock_gap}
            out[f"depth{int(depth)}_N{n}"] = {
                "surfaces": per_surface, "spikes": spikes,
                "momentum_residual_collapse_end":
                    momentum_residual(col_store[depth], depth, dx),
            }
    return out


def classify(campaign):
    fine = campaign["depth24_N1024"]["surfaces"]
    d_fine = max(v["D"] for v in fine.values())
    clock_ok = all(v["lam_rho_relgap"] < D_TOL * 10 for v in fine.values())
    spikes = any(campaign[k]["spikes"]["flagged"] for k in campaign)
    dec_N = all(
        max(v["D"] for v in campaign[f"depth24_N{n2}"]["surfaces"].values())
        <= max(v["D"] for v in campaign[f"depth24_N{n1}"]["surfaces"].values())
        for n1, n2 in ((256, 512), (512, 1024))
    )
    d16 = max(v["D"] for v in campaign["depth16_N1024"]["surfaces"].values())
    depth_ok = d_fine <= d16 or d_fine < D_TOL
    if d_fine <= D_TOL and dec_N and depth_ok and clock_ok and not spikes:
        return ("CBE-B extended to the polarized-inhomogeneous sector "
                "(CBE-2I-EQ): equivalence survives gradients; CBE-A not "
                "claimable, unpolarized/spiky sector untested")
    if d_fine > D_MISMATCH and dec_N:
        return "CBE-C — inhomogeneous invariant mismatch"
    if d_fine > D_TOL:
        return "intermediate — numbers reported, no classification"
    return "CBE-D — numerical nonclosure"


def main():
    controls = {
        "manufactured": control_manufactured(),
        "homogeneous": control_homogeneous(),
        "constraint": control_constraint(),
    }
    gates = all(c["passes"] for c in controls.values())
    result = {"program": "BB-CBE-2I", "line": "E",
              "spec": "BB_CBE2I_INHOMOGENEOUS_EQUIVALENCE_SPEC_2026-08-04.md",
              "controls": controls, "controls_pass": bool(gates)}
    if gates:
        campaign = run_campaign()
        result["campaign"] = campaign
        result["classification"] = classify(campaign)
    else:
        result["classification"] = ("controls failed — implementation "
                                    "failure; campaign not started")
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    brief = {"controls_pass": gates,
             "classification": result["classification"]}
    if gates:
        brief["D_fine_by_surface"] = {
            k: v["D"] for k, v in
            result["campaign"]["depth24_N1024"]["surfaces"].items()}
    print(json.dumps(brief, indent=2))


if __name__ == "__main__":
    main()
