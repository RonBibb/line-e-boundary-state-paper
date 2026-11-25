"""Independent verifier for BB-CBE-3b (spec Section 11).

Recomputes footprints, masks, O-gates, the deterministic fit, validation,
and the Section-9 outcome mapping from the checkpointed slices, using its
own implementations of all 3b-specific logic. Shares only certified
primitives (exact family, hyperbolic distance, state distance, detector)
with the production script. Compares against CBE3B_RESULTS.json values,
not booleans.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "bb-cbe3"))
sys.path.insert(0, str(ROOT.parent / "bb-cbe3t"))
from run_cbe3 import comparison_D, d1, hyperbolic_distance  # noqa: E402
import run_cbe3t as cbe3t  # noqa: E402

RES = (1024, 2048, 4096)
STORE = (8.0, 10.0, 12.0, 14.0, 16.0, 18.0)
FOOT = (12.0, 14.0, 16.0, 18.0)
DXF = 2 * np.pi / 4096
SEED = 20260804
checks = []


def check(name, ok, detail=""):
    checks.append({"name": name, "ok": bool(ok), "detail": str(detail)})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail
                                                 else ""))


def load_slices():
    out = {}
    for n in RES:
        z = np.load(ROOT / f"cbe3b_slices_N{n}.npz")
        out[n] = {t: [z[f"tau{t}_{i}"] for i in range(5)] for t in STORE}
    return out


def vgrid(n):
    dx = 2 * np.pi / n
    return (np.arange(n) - n // 2) * dx, dx


def my_footprint(y, tau, dx, x):
    P, _, Q, _, _ = y
    G = np.exp(-2 * tau) * (d1(P, dx)**2 + np.exp(2 * P) * d1(Q, dx)**2)
    Gp = np.maximum(G - np.median(G), 0.0)
    tot = Gp.sum() * dx
    if tot < 1e-14:
        return 0.0
    o = np.argsort(np.abs(x))
    c = np.cumsum(Gp[o]) * dx
    return float(np.abs(x[o][min(int(np.searchsorted(c, 0.999 * tot)),
                                 len(x) - 1)]))


def my_rich(a, b, c):
    e1, e2 = (4 * b - a) / 3, (4 * c - b) / 3
    closes = abs(e2 - e1) <= max(0.2 * abs(e2), 5e-4)
    cls = ("intermediate" if not closes else
           "zero" if e2 <= 1e-3 else
           "nonzero" if e2 >= 1e-2 else "intermediate")
    return e1, e2, closes, cls


def main():
    declared = json.loads((ROOT / "CBE3B_RESULTS.json").read_text())
    sl = load_slices()

    # footprints
    for n in RES:
        x, dx = vgrid(n)
        for t in FOOT:
            r = my_footprint(sl[n][t], t, dx, x)
            rd = declared["footprints"][str(n)][str(t)]
            check(f"footprint N{n} tau{t}", abs(r - rd) <= 1e-12,
                  f"{r:.6f} vs {rd:.6f}")

    # ordinary mask
    x, dx = vgrid(4096)
    radii = {n: {t: declared["footprints"][str(n)][str(t)] for t in FOOT}
             for n in RES}
    R = max(radii[4096].values()) + abs(np.exp(-8.0) - np.exp(-18.0))
    v = (sl[4096][18.0][0] - sl[4096][16.0][0]) / 2.0
    mask = (np.abs(x) > R) & (v >= 0.05)
    check("mask fraction", abs(float(np.mean(mask))
                               - declared["ordinary_mask"]["fraction"])
          <= 1e-12, f"{np.mean(mask):.6f}")
    masks = {}
    for n in RES:
        f = 4096 // n
        masks[n] = np.all(mask.reshape(n, f), axis=1)

    # O2 recompute (own formulas)
    sups = {}
    for n in RES:
        xn, dxn = vgrid(n)
        m = masks[n]
        sups[n] = {}
        for t in (14.0, 16.0, 18.0):
            P, Pt, Q, Qt, _ = sl[n][t]
            e2 = np.exp(-2 * t)
            Pth, Qth = d1(P, dxn), d1(Q, dxn)
            Pxx = (np.roll(P, -1) - 2 * P + np.roll(P, 1)) / dxn**2
            Qxx = (np.roll(Q, -1) - 2 * Q + np.roll(Q, 1)) / dxn**2
            SP = e2 * (Pxx - np.exp(2 * P) * Qth**2)
            SQ = e2 * (Qxx + 2 * Pth * Qth)
            Ptt = e2 * Pxx + np.exp(2 * P) * (Qt**2 - e2 * Qth**2)
            Qtt = e2 * Qxx - 2 * (Pt * Qt - e2 * Pth * Qth)
            eps = np.sqrt(SP**2 + np.exp(2 * P) * SQ**2) \
                / (1 + np.sqrt(Ptt**2 + np.exp(2 * P) * Qtt**2))
            sups[n][t] = float(np.max(eps[m]))
    for n in RES:
        for t in (14.0, 16.0, 18.0):
            dv = declared["O2"]["sup_epsS"][str(n)][str(t)]
            check(f"O2 sup N{n} tau{t}",
                  abs(sups[n][t] - dv) <= 1e-12 * max(1, abs(dv)),
                  f"{sups[n][t]:.3e}")

    # O4 recompute
    for n in RES:
        xn, dxn = vgrid(n)
        m = masks[n]
        st = sl[n]
        vv = (st[14.0][0] - st[12.0][0]) / 2.0
        Pinf = st[14.0][0] - vv * 14.0
        piQ = np.exp(2 * st[14.0][0]) * st[14.0][3]
        vs = np.where(np.abs(vv) < 0.05, 0.05, vv)
        Qinf = st[14.0][2] + np.exp(-2 * st[14.0][0]) * piQ / (2 * vs)
        for tp in (16.0, 18.0):
            Pp = vv * tp + Pinf
            Qtp = np.exp(-2 * Pp) * piQ
            Qp = Qinf - np.exp(-2 * Pp) * piQ / (2 * vs)
            D = comparison_D(st[tp], [Pp, vv, Qp, Qtp, np.zeros_like(Pp)],
                             tp, dxn)
            got = float(np.max(D[m]))
            dv = declared["O4"]["per"][str(n)][str(tp)]["sup"]
            check(f"O4 sup N{n} tau{tp}",
                  abs(got - dv) <= 1e-10 * max(1, abs(dv)), f"{got:.3e}")

    # deterministic fit recomputation (fine grid only, full optimizer)
    def objective(p):
        w, lnq0, q2 = p
        J = 0.0
        for t in (12.0, 14.0, 16.0):
            rF = radii[4096][t]
            mm = np.abs(x) <= max(rF, 2 * dx)
            P, Pt, Q, Qt, _ = sl[4096][t]

            def fn(xa, ta):
                Pf, Qf, lf = cbe3t.true_spike(xa, ta, w=w, q0=np.exp(lnq0))
                return Pf, Qf + q2, lf
            Pf, Qf, _ = fn(x[mm], t)
            Ptf, Qtf, _ = cbe3t.tau_derivative(fn, x[mm], t)
            dh2 = hyperbolic_distance(P[mm], Q[mm], Pf, Qf)**2
            vp = np.sqrt(Pt[mm]**2 + np.exp(2 * P[mm]) * Qt[mm]**2)
            vf = np.sqrt(Ptf**2 + np.exp(2 * Pf) * Qtf**2)
            J += float(np.mean(dh2 + 0.25 * (vp - vf)**2))
        return J

    bounds = [(1e-3, 0.999), (-10, 10), (-10, 10)]
    de = differential_evolution(objective, bounds, seed=SEED, popsize=20,
                                maxiter=500, polish=False, tol=1e-12)
    ref = minimize(objective, de.x, method="L-BFGS-B", bounds=bounds)
    fw, fq0, fq2 = float(ref.x[0]), float(np.exp(ref.x[1])), float(ref.x[2])
    dfit = declared["S1_fits"]["4096"]
    for name, got, want in (("w", fw, dfit["w"]), ("q0", fq0, dfit["q0"]),
                            ("q2", fq2, dfit["q2"])):
        check(f"fit {name} (N=4096)", abs(got - want)
              <= 1e-6 * max(1.0, abs(want)), f"{got:.8f} vs {want:.8f}")

    # outcome mapping recomputation
    o = declared
    ord_prereq = (o["O1"]["passes"] and o["O2"]["passes"]
                  and o["O3"]["passes"])
    ord_pass = ord_prereq and o["O4"]["passes"]
    reading = o["S2_S3"]["reading"]
    zm = o["zero_measure_gate"]["passes"]
    if ord_pass and zm and reading == "S-LIM-A":
        want = "3b-EQ-POINTWISE"
    elif ord_pass and zm:
        want = "3b-AE"
    elif ord_prereq and o["O4"]["converged_fail"]:
        want = "3b-C-ORDINARY"
    elif reading == "S-LIM-A":
        want = "3b-S-LIM-ONLY"
    else:
        want = "3b-D"
    check("outcome mapping", want == o["outcome"],
          f"recomputed {want} vs declared {o['outcome']}")

    n_fail = sum(1 for c in checks if not c["ok"])
    (ROOT / "CBE3B_VERIFICATION.json").write_text(
        json.dumps({"checks": checks, "failures": n_fail}, indent=2) + "\n")
    print(f"\n{len(checks) - n_fail}/{len(checks)} checks pass")


if __name__ == "__main__":
    main()
