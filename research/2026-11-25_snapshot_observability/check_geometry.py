"""Reconstruct with the Einstein-consistent metric, using frozen field solvers."""
import importlib.util
import json
from pathlib import Path
import platform

import numpy as np

PAPER = Path(__file__).resolve().parents[2]
KEYS = ("Omega", "Sigma2", "s_P", "clock_rate")


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    item = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(item)
    return item


def observables(state, tau, dx, d1):
    P, Pt, phi, phit, lam = state
    e2 = np.exp(-2.0 * tau)
    scalar_energy = phit**2 + e2 * d1(phi, dx)**2
    lam_tau = Pt**2 + e2 * d1(P, dx)**2 + scalar_energy
    lapse = np.exp((-lam - 3.0 * tau) / 4.0)
    h_theta = (1.0 - lam_tau) / (4.0 * lapse)
    h_x = (Pt - 1.0) / (2.0 * lapse)
    h_y = (-Pt - 1.0) / (2.0 * lapse)
    theta = h_theta + h_x + h_y
    H = theta / 3.0
    rho = scalar_energy / (4.0 * lapse**2)
    shear_square = sum((rate - H)**2 for rate in (h_theta, h_x, h_y))
    return {
        "Omega": rho / (3.0 * H**2),
        "Sigma2": shear_square / (6.0 * H**2),
        "s_P": (h_x - h_y) / theta,
        "clock_rate": theta * np.exp(-tau / 2.0),
        "lam_rho": np.divide(1.0, np.sqrt(rho), out=np.full_like(rho, np.inf), where=rho>0),
        "rho": rho,
        "H_theta": h_theta, "H_x": h_x, "H_y": h_y,
    }


def homogeneous_check(model):
    cases = {}
    for p, v in ((1.0, 0.0), (0.0, float(np.sqrt(3))), (0.37, 0.85)):
        tau, n = 2.0, 16
        hom = [np.full(n, a) for a in (p*tau, p, v*tau, v, (p*p+v*v)*tau)]
        o = model.invariants(hom, tau, model.dtheta_of(n))
        h0, h1, h2 = (o[key] for key in ("H_theta","H_x","H_y"))
        err = float(np.max(np.abs((h0*h1 + h0*h2 + h1*h2 - o["rho"])
                                  / (h0*h0 + h1*h1 + h2*h2))))
        assert err < 2e-15
        assert np.max(np.abs(o["Omega"] + o["Sigma2"] - 1)) < 2e-15
        if p == 0:
            assert np.max(np.abs(o["Sigma2"])) < 2e-15
        cases[str((p,v))] = {"relative_ADM_error": err, "Omega": float(o["Omega"][0]),
                              "Sigma2": float(o["Sigma2"][0])}
    return cases


def compare(m, ref, got, tau, dx):
    actual = m.invariants(got, tau, dx)
    truth = m.invariants(ref, tau, dx)
    diffs = [actual[k] - truth[k] for k in KEYS]
    return {
        "D": float(np.max(np.sqrt(sum(delta**2 for delta in diffs)))),
        "shape_D": float(np.max(np.sqrt(sum(delta**2 for delta in diffs[:3])))),
        "density_relative_error": float(np.max(np.abs(actual["lam_rho"] / truth["lam_rho"] - 1.0))),
    }


def holdout_restart(m, depth, dx, s0, s1, metric=None):
    """Only pre-endpoint data are passed to the main estimator."""
    t1 = depth - 2.0
    vp, vf = ((s1[i] - s0[i]) / 2.0 for i in (0, 2))
    p_inf, f_inf = s1[0] - t1*vp, s1[2] - t1*vf
    q = vp**2 + vf**2
    primitive = m.periodic_antiderivative(2.0 * (vp*m.d1(p_inf,dx) + vf*m.d1(f_inf,dx)), dx)
    c = float(s1[4][0] - t1*q[0]) if metric is None else metric
    return [vp*depth+p_inf, vp, vf*depth+f_inf, vf, q*depth+primitive+c], c


def production(m):
    out = {}
    for depth in (8.0, 12.0, 16.0, 24.0):
        for n in (256, 512, 1024):
            surfaces = (2.0, 4.0)
            dx = m.dtheta_of(n)
            end, _, snapshots = m.evolve(m.production_data(n), 0.0, depth, n,
                                          store_taus=tuple(sorted(set(surfaces+(depth-4, depth-2)))))
            restart, c = holdout_restart(m, depth, dx, snapshots[depth-4], snapshots[depth-2])
            arms = {"pre_endpoint_metric": (restart, c)}
            # Counterfactuals change one additive constant and leave all reduced
            # derivatives unchanged, so one integration supplies all three arms.
            endpoint_c = c + float(end[4][0] - restart[4][0])
            arms["imported_endpoint_metric"] = (restart, endpoint_c)
            arms["zero_asymptotic_constant"] = (restart, 0.0)
            _, _, recovered = m.evolve(restart, depth, 0.0, n, store_taus=surfaces)
            scores = {}
            for name, (_, adopted_c) in arms.items():
                rows = {}
                for tau in surfaces:
                    shifted = recovered[tau][:4] + [recovered[tau][4] + adopted_c - c]
                    rows[str(tau)] = compare(m, snapshots[tau], shifted, tau, dx)
                scores[name] = {"asymptotic_constant": adopted_c, "surfaces": rows,
                                "D_max": max(x["D"] for x in rows.values())}
            out[f"{int(depth)}_{n}"] = scores
    return out


def bessel(m, fg):
    out = {}
    for depth in (8.0, 12.0, 16.0):
        for n in (256, 512, 1024):
            dx = m.dtheta_of(n)
            truth = fg.exact_state(m, 2.0, n)
            exact_input, _, _, _ = fg.construct_arm(m, depth, n, "FUCHSIAN")
            s0, s1 = (fg.exact_state(m, depth-offset, n) for offset in (4.0,2.0))
            holdout, _ = holdout_restart(m, depth, dx, s0, s1)
            arms = {}
            for name, restart in (("analytic_coefficients", exact_input), ("pre_endpoint_holdout", holdout)):
                recovered, _, _ = m.evolve(restart, depth, 2.0, n)
                arms[name] = {
                    **compare(m, truth, recovered, 2.0, dx),
                    "P_error": float(np.max(np.abs(truth[0]-recovered[0]))),
                    "phi_error": float(np.max(np.abs(truth[2]-recovered[2]))),
                    "lambda_error": float(np.max(np.abs(truth[4]-recovered[4]))),
                }
            out[f"{int(depth)}_{n}"] = arms
    return out


if __name__ == "__main__":
    model = load(PAPER / "reproduction" / "inputs" / "bb-cbe2i" / "run_cbe2i.py", "e1_smooth_source")
    fg = load(PAPER / "verification" / "run_e1_fg1.py", "e1_exact_source")
    model.invariants = lambda state, tau, dx: observables(state, tau, dx, model.d1)
    print(json.dumps({"platform": platform.platform(), "homogeneous": homogeneous_check(model),
                      "production": production(model), "bessel": bessel(model, fg)},
                     indent=2, allow_nan=False))
