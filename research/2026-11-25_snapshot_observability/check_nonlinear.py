"""Use identical finite-depth data to test two vacuum Gowdy restart rules."""
import importlib.util
import json
from pathlib import Path
import platform

import numpy as np
from scipy.optimize import brentq

INPUTS = Path(__file__).resolve().parents[2] / "reproduction" / "inputs"
TIMES = (2.0, 4.0, 8.0)
WINDOW = 2.0
MAX_ITERATIONS = 12
FIT_TOLERANCE = 1e-10


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, INPUTS / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def use_fourth_order_rhs(model):
    def d1(f, h):
        return (np.roll(f, 2) - 8.0*np.roll(f,1) + 8.0*np.roll(f,-1) - np.roll(f,-2))/(12.0*h)

    def d2(f, h):
        return (-np.roll(f,2) + 16.0*np.roll(f,1) - 30.0*f
                + 16.0*np.roll(f,-1) - np.roll(f,-2))/(12.0*h*h)

    def rhs(tau, state, dx):
        P, V, Q, R, _ = state
        e2, eP2 = np.exp(-2.0*tau), np.exp(2.0*P)
        P_x, Q_x = d1(P, dx), d1(Q, dx)
        return [V, e2*d2(P,dx) + eP2*(R*R-e2*Q_x*Q_x),
                R, e2*d2(Q,dx) - 2.0*(V*R-e2*P_x*Q_x),
                V*V + e2*P_x*P_x + eP2*(R*R+e2*Q_x*Q_x)]

    model.rhs = rhs


def infer_velocity(evolve, depth, observed):
    P_old, P, Q, R = observed
    zero_metric = np.zeros_like(P)
    secant = (P-P_old)/WINDOW
    V = secant.copy()
    mismatches = []
    for _ in range(MAX_ITERATIONS + 1):
        earlier, _, _ = evolve([P,V,Q,R,zero_metric], depth, depth-WINDOW, len(P))
        # The observed earlier point minus the trial earlier point equals
        # minus WINDOW times the missing velocity in the zero-force limit.
        defect = (earlier[0]-P_old)/WINDOW
        residual = float(np.max(np.abs(defect)))
        mismatches.append(residual)
        if residual < FIT_TOLERANCE:
            break
        if len(mismatches) <= MAX_ITERATIONS:
            V = V + defect
    return secant, V, mismatches


def score(model, dx, trial, references, mask):
    rows = {}
    for tau, got in trial.items():
        truth = references[tau]
        D = model.comparison_D(truth, got, tau, dx)
        rows[str(tau)] = {
            "all": float(np.max(D)),
            "detector_mask": float(np.max(D[mask])) if np.any(mask) else None,
            "complement": float(np.max(D[~mask])) if np.any(~mask) else None,
        }
    return rows


def benchmark(scoring, model, end, earlier, references, depth, mask):
    P, V_truth, Q, R, _ = end
    observed = (earlier[0], P, Q, R)
    secant, fitted, mismatches = infer_velocity(model.evolve, depth, observed)
    n, dx = len(P), model.dtheta_of(len(P))
    zero_metric = np.zeros_like(P)
    arms = {}
    for name, velocity in (("secant",secant),("matched_window",fitted),("privileged",V_truth)):
        # Earlier reference solutions and V_truth are available only here.
        _, _, states = model.evolve([P,velocity,Q,R,zero_metric],
                                     depth, min(references), n, store_taus=tuple(references))
        arms[name] = {"velocity_error": float(np.max(np.abs(velocity-V_truth))),
                      "scores": score(scoring, dx, states, references, mask)}
    return {"fit_residuals": mismatches, "converged": mismatches[-1] < FIT_TOLERANCE, "arms": arms}


def centered_template_fit(x, t, P, v, Q):
    # exp(-P) is quadratic in x on each Lim slice. Calculate the center from
    # the sampled field alone, then evaluate both time samples there.
    center_guess = float(x[np.argmax(P)])
    dx = float(x[1]-x[0])
    selected = np.abs(x-center_guess) <= 3.25*dx
    coordinate = (x[selected]-center_guess)/dx
    def quadratic(samples):
        return np.polynomial.polynomial.polyfit(coordinate, samples[selected], 2)
    coeff = quadratic(np.exp(-P))
    center = -coeff[1]/(2.0*coeff[2])
    evaluate = lambda c: np.polynomial.polynomial.polyval(center,c)
    p_end = -np.log(evaluate(coeff))
    p_old = -np.log(evaluate(quadratic(np.exp(-P+WINDOW*v))))
    v_center = (p_end-p_old)/WINDOW
    q_center = evaluate(quadratic(Q))
    def defect(w):
        return 2.0 - (np.log(np.cosh(w*t))-np.log(np.cosh(w*(t-WINDOW))))/WINDOW - v_center
    w = brentq(defect, 0.001, 0.999)
    q0 = 0.5*np.exp(2.0*t-np.log(np.cosh(w*t))-p_end)
    q2 = q_center + q0*w*np.exp(-2.0*t)
    return (w,q0,q2,center_guess+dx*center)


def blind_local_control(model, exact_spike):
    rows = {}
    t, target = 6.0, 2.0
    for n in (1024, 2048, 4096):
        x = np.arange(-n//2,n//2)*model.dtheta_of(n)
        for w,q0,q2,direction in ((0.35,0.30,0.17,1.0),(0.85,0.80,-0.11,-1.0)):
            dx = x[1]-x[0]
            true_center = (int(0.4/dx) + 0.5)*dx*direction
            P,Q,_ = exact_spike.true_spike(x-true_center,t,w,q0)
            P_old,_,_ = exact_spike.true_spike(x-true_center,t-WINDOW,w,q0)
            w_hat,q0_hat,q2_hat,x_hat = centered_template_fit(x,t,P,(P-P_old)/WINDOW,Q+q2)
            truth_P,truth_Q,_ = exact_spike.true_spike(x-true_center,target,w,q0)
            got_P,got_Q,_ = exact_spike.true_spike(x-x_hat,target,w_hat,q0_hat)
            mask = np.abs(x-x_hat) <= 0.8
            dH = model.hyperbolic_distance(truth_P[mask],truth_Q[mask]+q2,
                                          got_P[mask],got_Q[mask]+q2_hat)
            rows[f"{w}_{n}"] = {"center_error": abs(true_center-x_hat),
                                "parameter_errors": [abs(w-w_hat),abs(q0-q0_hat),abs(q2-q2_hat)],
                                "target_distance": float(np.max(dH))}
    return rows


if __name__ == "__main__":
    legacy = load("bb-cbe3/run_cbe3.py", "original")
    independent_stencil = load("bb-cbe3/run_cbe3.py", "fourth_order")
    use_fourth_order_rhs(independent_stencil)
    # Used only for a known-family local algebraic control.
    exact_spike = load("bb-cbe3t/run_cbe3t.py", "lim_source")

    output = {"platform": platform.platform(), "production": {},
              "exact_periodic": {}, "blind_local_lim": blind_local_control(legacy, exact_spike)}
    for n in (256,512,1024):
        depth = 8.0
        exact = lambda tau: legacy.exact_transformed_state(n,tau,0.40,0.80,3,1.0)[0]
        output["exact_periodic"][str(n)] = benchmark(
            legacy, independent_stencil, exact(depth), exact(depth-WINDOW),
            {tau: exact(tau) for tau in (2.0,4.0)}, depth, np.zeros(n,dtype=bool))

    for n in (1024,2048,4096):
        # Keep the integration segment at tau=12 out of the tau=18 lineage.
        # This follows the corresponding original campaign setup exactly.
        for depth in (12.0,18.0):
            earlier = depth-WINDOW
            end,_,snapshots = legacy.evolve(legacy.production_data(n),0.0,depth,n,store_taus=TIMES+(earlier,))
            dx = legacy.dtheta_of(n)
            v = (end[0]-snapshots[earlier][0])/WINDOW
            components = legacy.spike_components(v,end[2],dx)
            mask = legacy.region_masks(n,components,v)["spike"]
            refs = {tau:snapshots[tau] for tau in TIMES}
            output["production"][f"{int(depth)}_{n}"] = {
                "masked_cells": int(np.sum(mask)),
                "same_stencil": benchmark(legacy, legacy, end, snapshots[earlier], refs, depth, mask),
                "fourth_order_restart": benchmark(legacy, independent_stencil, end, snapshots[earlier], refs, depth, mask),
            }
    print(json.dumps(output, indent=2, allow_nan=False))
