"""Reproduce two-mode rank-two collisions after the Gowdy-to-Ernst lift."""
import json
import math
import platform

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import j0, j1, y0, y1


def annulus_roots(ratio):
    def determinant(z):
        return j0(z/ratio)*y0(z) - y0(z/ratio)*j0(z)

    grid = np.linspace(0.001, 3.5*np.pi/(1-1/ratio), 1536)
    values = determinant(grid)
    brackets = [(a, b) for a, b, fa, fb in zip(grid, grid[1:], values, values[1:])
                if fa*fb < 0]
    return [brentq(determinant, a, b, xtol=2e-14) for a, b in brackets[:2]]


def eigenvalue_ratio(ratio):
    first, second = annulus_roots(ratio)
    return second/first


def collision():
    log_ratio = brentq(lambda trial: eigenvalue_ratio(math.exp(trial))-13/6,
                       math.log(1.1), 100, xtol=1e-12)
    ratio = math.exp(log_ratio)
    first, second = annulus_roots(ratio)
    b = first/6
    a = b/ratio
    g, q, P0 = 0.55, 0.37, 0.2
    sigma = math.hypot(g, q)
    angular = np.linspace(0, 2*np.pi, 2048, endpoint=False)
    modes = [(6, 0., 0.02), (13, np.pi*np.sqrt(2), 0.014)]
    constant_current = q*math.exp(P0)

    def phase(t):
        w = np.zeros_like(angular)
        wt, wx, ernst_q = w.copy(), w.copy(), w.copy()
        for k, shift, amplitude in modes:
            z = k*a
            profile = np.pi/2*(j0(z)*y0(k*t)-y0(z)*j0(k*t))
            derivative = np.pi/2*k*(-j0(z)*y1(k*t)+y0(z)*j1(k*t))
            sine, cosine = np.sin(k*angular+shift), np.cos(k*angular+shift)
            w += amplitude*profile*sine
            wt += amplitude*derivative*sine
            wx += amplitude*k*profile*cosine
            ernst_q -= constant_current*t*amplitude*derivative*cosine/k
        return w, wt, wx, ernst_q

    def ernst(t, g):
        w, wt, wx, Q = phase(t)
        u = np.cosh(sigma*w)+g*np.sinh(sigma*w)/sigma
        P = P0+np.log(u)
        P_w = (sigma*np.sinh(sigma*w)+g*np.cosh(sigma*w))/u
        Pt = -P_w*wt-1/t
        Ptheta = -P_w*wx
        Qt = constant_current*t*wx
        Qtheta = constant_current*t*wt
        return -P-math.log(t), Q, Pt, Qt, Pt*Qtheta-Ptheta*Qt, Ptheta, Qtheta

    interior = math.sqrt(a*b)
    left, right = ([ernst(t, sign*g) for t in (a, interior, b)] for sign in (1., -1.))
    absolute = lambda vector: float(np.max(np.abs(vector)))
    record_gap = max(absolute(left[i][field]-right[i][field])
                     for i in (0, 2) for field in (0, 1, 3))
    assert record_gap < 1e-10
    initial_gap = absolute(left[0][2]-right[0][2])
    rank_witness = min(absolute(left[0][4]), absolute(right[0][4]))
    assert initial_gap > 0.01 and rank_witness > 0.001
    def lambda_log_derivative(log_t, branch):
        t = math.exp(log_t)
        P, _, Pt, Qt, _, Ptheta, Qtheta = ernst(t, branch*g)
        energy = Pt**2+Ptheta**2+np.exp(2*P)*(Qt**2+Qtheta**2)
        return float(-t*t*np.mean(energy))

    def phase_energy_log_integrand(log_t):
        t = math.exp(log_t)
        _, wt, wtheta, _ = phase(t)
        return float(t*t*constant_current**2*np.mean(wt**2+wtheta**2))

    integration_window = (math.log(a), math.log(b))
    I = quad(phase_energy_log_integrand, *integration_window, epsabs=1e-12)[0]
    shifts = [quad(lambda log_t: lambda_log_derivative(log_t, branch),
                   *integration_window, epsabs=1e-12)[0] for branch in (1., -1.)]
    inferred_kappa = math.sqrt((-shifts[0]-math.log(ratio))/I)
    psi = constant_current*phase(interior)[0]
    index = int(np.argmax(np.abs(psi)))
    inferred = [(math.exp(-sol[1][0][index]-math.log(interior)-P0)
                 - math.cosh(inferred_kappa*psi[index]))
                *inferred_kappa/math.sinh(inferred_kappa*psi[index]) for sol in (left, right)]
    assert max(abs(estimate-sign*g/constant_current)
               for estimate, sign in zip(inferred,(1., -1.))) < 1e-9
    return {"platform": platform.platform(), "ratio": ratio, "a": a, "b": b,
            "outer_bessel_roots": [first, second],
            "inner_bessel_frequencies": [6*a, 13*a],
            "initial_ernst_Q_peak_to_peak": float(np.ptp(left[0][1])),
            "both_snapshots_P_Q_Qt_gap": record_gap,
            "both_snapshots_lambda_gap_from_transform": 2*max(
                absolute(left[i][0]-right[i][0]) for i in (0, 2)),
            "initial_Pt_gap": initial_gap,
            "smaller_maximum_absolute_rank_determinant_at_a": rank_witness,

            "interior_P_gap": absolute(left[1][0]-right[1][0]),
            "terminal_metric_shifts_from_quadrature": shifts,
            "expected_branch_parameters": [g/constant_current, -g/constant_current],
            "recovered_branch_parameters": inferred}


if __name__ == "__main__":
    print(json.dumps(collision(), indent=2, allow_nan=False))
