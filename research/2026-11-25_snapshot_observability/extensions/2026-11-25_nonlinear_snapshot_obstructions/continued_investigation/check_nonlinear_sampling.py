"""Delayed displacement recovery on and near a rank-two Gowdy collision."""
from dataclasses import dataclass
import json
import math
import platform

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import j0, j1, y0, y1


@dataclass(frozen=True)
class Collision:
    a: float = 1.6
    mode: int = 24
    amplitude: float = 0.03
    P0: float = 0.2
    g: float = 0.55
    q: float = 0.37

    @property
    def current(self):
        return math.exp(self.P0)*self.q

    def radial_phase(self, t):
        k = self.mode
        start = k*self.a
        value = math.pi/2*(j0(start)*y0(k*t)-y0(start)*j0(k*t))
        deriv = math.pi*k/2*(-j0(start)*y1(k*t)+y0(start)*j1(k*t))
        return value, deriv

    @property
    def gap(self):
        k = self.mode
        return brentq(lambda d: self.radial_phase(self.a+d)[0],
                      0.8*math.pi/k, 1.2*math.pi/k, xtol=1e-15)

    def exact_fields(self, t, theta, branch):
        u, ut = self.radial_phase(t)
        w = self.amplitude*u*np.cos(self.mode*theta)
        wt = self.amplitude*ut*np.cos(self.mode*theta)
        g = branch*self.g
        speed = math.hypot(g, self.q)
        factor = np.cosh(speed*w)+g*np.sinh(speed*w)/speed
        P = -self.P0-np.log(factor)-math.log(t)
        P_w = (speed*np.sinh(speed*w)+g*np.cosh(speed*w))/factor
        Q = self.current*t*self.amplitude*ut*np.sin(self.mode*theta)/self.mode
        Qt = -self.current*t*self.amplitude*self.mode*u*np.sin(self.mode*theta)
        return np.stack([P, -P_w*wt-1/t, Q, Qt])


@dataclass(frozen=True)
class Schedule:
    a: float
    delays: np.ndarray
    mode_slots: np.ndarray
    weights: np.ndarray
    gains: np.ndarray

    @property
    def K(self):
        return len(self.mode_slots)-1


def schedule(a, window, K, scheme="delayed"):
    modes = np.arange(K+1)
    base_delay = window/4
    limit = math.ceil(2*math.pi/window)
    mode_slots = np.zeros(K+1, dtype=int)
    times = [base_delay]
    if scheme == "delayed":
        low = limit
        while low <= K:
            indices = np.arange(low, min(2*low, K+1))
            added = base_delay+math.pi/(4*low)
            use_added = abs(np.sin(indices*added)) >= abs(np.sin(indices*base_delay))
            if np.any(use_added):
                mode_slots[indices[use_added]] = len(times)
                times.append(added)
            low *= 2
    return retime(a, np.array(times), mode_slots)


def retime(a, delays, mode_slots):
    modes = np.arange(len(mode_slots))
    per_mode = delays[mode_slots]
    weights = np.maximum(modes, 1/per_mode)
    gains = weights*per_mode*np.sinc(modes*per_mode/math.pi)
    return Schedule(a, delays, mode_slots, weights, gains)


def theta_grid(n):
    return 2*math.pi*np.arange(n)/n


def series(values):
    return np.fft.rfft(values, axis=-1)/values.shape[-1]


def to_grid(coefficients, n):
    return np.fft.irfft(n*coefficients, n=n)


def h1_norm(coefficients):
    modes = np.arange(len(coefficients))
    return float(np.sqrt(abs(coefficients[0])**2+
                         2*np.sum((1+modes[1:]**2)*abs(coefficients[1:])**2)))


def spectral_operator(n):
    k = np.arange(n//2+1)
    retained = k <= n//3
    def evaluate(fields):
        coefficients = np.fft.rfft(fields, axis=-1)
        first = np.fft.irfft(1j*k*coefficients, n=n, axis=-1)
        second = np.fft.irfft(-k*k*coefficients, n=n, axis=-1)
        return first, second

    def filter_products(products):
        spectra = np.fft.rfft(products, axis=-1)
        spectra[..., ~retained] = 0
        return np.fft.irfft(spectra, n=n, axis=-1)

    return evaluate, filter_products


def finite_difference_operator(n):
    step = 2*math.pi/n

    def evaluate(fields):
        left, far_left = np.roll(fields, 1, axis=-1), np.roll(fields, 2, axis=-1)
        right, far_right = np.roll(fields, -1, axis=-1), np.roll(fields, -2, axis=-1)
        first = (far_left-8*left+8*right-far_right)/(12*step)
        second = (-far_left+16*left-30*fields+16*right-far_right)/(12*step**2)
        return first, second

    return evaluate, lambda products: products


def evolve(start, a, delays, discretization, tolerance):
    n = start.shape[-1]
    spatial, filter_products = (
        spectral_operator(n) if discretization == "fourier"
        else finite_difference_operator(n)
    )

    def rhs(t, flat):
        P, Pt, Q, Qt = flat.reshape(4,n)
        (Px, Qx), (Pxx, Qxx) = spatial(np.stack([P,Q]))
        P_source, Q_source = filter_products(np.stack([
            np.exp(2*P)*(Qt**2-Qx**2), -2*(Pt*Qt-Px*Qx)
        ]))
        return np.stack([Pt, Pxx-Pt/t+P_source, Qt, Qxx-Qt/t+Q_source]).ravel()

    ordered, indices = np.unique(delays, return_inverse=True)
    time_points = a+ordered
    solved = solve_ivp(rhs, (a,time_points[-1]), start.ravel(), t_eval=time_points,
                       method="DOP853", atol=0.025*tolerance, rtol=tolerance)
    if not solved.success:
        raise RuntimeError(solved.message)
    return solved.y.T[indices].reshape(len(delays),4,n)


def records(initial_P, evolved, spec):
    modes = np.arange(spec.K+1)
    initial_piece = np.cos(modes*spec.delays[spec.mode_slots])*initial_P
    selected = series(evolved[:,0])[spec.mode_slots,modes]
    return spec.weights*(selected-initial_piece)


def measure(start, spec, discretization, tolerance):
    measured = evolve(start,spec.a,spec.delays,discretization,tolerance)
    return records(series(start[0])[:spec.K+1], measured,spec)


def distributed_velocity(K, collision):
    rng = np.random.default_rng(171)
    coefficients = np.zeros(K+1, dtype=complex)
    coefficients[0] = -1/collision.a
    coefficients[collision.mode] = -0.5*collision.g*collision.amplitude/collision.a
    features = (3, 19, 4*collision.mode, 8*collision.mode, 11*collision.mode, K-3)
    for k in features:
        if k <= K:
            components = (rng.normal()+1j*rng.normal())*0.042/(1+k)
            coefficients[k] += components
    return coefficients


def decode(collision, spec, data, n, truth, tolerance, iterations=11):
    initial = collision.exact_fields(collision.a,theta_grid(n),0)
    guess = np.zeros(spec.K+1,dtype=complex)
    guess[0] = -1/collision.a
    history = []
    for step in range(iterations+1):
        initial[1] = to_grid(guess,n)
        predicted = measure(initial,spec,"finite_difference",tolerance)
        residual = data-predicted
        history.append({"iteration": step, "residual": h1_norm(residual),
                        "velocity_error": h1_norm(guess-truth)})
        if h1_norm(residual) < 5e-9:
            break
        guess += residual/spec.gains
    return history


def pair_experiment(collision, spec, n, tolerance):
    theta = theta_grid(n)
    delays = np.append(spec.delays,collision.gap)
    signed_initial = [collision.exact_fields(collision.a,theta,sgn) for sgn in (1.,-1.)]
    histories = [evolve(initial,collision.a,delays,"fourier",tolerance)
                 for initial in signed_initial]
    maximum_pde_error = max(
        float(np.max(np.abs(history[index]-collision.exact_fields(collision.a+d,theta,sgn))))
        for sgn,history in zip((1.,-1.),histories) for index,d in enumerate(delays)
    )
    signed_record = [records(series(initial[0])[:spec.K+1], history[:-1],spec)
                     for initial, history in zip(signed_initial,histories)]
    velocity_gap = h1_norm(series(signed_initial[0][1]-signed_initial[1][1]))
    mixed_gap = h1_norm(signed_record[0]-signed_record[1])
    endpoint_gap = float(np.max(np.abs(histories[0][-1,(0,2,3)]-histories[1][-1,(0,2,3)])))
    return {"n": n, "analytic_versus_independent_pde_max_difference": maximum_pde_error,
            "original_full_orbit_snapshot_gap": endpoint_gap, "initial_velocity_H1_gap": velocity_gap,
            "delayed_record_H1_gap": mixed_gap, "delayed_gain_on_exact_collision": mixed_gap/velocity_gap}


def campaign(K, reference_n, inverse_grids, tolerance, collision):
    spec = schedule(collision.a,collision.gap,K)
    baseline = schedule(collision.a,collision.gap,K,"fixed")
    collision_result = pair_experiment(collision,spec,reference_n,tolerance)
    truth = distributed_velocity(K,collision)
    initial = collision.exact_fields(collision.a,theta_grid(reference_n),0)
    initial[1] = to_grid(truth,reference_n)
    clean = measure(initial,spec,"fourier",tolerance)
    rng = np.random.default_rng(282)
    disturbance = rng.normal(size=K+1)+1j*rng.normal(size=K+1)
    disturbance[0] = disturbance[0].real
    disturbance *= 2e-5/h1_norm(disturbance)
    jitter = (0.03/K)*np.cos(0.7+np.arange(len(spec.delays)))
    shifted = retime(spec.a,spec.delays+jitter,spec.mode_slots)
    shifted_fields = evolve(initial,shifted.a,shifted.delays,"fourier",tolerance)
    unknown_jitter = records(series(initial[0])[:K+1],shifted_fields,spec)
    known_jitter = records(series(initial[0])[:K+1],shifted_fields,shifted)
    results = {"noiseless": {str(n):decode(collision,spec,clean,n,truth,tolerance) for n in inverse_grids},
               "additive_noise": decode(collision,spec,clean+disturbance,inverse_grids[-1],truth,tolerance),
               "unknown_clock_offsets": decode(collision,spec,unknown_jitter,inverse_grids[-1],truth,tolerance),
               "known_clock_offsets": decode(collision,shifted,known_jitter,inverse_grids[-1],truth,tolerance)}
    return {"bandwidth": K, "reference_grid": reference_n, "inverse_grids": inverse_grids,
            "slots": spec.delays.tolist(), "counts_per_slot": np.bincount(spec.mode_slots).tolist(),
            "free_minimum_gain_delayed": float(np.min(abs(spec.gains))),
            "free_minimum_gain_fixed": float(np.min(abs(baseline.gains))),
            "free_weakest_mode_fixed": int(np.argmin(abs(baseline.gains))),
            "clock_offsets": jitter.tolist(), "nominal_noise_H1": h1_norm(disturbance),
            "weighted_error_from_unknown_clock_H1": h1_norm(unknown_jitter-clean),
            "collision": collision_result, "decode": results}


if __name__ == "__main__":
    model = Collision()
    tol = 7e-12
    results = [campaign(96,768,(768,1536),tol,model),
               campaign(288,1536,(768,1536,3072),tol,model)]
    print(json.dumps({"platform": platform.platform(), "a": model.a,
                      "resonant_gap": model.gap, "characteristic_mode": model.mode,
                      "campaigns": results}, indent=2, allow_nan=False))
