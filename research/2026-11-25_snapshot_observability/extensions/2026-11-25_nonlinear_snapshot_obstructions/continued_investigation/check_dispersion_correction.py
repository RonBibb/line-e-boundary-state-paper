"""Separate principal stencil dispersion from nonlinear residual model error."""
import json
import math
import platform

import numpy as np
from scipy.special import j0,j1,y0,y1

from check_nonlinear_sampling import (Collision,distributed_velocity,measure,
                                      schedule,series,theta_grid,to_grid,h1_norm)


def fourth_order_frequencies(k,n):
    h = 2*math.pi/n
    sine = np.sin(k*h/2)
    return 2*sine/h*np.sqrt(1+sine*sine/3)


def damped_transfer(a, later, frequencies):
    C,S = np.ones_like(frequencies,dtype=float),a*np.log(later/a)
    mask = frequencies>0
    k = frequencies[mask]
    begin,end = k*a,k*later[mask]
    det = 2/(math.pi*a)
    C[mask] = k*(-j0(end)*y1(begin)+y0(end)*j1(begin))/det
    S[mask] = (-j0(end)*y0(begin)+y0(end)*j0(begin))/det
    return C,S


def initial_state(problem,n,case):
    angular = theta_grid(n)
    if case=="rank_two":
        return problem.exact_fields(problem.a,angular,0.)
    return np.stack([
        -0.5+0.37*np.cos(3*angular)+0.12*np.cos(17*angular),
        np.full(n,-1/problem.a),
        0.2+0.20*np.cos(4*angular)+0.10*np.cos(27*angular),
        0.5+0.42*np.cos(2*angular)+0.06*np.cos(29*angular),
    ])


def execute(K,reference_grid,inverse_grid,case):
    problem = Collision()
    spec = schedule(problem.a,problem.gap,K)
    record_times = spec.a+spec.delays[spec.mode_slots]
    true_velocity = distributed_velocity(K,problem)
    reference = initial_state(problem,reference_grid,case)
    reference[1] = to_grid(true_velocity,reference_grid)
    record = measure(reference,spec,"fourier",9e-12)
    initial = initial_state(problem,inverse_grid,case)
    initial_P = series(initial[0])[:K+1]
    modes = np.arange(K+1)
    discrete = fourth_order_frequencies(modes,inverse_grid)
    continuum_C,continuum_S = damped_transfer(spec.a,record_times,modes)
    discrete_C,discrete_S = damped_transfer(spec.a,record_times,discrete)
    linear_compatibility = spec.weights*continuum_S
    discrete_gain = spec.weights*discrete_S
    def predict(guess, corrected):
        initial[1] = to_grid(guess,inverse_grid)
        uncorrected = measure(initial,spec,"finite_difference",9e-12)
        if not corrected:
            return uncorrected
        return uncorrected+spec.weights*((continuum_C-discrete_C)*initial_P+
                                         (continuum_S-discrete_S)*guess)
    histories = {}
    for method in ("unmodified","dispersion_corrected"):
        corrected = method=="dispersion_corrected"
        guess = np.zeros(K+1,dtype=complex)
        guess[0] = -1/spec.a
        history = []
        for iteration in range(10):
            residual = record-predict(guess,corrected)
            error = h1_norm(guess-true_velocity)
            history.append({"iteration":iteration,"record_residual":h1_norm(residual),
                            "velocity_error":error})
            if h1_norm(residual)<1e-9 or error>2:
                break
            guess += residual/linear_compatibility
        histories[method] = history
    return {"case":case,"K":K,"reference_grid":reference_grid,"inverse_grid":inverse_grid,
            "minimum_continuum_linear_gain":float(np.min(abs(linear_compatibility))),
            "minimum_discrete_linear_gain":float(np.min(abs(discrete_gain))),
            "maximum_linear_gain_discrepancy":float(np.max(abs(
                linear_compatibility-discrete_gain))),
            "raw_model_discrepancy_at_truth":h1_norm(record-predict(true_velocity,False)),
            "corrected_model_discrepancy_at_truth":h1_norm(record-predict(true_velocity,True)),
            "histories":histories}


if __name__=="__main__":
    result = {"platform":platform.platform(),"runs":[execute(K,r,n,case)
               for K,r,n,case in ((96,768,768,"rank_two"),(288,1536,768,"rank_two"),
                                  (288,1536,1536,"rank_two"),(288,1536,3072,"rank_two"),
                                  (288,1536,768,"coupled"),(288,1536,1536,"coupled"))]}
    print(json.dumps(result,indent=2,allow_nan=False))
