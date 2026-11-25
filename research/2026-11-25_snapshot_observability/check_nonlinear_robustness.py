"""Different IVP, different forward stencil and perturbed observations."""
import json
import platform

import numpy as np

from check_nonlinear import benchmark, infer_velocity, load, use_fourth_order_rhs

TIME = 12.0
EARLIER = 10.0
TARGETS = (2.0, 4.0, 8.0)
NOISE = 1e-7


def vacuum_initial(n, model):
    x = 2.0*np.pi*np.arange(n)/n
    Q = 0.9*np.cos(x)+0.1*np.sin(2*x)
    R = 0.025*np.cos(3*x)
    Q_x = -0.9*np.sin(x)+0.2*np.cos(2*x)
    lam, mean = model.periodic_antiderivative(2.0*R*Q_x, model.dtheta_of(n))
    assert abs(mean) < 1e-14
    return [np.zeros(n), 2.0+1.4*np.cos(x), Q, R, lam]


if __name__ == "__main__":
    inverse = load("bb-cbe3/run_cbe3.py", "inverse_second_order")
    direct = load("bb-cbe3/run_cbe3.py", "truth_fourth_order")
    use_fourth_order_rhs(direct)
    output = {"platform": platform.platform(), "directional_noise": {},
              "distinct_IVP": {}}
    for n in (512, 1024):
        x = np.arange(n)*inverse.dtheta_of(n)
        dx = inverse.dtheta_of(n)
        end, _, tape = direct.evolve(vacuum_initial(n,direct), 0.0, TIME, n,
                                      store_taus=TARGETS+(EARLIER,))
        scored = benchmark(inverse, inverse, end, tape[EARLIER],
                           {t:tape[t] for t in TARGETS}, TIME, np.zeros(n,dtype=bool))
        scored["truth_end_momentum_residual_with_second_order_derivative"] = inverse.momentum_residual(end, TIME, dx)
        output["distinct_IVP"][str(n)] = scored

        if n == 1024:
            P, _, Q, R, _ = end
            zero_metric = np.zeros(n)
            observed = (tape[EARLIER][0],P,Q,R)
            _, nominal, _ = infer_velocity(inverse.evolve, TIME, observed)
            y_nominal, _, _ = inverse.evolve([P,nominal,Q,R,zero_metric], TIME, 2.0, n)
            for mode in (1,128):
                _, shifted, fit = infer_velocity(
                    inverse.evolve, TIME, (observed[0]+NOISE*np.cos(mode*x), P,Q,R))
                y_shifted, _, _ = inverse.evolve([P,shifted,Q,R,zero_metric], TIME, 2.0, n)
                diff = inverse.comparison_D(y_nominal,y_shifted,2.0,dx)
                output["directional_noise"][str(mode)] = {
                    "input_P_noise": NOISE,
                    "fit_residuals": fit,
                    "velocity_gain": float(np.max(np.abs(shifted-nominal))/NOISE),
                    "target_gain": float(np.max(diff)/NOISE),
                }
    print(json.dumps(output,indent=2,allow_nan=False))
