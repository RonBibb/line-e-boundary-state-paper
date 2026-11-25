"""Replay the regular-time holonomy and asymptotic consistency in the same family."""
import json
import math
import platform

import numpy as np
from scipy.integrate import simpson, solve_ivp

from check_nonlinear_sampling import evolve, spectral_operator, theta_grid

SPEED = 0.36
AREAL_START = 1.0
DENOMINATOR = 12
DELAY = math.pi / DENOMINATOR


def norms(values, order):
    coefficients = np.fft.fft(values, axis=-1) / values.shape[-1]
    k = np.fft.fftfreq(values.shape[-1]) * values.shape[-1]
    return float(np.sqrt(np.sum(abs(coefficients)**2 * (1+k*k)**order)))


def initial_fields(n, amplitude, perturbation=0.0):
    theta = theta_grid(n)
    P = amplitude * np.cos(theta)
    Q = amplitude * np.sin(theta)
    V = SPEED + perturbation * np.cos(2*theta)
    X = perturbation * np.sin(3*theta)
    return np.stack([P, -V/AREAL_START, Q, -np.exp(-P)*X/AREAL_START])


def holonomy(n, amplitude):
    delays = np.linspace(0.0, DELAY, 161)
    t = AREAL_START + delays
    P, _, Q, Qt = evolve(
        initial_fields(n, amplitude), AREAL_START, delays, "fourier", 8e-14
    ).transpose(1, 0, 2)
    differentiate, _ = spectral_operator(n)
    Qx, _ = differentiate(Q)
    omega_t, omega_theta = np.exp(P)*Qt, np.exp(P)*Qx
    k = np.arange(n//2+1)
    offset = (DELAY-delays)[:, None] * k[None, :]

    def shifted(values, distance):
        return np.fft.irfft(
            np.fft.rfft(values, axis=-1) * np.exp(1j*distance),
            n=n, axis=-1
        )

    diff = shifted(omega_t-omega_theta, offset) - shifted(omega_t+omega_theta, -offset)
    theta_angle = simpson(diff, x=t, axis=0)
    total = sum(np.roll(theta_angle, j*n//DENOMINATOR) for j in range(DENOMINATOR))

    def mode_equation(t, data):
        p, pt, q, qt = data
        return [pt, -p-pt/t, qt, -q-(1-2*SPEED)*qt/t]

    p1, _, q1, _ = solve_ivp(
        mode_equation, (AREAL_START, AREAL_START+DELAY), [1, 0, 1, 0],
        method="DOP853", t_eval=t, atol=2e-15, rtol=4e-14
    ).y
    prediction = -DENOMINATOR*simpson((t/AREAL_START)**(-SPEED)*p1*q1, x=t)
    return {
        "amplitude": amplitude,
        "predicted_angle_divided_by_amplitude_squared": prediction,
        "observed_mean_divided_by_amplitude_squared": float(np.mean(total)/amplitude**2),
        "smallest_absolute_angle": float(np.min(abs(total))),
        "angle_variation": float(np.ptp(total))
    }


def integrate_tail(data):
    n = data.shape[-1]
    differentiate, product_filter = spectral_operator(n)
    P_initial, Pt, Q_initial, Qt = data
    Px, _ = differentiate(P_initial)
    Qx, _ = differentiate(Q_initial)
    V_initial = -AREAL_START * Pt
    W_initial = AREAL_START * Px
    X_initial = -AREAL_START * np.exp(P_initial) * Qt
    Y_initial = AREAL_START * np.exp(P_initial) * Qx
    momentum = 2*(V_initial*Px+np.exp(P_initial)*X_initial*Qx)
    k = np.arange(n//2+1)
    integrated = np.fft.rfft(momentum)
    integrated[0] = 0
    integrated[1:] /= 1j*k[1:]
    lambda_initial = np.fft.irfft(integrated, n=n)
    lambda_initial -= lambda_initial[0]
    Pi_initial = np.exp(P_initial)*X_initial

    def rhs(r, flat):
        U, W, X, Y, P_residual, Q, _, _ = flat.reshape(8, n)
        P = P_residual + SPEED*r
        V = SPEED + U
        (Ux, Wx, Xx, Yx, Qx), _ = differentiate(np.stack((U,W,X,Y,Q)))
        e = AREAL_START*np.exp(-r)
        nonlinear = product_filter(np.stack((
            X*X-Y*Y, -U*X+W*Y, U*Y-W*X,
            np.exp(-P)*X, U*U+W*W+X*X+Y*Y,
            np.exp(2*(P-r))*AREAL_START**2*Qx
        )))
        derivative_flux, _ = differentiate(nonlinear[5])
        return np.stack((
            e*Wx+nonlinear[0], e*Ux-W, e*Yx-SPEED*X+nonlinear[1],
            e*Xx+(SPEED-1)*Y+nonlinear[2], U, nonlinear[3],
            2*SPEED*U+nonlinear[4], derivative_flux
        )).ravel()

    initial = np.stack((
        V_initial-SPEED, W_initial, X_initial, Y_initial,
        P_initial, Q_initial, lambda_initial, Pi_initial
    ))
    sample_times = np.array([0., 4., 8., 12., 16., 20., 24., 28.])
    solved = solve_ivp(
        rhs, (0.0, sample_times[-1]), initial.ravel(), t_eval=sample_times,
        method="DOP853", atol=6e-16, rtol=5e-13
    )
    if not solved.success:
        raise RuntimeError(solved.message)
    solved = solved.y.T.reshape(-1, 8, n)
    U, W, X, Y, P_residual, Q, lambda_residual, Pi = solved.transpose(1, 0, 2)
    v = SPEED + U
    p = P_residual - sample_times[:, None]*U
    lam = lambda_residual - sample_times[:, None]*(2*SPEED*U+U*U)
    (p_x, Qx, lam_x), _ = differentiate(np.stack((p[-1],Q[-1],lam[-1])))
    expected_momentum = 2*(v[-1]*p_x+Pi[-1]*Qx)
    state_records = np.stack([v,p,Q,lam,Pi], axis=1)
    orders = [2,2,3,2,1]

    def summarize(difference):
        return [norms(field,order) for field,order in zip(difference,orders)]

    return state_records, {
        "initial_H3_distance_to_vertical": norms(initial[:4],3),
        "max_normalized_H3_state_increase": max(norms(record[:4],3) for record in solved)
            / norms(initial[:4],3) - 1.0,
        "initial_momentum_constraint_defect": float(abs(np.mean(momentum))),
        "sample_tail_depths": sample_times.tolist(),
        "damped_H2_norms": [norms(record[1:4],2) for record in solved],
        "tail_errors_against_depth_28": [summarize(fields-state_records[-1])
                                         for fields in state_records[:-1]],
        "limiting_momentum_identity_H1_error": norms(lam_x-expected_momentum,1),
        "Pi_direct_versus_conservation_H1_error": norms(
            Pi[-1] - np.exp(P_residual[-1] + SPEED*sample_times[-1])*X[-1],1
        )
    }


def main():
    amplitude = 0.00025
    n = 192
    reference_data = initial_fields(n,amplitude)
    reference, report = integrate_tail(reference_data)
    coarser, coarse_report = integrate_tail(initial_fields(n//2,amplitude))
    orders = [2,2,3,2,1]
    delays = [0.02, DELAY]

    def measured(data):
        return evolve(data, AREAL_START, delays,"fourier",8e-14)[:,[0,2]]

    measured_reference = measured(reference_data)
    comparisons = []
    for perturbation in [0.00003, 0.000015, 0.0000075]:
        perturbed_data = initial_fields(n,amplitude,perturbation)
        observed_difference = norms(measured(perturbed_data)-measured_reference,4)
        perturbed, perturbed_report = integrate_tail(perturbed_data)
        output_differences = [norms(field,order) for field,order in zip(
            perturbed[-1]-reference[-1], orders
        )]
        comparisons.append({
            "velocity_perturbation_amplitude": perturbation,
            "two_later_snapshot_H4_difference": observed_difference,
            "asymptotic_differences": output_differences,
            "largest_individual_output_over_input_ratio": max(output_differences)/observed_difference,
            "initial_H3_distance_to_vertical": perturbed_report["initial_H3_distance_to_vertical"]
        })
    algebra_bound = 2**3*math.sqrt(1+2*(1+1/3))
    radius = min(SPEED,1-SPEED)/(16*algebra_bound)
    for entry in comparisons:
        assert entry["initial_H3_distance_to_vertical"] < radius
    assert report["initial_H3_distance_to_vertical"] < radius
    assert report["limiting_momentum_identity_H1_error"] < 1e-8
    holonomy_checks = [holonomy(n,epsilon) for epsilon in (0.00025,0.002,0.02)]
    for entry in holonomy_checks:
        assert entry["smallest_absolute_angle"] > 0
        assert abs(entry["predicted_angle_divided_by_amplitude_squared"]
                   - entry["observed_mean_divided_by_amplitude_squared"]) < 0.003
    print(json.dumps({
        "platform": platform.platform(),
        "reference_speed": SPEED,
        "analytic_H3_radius": radius,
        "reference_solution": report,
        "output_field_order": ["v","p","q","lambda_infinity","Pi_infinity"],
        "depth_28_grid_comparison": [norms(field,order) for field,order in zip(
            reference[-1,:,::2]-coarser[-1], orders
        )],
        "coarse_momentum_identity_H1_error": coarse_report["limiting_momentum_identity_H1_error"],
        "positive_angle_test": holonomy_checks,
        "perturbation_comparisons": comparisons
    }, indent=2))


if __name__ == "__main__":
    main()
