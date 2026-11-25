"""Finite numerical checks of the sharp polarized snapshot estimates."""
import json
import platform

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import gamma, j0, j1, jv, y0, y1, yv


def transfer(a, t, k):
    """Return C,S for y''+(k²+1/(4t²))y=0, (y(a), y'(a)/k)."""
    scale = np.sqrt(t / a)
    ja, ya = j0(k * a), y0(k * a)
    jt, yt = j0(k * t), y0(k * t)
    wa = np.pi * a / 2.0
    # h=sqrt(t)*f. The extra sqrt(t) derivative enters C.
    s = k * wa * scale * (ja * yt - ya * jt)
    c = wa * scale * k * (y1(k * a) * -jt + j1(k * a) * yt) - s / (2.0 * a * k)
    return c, s


def ode_check(a, b, k):
    analytic = np.array(transfer(a, b, k))
    numeric = []
    for y, z in ((1.0, 0.0), (0.0, 1.0)):
        result = solve_ivp(
            lambda t, state: [state[1], -(k**2 + 1.0 / (4.0 * t**2)) * state[0]],
            (a, b), (y, k * z), method="DOP853", atol=3e-13, rtol=3e-13,
        )
        assert result.success
        numeric.append(result.y[0, -1])
    return float(np.max(np.abs(analytic - numeric)))


def target_gains(a, b, target, k):
    cb, sb = transfer(a, b, k)
    ct, st = transfer(a, target, k)
    ratio = st / sb
    # Euclidean operator norm of one row acting on the weighted snapshots at a,b.
    return np.hypot(ct - ratio * cb, ratio)


if __name__ == "__main__":
    a, q = np.pi / 4, 16
    b = a + 3.0 * np.pi / q
    on, off = a + 4.0 * np.pi / q, a + 3.5 * np.pi / q
    area = lambda t: (1.0 / a - 1.0 / t) / 4.0
    assert 0.0 < area(b) < np.pi / q
    max_k = 65536
    modes = np.arange(1, max_k + 1)
    cb, sb = transfer(a, b, modes)
    on_gains = target_gains(a, b, on, modes)
    off_gains = target_gains(a, b, off, modes)
    # With y(a)=0, y'(a)=k, the Bessel coefficient is B=pi*k*sqrt(a)*J0(ka)/2.
    slope_gains = np.abs(modes * np.sqrt(a) * j0(modes * a) / sb)
    velocity_gains = np.hypot(cb, 1.0) / np.abs(sb)
    resonant = (modes % q == 0)
    odd_witness = modes[(modes % (2 * q) == q) & (modes > max_k // 2)]
    selected = odd_witness - 1
    ode_errors = {str(k): ode_check(a, b, k) for k in (1, 15, 16, 32, 49, 144)}
    assert max(ode_errors.values()) < 1e-10, ode_errors
    assert np.min(np.abs(modes * sb)) > 0.01
    leading = area(b) / 2.0
    asymptotic_error = np.max(np.abs(np.abs(modes[resonant][-32:] * sb[resonant][-32:]) - leading))
    assert asymptotic_error < 1e-5
    assert max(on_gains) < 8
    assert min(off_gains[selected] / odd_witness) > 1.0
    assert min(slope_gains[selected] / odd_witness**1.5) > 1.0
    band = np.floor(np.log2(modes)).astype(int)
    d = np.pi / (4.0 * (2.0 + 1.0 / (4.0 * a**2)))
    band_time = a + d * np.exp2(-band)
    band_c, band_s = transfer(a, band_time, modes)
    band_velocity_gain = np.hypot(band_c, 1.0) / np.abs(band_s)
    assert np.min(np.abs(band_s)) > 0.2
    assert np.max(band_velocity_gain) < 7.0
    sizes = np.array([16, 48, 80, 144, 272, 528, 1040, 2064, 4112, 8208, 16400, 32784, 65520])
    c_off, s_off = transfer(a, off, sizes)
    band_ratio = s_off / band_s[sizes-1]
    off_band_gain = np.hypot(c_off - band_ratio * band_c[sizes-1], band_ratio)
    epd = {}
    for alpha in (0.4, 1.6):
        nu = (alpha-1.0)/2.0
        ks = np.array([1024,4096,16384,65536])
        s_b = (np.pi/2.0*ks*np.sqrt(a*b)
               *(jv(nu,ks*a)*yv(nu,ks*b)-yv(nu,ks*a)*jv(nu,ks*b)))
        # Witness y(a)=0, y'(a)=k. Only the divergent/subleading power
        # coefficient of Y_nu contributes to the limiting flux.
        b_coeff = np.pi/2.0*np.sqrt(a)*ks*jv(nu,ks*a)
        power_coeff = -gamma(nu)*2.0**nu/np.pi * b_coeff * ks**(-nu)
        gain = np.abs(-2.0*nu*power_coeff/s_b)
        slope = float(np.log(gain[-1]/gain[-2])/np.log(ks[-1]/ks[-2]))
        assert abs(slope-(2-alpha/2.0)) < 1e-4
        epd[str(alpha)] = {"predicted_power": 2-alpha/2.0, "observed_power": slope}
    timing = {}
    for k in (128, 512, 2048, 8192):
        expected = leading / k**2
        shifted_b = brentq(lambda t: float(transfer(a, t, k)[1]),
                           b - 2*expected, b, xtol=5e-16)
        timing[str(k)] = {
            "k_squared_times_sample_shift": float(k**2*(b-shifted_b)),
            "abs_S_at_blind_time": abs(float(transfer(a, shifted_b, k)[1])),
            "stable_target_gain_with_k_minus_2_shift": float(target_gains(a,b,on+1/k**2,k)),
            "stable_target_gain_with_k_minus_1_shift": float(target_gains(a,b,on+1/(2*k),k)),
        }
        assert abs(timing[str(k)]["k_squared_times_sample_shift"] - leading) < 1e-4
    # For deliberately resonant sample times, an exact kernel can occur. Locate
    # a root of S_20 on a bracket near the first zero of the leading sine.
    sample_a, sample_k = 0.85, 20
    near = sample_a + np.pi / sample_k
    sample_b = brentq(lambda t: float(transfer(sample_a, t, sample_k)[1]),
                      near - 0.01, near + 0.01, xtol=1e-15)
    sample_target = sample_a + np.pi / (2 * sample_k)
    assert abs(transfer(sample_a, sample_b, sample_k)[1]) < 1e-12
    assert abs(transfer(sample_a, sample_target, sample_k)[1]) > 0.99
    print(json.dumps({
        "platform": platform.platform(), "a": a, "b": b, "on_target": on,
        "off_target": off, "condition": area(b), "condition_threshold": np.pi / q,
        "modes_checked": max_k, "ode_max_errors": ode_errors,
        "zero_mode": "f=A+B log(t); y(a), y(b) determine A and B for a!=b",
        "min_k_abs_S_b": float(np.min(np.abs(modes * sb))),
        "asymptotic_k_abs_S_b": leading,
        "observed_error_in_last_32_resonant_modes": float(asymptotic_error),
        "on_target_max_gain": float(max(on_gains)),
        "off_target_gain_by_mode": {str(k): float(off_gains[k-1]) for k in (16, 144, 1008, 10000, 65520)},
        "min_off_target_gain_per_k_on_upper_witnesses": float(min(off_gains[selected] / odd_witness)),
        "min_velocity_gain_per_k_on_upper_witnesses": float(min(velocity_gains[selected] / odd_witness)),
        "min_singularity_slope_gain_per_k_pow_1p5_on_upper_witnesses": float(min(slope_gains[selected] / odd_witness**1.5)),
        "gain_curve": [
            {"k": int(k), "fixed_aligned": float(on_gains[k-1]),
             "fixed_offset": float(off_gains[k-1]), "band_offset": float(gain)}
            for k, gain in zip(sizes, off_band_gain)
        ],
        "band_schedule": {
            "distinct_additional_times": int(len(np.unique(band))),
            "maximum_velocity_gain": float(np.max(band_velocity_gain)),
            "minimum_abs_S": float(np.min(np.abs(band_s))),
        },
        "epd_flux_check": epd,
        "timing_perturbation_checks": timing,
        "exact_resonance_numerical_example": {
            "k": sample_k, "a": sample_a, "b": sample_b,
            "S_b": float(transfer(sample_a, sample_b, sample_k)[1]),
            "target": sample_target, "S_target": float(transfer(sample_a, sample_target, sample_k)[1]),
        },
    }, indent=2, allow_nan=False))
