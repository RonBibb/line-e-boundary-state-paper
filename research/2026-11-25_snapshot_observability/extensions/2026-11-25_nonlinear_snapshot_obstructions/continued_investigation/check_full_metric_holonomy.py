"""Compare exact circular-loop orbit and metric transfers at an irrational delay."""
import json
import math
import platform

import mpmath as mp
import numpy as np
from scipy.integrate import quad

from check_both_velocity_holonomy import radial_background, vector_transfer


def compare(radius, tolerance, d=math.pi*math.sqrt(2)/24):
    a, winding = 1.6, 4
    c = a+d
    radial = radial_background(a,c,winding,radius,tolerance)
    terminal_radius, radial_rate = radial(c)
    omega = winding*math.cosh(terminal_radius)
    angular_speed = winding*math.sinh(terminal_radius)
    rotation = quad(lambda t: winding*math.cosh(radial(t)[0]), a,c,epsabs=1e-13)[0]
    metric_gramian_eigenvalue = angular_speed**2*(1-abs(math.cos(2*rotation)))
    selected = []
    for start,stop in ((128,256),(256,512),(512,1024),(1024,2048),(2048,4096)):
        bands = np.arange(start,stop)
        principal = np.minimum(abs(np.sin(bands*d+rotation)),
                               abs(np.sin(bands*d-rotation)))
        selected.append(int(bands[np.argmin(principal)]))
    modes = np.array([*range(129),*selected])
    z, zt = vector_transfer(a,c,modes,winding,radial,tolerance,return_terminal=True)
    norm_weight = math.sqrt(c/a)
    orbit = norm_weight*np.sqrt(1+modes**2)[:,None,None]*z
    gradient = -2*norm_weight*(angular_speed*zt[:,1]
        + radial_rate*(1j*modes[:,None]*z[:,0]-omega*z[:,1]))
    full_metric = np.concatenate((orbit,gradient[:,None,:]),axis=1)
    all_orbit = np.linalg.svd(orbit,compute_uv=False)[:,-1]
    all_metric = np.linalg.svd(full_metric,compute_uv=False)[:,-1]
    assert metric_gramian_eigenvalue > 0
    assert min(all_metric) > min(all_orbit)*4
    return {
        "circular_radius": radius, "areal_delay": d,
        "gramian_minimum_eigenvalue": metric_gramian_eigenvalue,
        "minimum_orbit_gain_through_128": float(min(all_orbit[:129])),
        "minimum_metric_gain_through_128": float(min(all_metric[:129])),
        "high_modes": [
            {"mode": int(k), "orbit_gain": float(o), "full_metric_gain": float(m)}
            for k,o,m in zip(modes[129:],all_orbit[129:],all_metric[129:])
        ]
    }


def certify_short_delay_constants():
    """Outward-rounded bounds for the all-mode reverse estimate."""
    mp.iv.dps = 45
    interval = mp.iv.mpf
    initial_radius, initial_time = interval("0.32"), interval("1.6")
    max_delay = interval("0.02")
    positive, negative = mp.iv.exp(initial_radius), mp.iv.exp(-initial_radius)
    speed = 2*(positive-negative)
    connection = 2*(positive+negative)
    angular_variation_per_delay = speed**2
    minimum_radius = initial_radius-speed*max_delay
    potential = (2*connection*angular_variation_per_delay*max_delay
                 + speed**2 + 1/(4*initial_time**2))
    slope = 2*angular_variation_per_delay
    perturbation_factor = 2*(potential+slope)
    correction = 2*perturbation_factor*max_delay*mp.iv.exp(
        perturbation_factor*max_delay**2
    )
    free_gain = 4*(1-(connection*max_delay)**2/6)/(1+connection*max_delay)
    assert minimum_radius.a > interval("0.29")
    assert speed.b < interval("1.31")
    assert connection.b < interval("4.21")
    assert potential.b < interval("2.11")
    assert slope.b < interval("3.44")
    assert correction.b < interval("0.45")
    assert free_gain.a > interval("3.68")
    return {
        "certifies_delays_up_to": "0.02",
        "radius_lower_interval": str(minimum_radius),
        "potential_upper_interval": str(potential),
        "connection_perturbation_slope_upper_interval": str(slope),
        "reverse_error_over_delay_upper_interval": str(correction),
        "frozen_gain_over_delay_lower_interval": str(free_gain)
    }


if __name__ == "__main__":
    print(json.dumps({
        "platform": platform.platform(),
        "short_delay_interval_check": certify_short_delay_constants(),
        "normalization": "sqrt(c/a) on orbit, sqrt(c/a)/c on metric gradient",
        "experiments": [compare(radius,2e-10) for radius in (0.08,0.16,0.32)],
        "strict_tolerance_recheck": compare(0.32,2e-12),
        "explicit_short_delay": compare(0.32,2e-12,d=0.02)
    },indent=2,allow_nan=False))
