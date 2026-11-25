"""Weak modes of the two-full-metric map at a circular-loop transport resonance."""
import json
import math
import platform

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from check_both_velocity_holonomy import radial_background, vector_transfer


def resonant_modes(tolerance):
    a, winding, radius = 4.0, 4, 0.32
    limit = math.pi/winding
    radial = radial_background(a,a+limit,winding,radius,tolerance)
    def accumulated_connection(delay):
        return quad(lambda t: winding*math.cosh(radial(t)[0]),
                    a,a+delay,epsabs=2e-14,epsrel=2e-14)[0]
    delay = brentq(lambda d: accumulated_connection(d)-math.pi,0.8*limit,limit,
                   xtol=1e-14)
    c = a+delay
    modes = []
    for left,right in ((64,128),(128,256),(256,512),(512,1024),(1024,2048),(2048,4096)):
        k = np.arange(left,right)
        modes.append(int(k[np.argmin(abs(np.sin(k*delay)))]))
    k = np.array(modes)
    z, zt = vector_transfer(a,c,k,winding,radial,tolerance,return_terminal=True)
    R, Rt = radial(c)
    omega = winding*math.cosh(R)
    orbit = (math.sqrt(c/a)*np.sqrt(1+k*k))[:,None,None]*z
    terminal_row = -2*math.sqrt(c/a)*(winding*math.sinh(R)*zt[:,1]
        + Rt*(1j*k[:,None]*z[:,0]-omega*z[:,1]))
    initial_row = np.zeros_like(terminal_row)
    initial_row[:,1] = -2*winding*math.sinh(radius)
    full = np.concatenate((orbit,initial_row[:,None,:],terminal_row[:,None,:]),axis=1)
    gain = np.linalg.svd(full,compute_uv=False)[:,-1]
    assert abs(accumulated_connection(delay)-math.pi)<1e-12
    assert R < -0.2
    return {
        "initial_time": a, "initial_radius": radius, "terminal_signed_radius": float(R),
        "areal_delay": delay, "connection_angle": accumulated_connection(delay),
        "probes": [{"mode":m,"free_phase_gain":float(abs(np.sin(m*delay))),
                    "both_full_metrics_gain":float(g)} for m,g in zip(modes,gain)]
    }


if __name__ == "__main__":
    print(json.dumps({"platform":platform.platform(),
                      "resonant_loop":resonant_modes(5e-11)},indent=2,allow_nan=False))
