"""Checks of small-amplitude recovery, metric regularity, and twisted examples."""
from dataclasses import replace
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np
from scipy.special import j0,j1,y0,y1
from scipy.optimize import brentq

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from check_obstruction import collision
from check_nonlinear_sampling import Collision, theta_grid
from check_null_transport import tangent_at,l2


def rank_two_amplitude():
    a, k0, k = 1.6, 4, 512
    d = math.pi/16
    n = 4*k+128
    theta = theta_grid(n)
    potential = np.sin(k*theta)
    norm = l2(potential)
    rows = {}
    reference = None
    for amplitude in (0.,0.02,0.04,0.08,0.16,0.32):
        background = Collision(a=a,mode=k0,amplitude=amplitude)
        tangent = tangent_at(background,(d,),(k*np.cos(k*theta))[None,:],7e-11)[0]
        P = background.exact_fields(a+d,theta,1.)[0]
        p = np.sqrt((a+d)/a)*tangent[0,0]
        q = np.sqrt((a+d)/a)*np.exp(P)*tangent[2,0]
        if amplitude==0:
            reference = p
            continue
        phase,_ = background.radial_phase(a+d)
        expected = background.q*phase*np.cos(k0*theta)*potential
        rows[str(amplitude)] = {
            "P_difference_over_amplitude_squared":l2(p-reference)/(norm*amplitude**2),
            "Q_gain_over_amplitude":l2(q)/(norm*amplitude),
            "Q_principal_linear_coefficient":l2(expected)/norm,
            "full_orbit_gain":l2(np.stack([p,q]))/norm}
    return {"a":a,"delay":d,"seed_mode":k0,"resonant_probe":k,
            "rank_two_collision_later_delay":Collision(a=a,mode=k0).gap,"rows":rows}


def modulus_endpoint():
    a,b = math.pi/4,7*math.pi/16
    Kmax = 131072
    modes = np.arange(16,Kmax+1,16)
    z0,z1 = modes*a,modes*b
    det = j0(z0)*y0(z1)-y0(z0)*j0(z1)
    coeff_A = -y0(z0)/det
    coeff_B = j0(z0)/det
    intercept = coeff_A+2*coeff_B*(np.log(modes/2)+np.euler_gamma)/math.pi
    v = 0.65
    multiplier = 2*v*intercept
    norm_squares = {}
    for alpha in (0.,1.,1.5,2.):
        weighted = 2*multiplier**2/((1+modes**2)**2*np.log(2+modes)**(2*alpha))
        running = np.cumsum(weighted)
        norm_squares[str(alpha)] = {str(K):float(running[np.searchsorted(modes,K,"right")-1])
                                    for K in (256,1024,4096,16384,65536,131072)}
    return {"starts":[a,b],"base_velocity":v,
            "mode_multiplier_over_k_to_three_halves_log_k": {
                str(k):float(multiplier[np.searchsorted(modes,k)]/(k**1.5*math.log(k)))
                for k in (256,4096,65536,131072)},
            "truncated_linearized_modulus_operator_norm_squared":norm_squares}


def twisted_collision():
    branches = collision(3.,0.95)
    q = branches["fixed_q"]
    k, eps = 4, 0.02
    ratio = math.exp(branches["log_gap"])
    current, P0 = 1/(2*math.pi), math.log(1/(2*math.pi*q))
    determinant = lambda z:j0(z)*y0(ratio*z)-y0(z)*j0(ratio*z)
    z = brentq(determinant,0.14,0.16,xtol=1e-15)
    a,b = z/k,ratio*z/k
    theta = theta_grid(1024)

    def state(t, branch, angular):
        u = math.pi/2*(j0(z)*y0(k*t)-y0(z)*j0(k*t))
        up = math.pi*k/2*(-j0(z)*y1(k*t)+y0(z)*j1(k*t))
        w = math.log(t/a)+eps*u*np.cos(k*angular)
        wt = 1/t+eps*up*np.cos(k*angular)
        g,sigma = branch["initial_P_s"],branch["speed"]
        factor = np.cosh(sigma*w)+g*np.sinh(sigma*w)/sigma
        wp = (sigma*np.sinh(sigma*w)+g*np.cosh(sigma*w))/factor
        P = -P0-np.log(factor)-math.log(t)
        Q = current*(angular+eps*t*up*np.sin(k*angular)/k)
        Pt = -wp*wt-1/t
        Qt = -current*t*eps*k*u*np.sin(k*angular)
        return np.stack([P,Q,Pt,Qt])

    left,right = ([state(t,branch,theta) for t in (a,math.sqrt(a*b),b)]
                  for branch in branches["branches"])
    record_gap = max(float(np.max(np.abs(left[i][part]-right[i][part])))
                     for i,part in ((0,0),(0,1),(0,3),(2,0)))
    monodromy = state(a,branches["branches"][0],theta+2*math.pi)[1]-left[0][1]
    bounds = [[1-x["speed"]*(1+abs(eps*j0(z))),1-x["speed"]*(1-abs(eps*j0(z)))]
              for x in branches["branches"]]
    assert record_gap < 1e-12 and min(row[0] for row in bounds)>0
    assert max(row[1] for row in bounds)<1
    return {"a":a,"b":b,"marking_shift":float(np.mean(monodromy)),
            "marking_shift_max_error":float(np.max(abs(monodromy-1))),
            "maximum_four_observable_gap":record_gap,
            "initial_Pt_gap":float(np.max(np.abs(left[0][2]-right[0][2]))),
            "strict_asymptotic_speed_ranges":bounds}


if __name__ == "__main__":
    output = {"platform":platform.platform(),"rank_two_amplitude":rank_two_amplitude(),
              "modulus_endpoint":modulus_endpoint(),"twisted_collision":twisted_collision()}
    print(json.dumps(output,indent=2,allow_nan=False))
