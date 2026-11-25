"""Two exact geodesic-image Gowdy histories with identical mixed data."""
import json
import math
import platform

import mpmath as mp
import numpy as np
from scipy.optimize import brentq
from scipy.special import j0, j1, y0, y1

from check_obstruction import bessel, collision


def collision_fields():
    fields=collision(3.0,0.95)
    k=4
    epsilon=0.02
    ratio=math.exp(fields["log_gap"])
    # First positive Dirichlet eigenvalue for the two-radii determinant.
    determinant=lambda z:j0(z)*y0(ratio*z)-y0(z)*j0(ratio*z)
    nodes=np.linspace(1e-4,0.5,1024)
    ends=[(lo,hi) for lo,hi in zip(nodes,nodes[1:]) if determinant(lo)*determinant(hi)<0]
    assert ends
    z=brentq(determinant,*ends[0],xtol=5e-16)
    a=z/k
    b=ratio*a
    theta=np.linspace(0,2*math.pi,512,endpoint=False)
    p_wave=np.cos(k*theta)
    def evaluated(t,guess):
        profile=(math.pi/2)*(j0(z)*y0(k*t)-y0(z)*j0(k*t))
        deriv=(math.pi/2)*k*(-j0(z)*y1(k*t)+y0(z)*j1(k*t))
        w=math.log(t/a)+epsilon*profile*p_wave
        ws=1.0+epsilon*t*deriv*p_wave
        g,v=guess["initial_P_s"],guess["speed"]
        u=np.cosh(v*w)+(g/v)*np.sinh(v*w)
        P=np.log(u)
        Q=fields["fixed_q"]*np.sinh(v*w)/(v*u)
        Qs=fields["fixed_q"]*ws/u**2
        return P,Q,Qs
    left=[evaluated(t,fields["branches"][0]) for t in (a,math.sqrt(a*b),b)]
    right=[evaluated(t,fields["branches"][1]) for t in (a,math.sqrt(a*b),b)]
    measured=max(*(float(np.max(np.abs(x-y))) for x,y in zip(left[0],right[0])),
                 float(np.max(np.abs(left[-1][0]-right[-1][0]))))
    gap=float(np.max(np.abs(left[1][0]-right[1][0])))
    assert measured<1e-12 and gap>0.3
    assert max(v["speed"] for v in fields["branches"])*(1+epsilon)<1
    return {"a":a,"b":b,"k":k,"epsilon":epsilon,"ratio":ratio,"determinant":float(determinant(z)),
            "two_snapshot_observation_max_difference":measured,
            "interior_P_max_difference":gap,"initial_Qs_peak_to_peak":float(np.ptp(left[0][2])),
            "terminal_Q_difference":float(np.max(np.abs(left[-1][1]-right[-1][1]))),
            "uniform_asymptotic_speed_bounds":[
                [row["speed"]*(1-abs(epsilon*j0(z))),row["speed"]*(1+abs(epsilon*j0(z)))]
                for row in fields["branches"]
            ]}


def full_metric_collision():
    ratio = math.e
    k, epsilon = 4, 0.02
    P0, Q0 = 0.2, 0.31
    g, q = 0.55, 0.37
    sigma = math.hypot(g, q)

    def determinant(z):
        return j0(z)*y0(ratio*z) - y0(z)*j0(ratio*z)

    nodes = np.linspace(1e-4, 4*math.pi/(ratio-1), 1024)
    bracket = next((lo, hi) for lo, hi in zip(nodes, nodes[1:])
                   if determinant(lo)*determinant(hi) < 0)
    z = brentq(determinant, *bracket, xtol=5e-16)
    a, b = z/k, ratio*z/k
    c = (a+b)/2
    theta = np.linspace(0, 2*math.pi, 512, endpoint=False)
    cosine, sine = np.cos(k*theta), np.sin(k*theta)

    def evaluated(t, initial_P_velocity):
        profile = (math.pi/2)*(j0(z)*y0(k*t) - y0(z)*j0(k*t))
        deriv = (math.pi/2)*k*(-j0(z)*y1(k*t) + y0(z)*j1(k*t))
        w = epsilon*profile*cosine
        w_t = epsilon*deriv*cosine
        w_theta = -epsilon*k*profile*sine
        u = np.cosh(sigma*w) + initial_P_velocity*np.sinh(sigma*w)/sigma
        P = P0 + np.log(u)
        Q = Q0 + math.exp(-P0)*q*np.sinh(sigma*w)/(sigma*u)
        P_w = (sigma*np.sinh(sigma*w) + initial_P_velocity*np.cosh(sigma*w))/u
        Q_w = math.exp(-P0)*q/u**2
        speed2 = P_w**2 + (np.exp(P)*Q_w)**2
        lambda_t = -t*speed2*(w_t**2+w_theta**2)
        lambda_theta = -2*t*speed2*w_t*w_theta
        return P, Q, P_w*w_t, Q_w*w_t, lambda_t, lambda_theta

    left = [evaluated(t, g) for t in (a, c, b)]
    right = [evaluated(t, -g) for t in (a, c, b)]
    def difference(at, field):
        return float(np.max(np.abs(left[at][field] - right[at][field])))

    measured = max(difference(at, field) for at in (0, 2) for field in (0, 1, 3))
    constraint_difference = max(difference(at, field) for at in range(3) for field in (4, 5))
    assert measured < 1e-12 and constraint_difference < 1e-13
    assert difference(0, 2) > 0.01 and difference(1, 0) > 0.001
    return {"a": a, "b": b, "interior_sample": c, "k": k, "epsilon": epsilon,
            "ratio": ratio, "determinant": float(determinant(z)),
            "both_snapshots_P_Q_Qt_max_difference": measured,
            "all_sampled_metric_gradient_max_difference": constraint_difference,
            "initial_Pt_max_difference": difference(0, 2),
            "third_snapshot_P_max_difference": difference(1, 0),
            "uniform_asymptotic_speed_upper_bound": sigma*epsilon*abs(j0(z))}


def safe_snapshot(a_reference, b_reference):
    mp.mp.dps = 60
    a, b = mp.mpf(a_reference), mp.mpf(b_reference)
    bound = 1/(4*a*a)
    L = int(mp.floor(max(mp.pi/(b-a), 2*bound+2)))+1
    c = a+mp.pi/L
    assert a < c < b
    limit = (1/a-1/c)/8
    resonant = [(k, float(k*abs(bessel(a,c,mp.mpf(1)/4,k))))
                for k in (L, 10*L, 100*L, 1000*L)]
    assert abs(resonant[-1][1]-float(limit)) < 2e-7
    return {"a": float(a), "b": float(b), "c": float(c), "L": L,
            "predicted_resonant_limit": float(limit),
            "k_times_absolute_transfer": dict(resonant)}


def detuning():
    mp.mp.dps=65
    a,b=mp.pi/4,7*mp.pi/16
    u,v=1/a,1/b
    C=(u*u+v*v)/6+u*v/24
    rows={}
    for k in (16,64,256,1024,4096,16384):
        A=bessel(a,b,mp.mpf(1)/4,k)
        B=bessel(a,b,-mp.mpf(1)/4,k)
        delta=-(A+B)/A
        rows[str(k)]={"A_over_B_shift":float(delta),"k2_times_shift":float(k*k*delta),
                      "detuned_scaled_null":float(abs(((1+delta)*A+B)*k**3))}
    assert abs(rows["16384"]["k2_times_shift"]-float(C))<1e-8
    return {"predicted_limit":float(C),"rows":rows}


if __name__ == "__main__":
    full = full_metric_collision()
    print(json.dumps({"platform": platform.platform(),
                      "inhomogeneous_exact_pair": collision_fields(),
                      "full_metric_collision": full,
                      "added_snapshot": safe_snapshot(full["a"], full["b"]),
                      "linearized_detunings": detuning()}, indent=2, allow_nan=False))
