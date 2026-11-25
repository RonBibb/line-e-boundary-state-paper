"""Checks of exact homogeneous collisions and linearized resonances."""
import json
import math
import platform

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def homogeneous(ell, g, q):
    speed = math.hypot(g, q)
    h = speed * ell
    if speed == 0:
        return 0.0, 0.0, 0.0
    u = math.cosh(h) + (g / speed) * math.sinh(h)
    return math.log(u), q * math.sinh(h) / (speed * u), -(speed**2) * ell


def collision(ell, max_speed):
    h = max_speed * ell
    left = -max_speed * math.tanh(h / 2)
    q = max_speed / math.cosh(h / 2)
    zero = lambda g: homogeneous(ell, g, q)[0]
    assert abs(zero(left)) < 2e-14
    assert zero(left + 1e-5) < -1e-8
    right = brentq(zero, left + 1e-5, 0.0, xtol=1e-14)
    assert right > left
    rows = []
    for g in (left, right):
        p, Q, lam = homogeneous(ell, g, q)
        def rhs(s, state):
            P, Pt, _, Qt, unused = state
            return (Pt, math.exp(2*P)*Qt*Qt, Qt, -2*Pt*Qt,
                    -(Pt*Pt + math.exp(2*P)*Qt*Qt))
        solution = solve_ivp(rhs, (0.0,ell), (0.0,g,0.0,q,0.0),
                             method="DOP853", rtol=1e-13, atol=1e-13)
        assert solution.success
        independent = solution.y[(0,2,4),-1]
        difference = float(np.max(np.abs(independent-np.array([p,Q,lam]))))
        assert difference < 1e-10
        rows.append({"initial_P_s": g, "speed": math.hypot(g,q), "P_b": p,
                     "Q_b": Q, "lambda_b_minus_a": lam,
                     "coordinate_ODE_error": difference})
    return {"log_gap": ell, "fixed_q": q, "branches": rows}


def bessel(a, b, beta, k):
    order = mp.sqrt(mp.mpf(1)/4-beta)
    return (mp.pi*k/2 * mp.sqrt(a*b)
            * (mp.besselj(order,k*a)*mp.bessely(order,k*b)
               -mp.bessely(order,k*a)*mp.besselj(order,k*b)))


def weighted_channels(a,b,speed,center,k):
    r0, r1 = speed*mp.log(a/center), speed*mp.log(b/center)
    tangent = mp.tanh(r0)*mp.tanh(r1)
    normal = 1/(mp.cosh(r0)*mp.cosh(r1))
    first = bessel(a,b,mp.mpf(1)/4,k) if k else None
    second = bessel(a,b,mp.mpf(1)/4-speed**2,k) if k else None
    return tangent,normal,first,second


def coordinates(a,b,speed,center,k):
    def rhs(t, v):
        r = speed*math.log(t/center)
        P = math.log(math.cosh(r))
        Pt = speed/t*math.tanh(r)
        Qt = speed/(t*math.cosh(r)**2)
        p,p_t,Q,Q_t = v
        return (p_t, -p_t/t-k*k*p + 2*math.exp(2*P)*Qt*(Qt*p+Q_t),
                Q_t, -Q_t/t-k*k*Q-2*(Pt*Q_t+Qt*p_t))
    solution = solve_ivp(rhs,(a,b),(0.0,1.0,0.0,0.0),method="DOP853",
                         rtol=2e-13,atol=2e-14)
    assert solution.success
    return solution.y[[0,2],-1]


def fourier_checks():
    mp.mp.dps = 65
    a = mp.pi/4
    b = 7*mp.pi/16
    ell = mp.log(b/a)
    speed = mp.sqrt(mp.mpf('0.5'))
    h = speed * ell
    # cosh(r_a+r_b)=2+cosh(h) sets sinh(r_a)sinh(r_b)=1.
    r_a = (mp.acosh(2+mp.cosh(h))-h)/2
    center = a*mp.exp(-r_a/speed)
    u,v = 1/a,1/b
    tangent,normal,_,_ = weighted_channels(a,b,speed,center,1)
    assert abs(tangent-normal) < mp.mpf('1e-60')
    c32 = -(u-v)*((u*u+v*v)/6+u*v/24)
    predicted = -2*normal*(mp.mpf(1)/4)**2*c32
    asymptotic = {}
    for k in (16,64,256,1024,4096,16384):
        _,_,first,second = weighted_channels(a,b,speed,center,k)
        val = tangent*first+normal*second
        # Tangential and normal output recovered from both target coordinates.
        speed0 = speed*mp.log(a/center)
        full_norm = mp.sqrt((mp.tanh(speed0)*first)**2
                            +(mp.sech(speed0)*second)**2)
        asymptotic[str(k)] = {"abs_partial_k3": float(abs(val)*k**3),
                             "full_target_norm_k": float(full_norm*k)}
    assert abs(asymptotic["16384"]["abs_partial_k3"]-float(predicted)) < 1e-9
    direct_rows = {}
    for s,t0,name in ((speed,center,"cancelled"),(mp.mpf('0.5'),a,"exact_kernel")):
        for k in (1,15,16,32,49):
            A,B,F,G = weighted_channels(a,b,s,t0,k)
            # e_P(a)=tanh(r_a) E_parallel - sech(r_a) E_normal.
            ra,rb = s*mp.log(a/t0),s*mp.log(b/t0)
            E_parallel = mp.tanh(ra)*F
            E_normal = -mp.sech(ra)*G
            expected = [mp.sqrt(a/b)/k*(A*F+B*G),
                        mp.sqrt(a/b)/k/mp.cosh(rb)
                        *(mp.sech(rb)*E_parallel+mp.tanh(rb)*E_normal)]
            direct = coordinates(float(a),float(b),float(s),float(t0),k)
            err = float(np.max(np.abs(direct-np.array(expected,dtype=float))))
            assert err < 5e-12, (name,k,err)
            direct_rows[f"{name}_{k}"] = {"independent_error": err,
                                         "P_trace": float(expected[0]),
                                         "Q_trace": float(expected[1])}
    return {"a": float(a),"b": float(b),"collision_center": float(center),
            "predicted_k3_limit": float(predicted), "asymptotic": asymptotic,
            "coordinate_linearization": direct_rows}


if __name__ == "__main__":
    threshold = brentq(lambda h: h*math.tanh(h/2)-2, 2.0,3.0)
    print(json.dumps({"platform":platform.platform(),"sharp_speed_times_log_gap":threshold,
                      "collisions":[collision(2.0,1.5),collision(3.0,0.95)],
                      "fourier":fourier_checks()},indent=2,allow_nan=False))
