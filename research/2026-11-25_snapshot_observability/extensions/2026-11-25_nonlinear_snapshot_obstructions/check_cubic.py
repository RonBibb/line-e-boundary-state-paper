"""Cloud-free periodic wave-map verification of the cubic terminal signal."""
import json
import math
import platform

import numpy as np
from scipy.integrate import solve_ivp

MODE = 2
GAP = math.pi / MODE
SPEED = 0.5
AMPLITUDES = (0.40, 0.20, 0.10, 0.05, 0.025)


def evolved(a, eps, n, rtol):
    theta = 2*math.pi*np.arange(n)/n
    wave_numbers = MODE*np.fft.rfftfreq(n,d=1/n)
    initial = np.zeros((5,n))
    initial[1] = eps*np.cos(theta)
    initial[3] = SPEED/a
    def rhs(t, flat):
        P, Pt, Q, Qt, _ = flat.reshape((5,n))
        fields = np.fft.rfft(np.stack([P,Q]),axis=-1)
        Px,Qx = np.fft.irfft(1j*wave_numbers*fields,n=n,axis=-1)
        Pxx,Qxx = np.fft.irfft(-wave_numbers**2*fields,n=n,axis=-1)
        e = np.exp(2*P)
        return np.stack([
            Pt,
            Pxx-Pt/t+e*(Qt**2-Qx**2),
            Qt,
            Qxx-Qt/t-2*(Pt*Qt-Px*Qx),
            -t*(Pt**2+Px**2+e*(Qt**2+Qx**2)),
        ]).ravel()
    result = solve_ivp(rhs,(a,a+GAP),initial.ravel(),method="DOP853",
                       rtol=rtol,atol=rtol*0.05)
    assert result.success
    return result.y[:,-1].reshape((5,n))


def cubic_coefficient(a):
    # Variations in geodesic Fermi coordinates L,N. A dot denotes d/dt.
    # N_1 = u cos(kx). The tangential coefficient L_2 has modes 0,2k.
    # The normal coefficient N_3 has modes k,3k. This ODE follows the
    # expansion in powers of eps without factorials.
    k=MODE
    def rhs(t,s):
        l0,l0t,l2,l2t,n1,n1t,n3,n3t = s
        delta=t-a
        lrate=SPEED/t
        u=math.sqrt(a/t)*math.sin(k*delta)/k
        ut=math.sqrt(a/t)*math.cos(k*delta)-u/(2*t)
        F=-u*ut*lrate
        fk=(lrate*lrate/2)*u**3+2*lrate*u*(l0t+0.5*l2t)
        f3=(lrate*lrate/6)*u**3+lrate*u*l2t
        wp=1/(4*t*t)
        return (l0t, -l0t/t+F, l2t, -l2t/t-(2*k)**2*l2+F,
                n1t,-n1t/t-(k*k-wp)*n1+fk,
                n3t,-n3t/t-((3*k)**2-wp)*n3+f3)
    result=solve_ivp(rhs,(a,a+GAP),np.zeros(8),method="DOP853",
                     rtol=1e-12,atol=1e-16)
    assert result.success
    response_P=2.0/math.cosh(SPEED*math.log((a+GAP)/a))
    return {"P_difference_cos1": response_P*float(result.y[4,-1]),
            "P_difference_cos3": response_P*float(result.y[6,-1]),
            "a2_N3_cos1": a*a*float(result.y[4,-1])}


def experiment(n,rtol):
    rows={}
    for eps in AMPLITUDES:
        plus=evolved(1.0,eps,n,rtol)
        minus=evolved(1.0,-eps,n,rtol)
        delta=(plus-minus)[[0,2,4]]
        p_modes=(2.0/n*np.fft.rfft(delta[0])).real
        sphere_pair=float(np.max(np.sqrt(np.sum(delta*delta,axis=0))))
        rows[str(eps)]={"P_difference_cos1_per_eps3":float(p_modes[1]/eps**3),
                        "P_difference_cos3_per_eps3":float(p_modes[3]/eps**3),
                        "P_Q_lambda_max_gap":sphere_pair}
    return rows


if __name__ == "__main__":
    result={"platform":platform.platform(),"a":1.0,"gap":GAP,"k":MODE,
            "coefficient_checks":{str(a):cubic_coefficient(a) for a in (1.,8.,32.,128.,512.)},
            "runs":{f"{n}_{tol}":experiment(n,tol)
                    for n,tol in ((48,1e-11),(96,1e-11),(192,1e-11),(96,2e-13))}}
    exact_large_a=-(9*SPEED**2*GAP)/(64*MODE**4)
    result["a2_N3_cos1_asymptotic"]=exact_large_a
    assert abs(result["coefficient_checks"]["512.0"]["a2_N3_cos1"]/exact_large_a-1)<0.02
    comparison=result["runs"]["96_2e-13"]
    expected=result["coefficient_checks"]["1.0"]["P_difference_cos1"]
    measured=comparison["0.025"]["P_difference_cos1_per_eps3"]
    assert abs(measured-expected)<2e-4*abs(expected),(expected,measured)
    print(json.dumps(result,indent=2,allow_nan=False))
