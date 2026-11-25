"""Check physical sizes in the archived finite-window trace."""
from pathlib import Path
import json
import numpy as np

INPUTS = Path(__file__).resolve().parents[5] / "reproduction" / "inputs"
results = {}
for n in (1024,2048,4096):
    dx = 2*np.pi/n
    def d1(x):
        return (np.roll(x,-1)-np.roll(x,1))/(2*dx)
    def d2(x):
        return (np.roll(x,-1)-2*x+np.roll(x,1))/dx**2
    with np.load(INPUTS/"bb-cbe3b"/f"cbe3b_slices_N{n}.npz") as archive:
        rows = {}
        for t in (8.,10.,12.,14.,16.,18.):
            P,V,Q,R,lam = (archive[f"tau{t}_{i}"] for i in range(5))
            a = np.exp(-t)
            B = np.exp(P)*R
            C = np.exp(P)*a*d1(Q)
            residual = a*a*d2(P)+B*B-C*C
            rows[str(t)] = {"range_P":[float(P.min()),float(P.max())],
                           "range_V":[float(V.min()),float(V.max())],
                           "max_B":float(np.max(np.abs(B))),
                           "max_C":float(np.max(np.abs(C))),
                           "max_acceleration":float(np.max(np.abs(residual))),
                           "max_spatial_P_contribution":float(np.max(np.abs(a*a*d2(P)))),
                           "max_weighted_Q_second":float(np.max(np.abs(a*a*np.exp(P)*d2(Q))))}
        secant=(archive["tau18.0_0"]-archive["tau16.0_0"])/2
        rows["secant_distance_from_depth18_velocity"] = float(np.max(np.abs(secant-archive["tau18.0_1"])))
        results[str(n)] = rows
print(json.dumps(results,indent=2,allow_nan=False))
