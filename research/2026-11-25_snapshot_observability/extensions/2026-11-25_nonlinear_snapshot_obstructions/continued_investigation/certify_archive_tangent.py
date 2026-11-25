"""Interval certificate for the centered-grid Gowdy observation derivative.

The stored depth-18 floats are the exact initial data of this certificate.
Interval arithmetic encloses the analytic bootstrap and comparison system;
there is no trajectory sampling in the stated bounds.
"""
import hashlib
import json
from pathlib import Path
import math
import platform

import mpmath as mp
import numpy as np

mp.iv.dps = 55
iv = mp.iv
I = iv.mpf

INPUTS = Path(__file__).resolve().parents[5] / "reproduction" / "inputs"
DURATION = I(2)
CAP_V_MIN, CAP_V_MAX = I("0.31"), I("1.425")
CAP_B, CAP_C, CAP_Z = I("0.035"), I("0.276"), I("0.965")
CAP_INITIAL_JUMP = I("0.18")


def rightmost(values):
    return max(value.b for value in values)


def leftmost(values):
    return min(value.a for value in values)


def upward(value):
    return math.nextafter(float(value.b), math.inf)


def downward(value):
    return math.nextafter(float(value.a), -math.inf)


def vmd_bound(V, C):
    u = 1-V
    speed = iv.sqrt(u**2+C**2)
    if abs(C).b == 0:
        return I(0), speed
    if u.b <= 0:
        return abs(C), speed
    pos = iv.exp(DURATION*speed)
    neg = 1/pos
    ch, sh = (pos+neg)/2, (pos-neg)/2
    end_u_numerator = u*ch-speed*sh
    global_growth = abs(C)*iv.exp(DURATION*max(I(0),u.b))
    if end_u_numerator.a >= 0:
        denominator = ch-u/speed*sh
        if denominator.a > 0:
            return min((abs(C)/denominator).b,global_growth.b), speed
    # The conserved speed and elementary growth estimate are
    # both valid even when the turning time is uncertain.
    return min(speed.b,global_growth.b), speed


def velocity_potential_reference(V, C):
    u = 1-V
    speed = iv.sqrt(u**2+C**2)
    if speed.b == 0:
        return I(1)
    pos = iv.exp(DURATION*speed)
    neg = 1/pos
    return 1+iv.ln((pos+neg-u/speed*(pos-neg))/2)/DURATION


def validate_initial_data(path, grid, velocity_radius, measured_uncertainty=None):
    with np.load(path) as archive:
        P, V, Q, R = [[I(float(value)) for value in archive[f"tau18.0_{i}"]]
                      for i in range(4)]
        previous_P = [I(float(value)) for value in archive["tau16.0_0"]]
    assert all(len(data)==grid for data in (P,V,Q,R,previous_P))
    center_V = V
    center_P, center_Q = P,Q
    secant = [(start-end)/DURATION for start,end in zip(P,previous_P)]
    if velocity_radius.b > 0:
        V = [I([(v-velocity_radius).a,(v+velocity_radius).b]) for v in V]
    if measured_uncertainty is not None:
        p_noise,q_noise,pi_noise = measured_uncertainty
        P = [I([(p-p_noise).a,(p+p_noise).b]) for p in center_P]
        orig_Q = Q
        Q = [I([(q-q_noise*iv.exp(-p)).a,(q+q_noise*iv.exp(-p)).b])
             for p,q in zip(center_P,orig_Q)]
        momenta = [iv.exp(2*p)*r for p,r in zip(center_P,R)]
        R = [I([(momentum-pi_noise*iv.exp(p)).a,(momentum+pi_noise*iv.exp(p)).b])*iv.exp(-2*perturbed_p)
             for p,perturbed_p,momentum in zip(center_P,P,momenta)]
    delta = (2*iv.pi)/grid
    B, C, Z, jumps, center_C = [], [], [], [], []
    for i in range(grid):
        after, before = (i+1)%grid, (i-1)%grid
        scaled = iv.exp(P[i]-18)/delta
        B.append(iv.exp(P[i])*R[i])
        C.append(scaled*(Q[after]-Q[before])/2)
        Z.append(scaled*(Q[after]-Q[i]))
        jumps.append(abs(P[after]-P[i]))
        center_C.append(iv.exp(center_P[i]-18)/delta*(center_Q[after]-center_Q[before])/2)
    vmd_sites = [vmd_bound(v,c) for v,c in zip(V,C)]
    return {
        "minimum_V": leftmost(V), "maximum_V": rightmost(V),
        "B": rightmost(map(abs,B)), "C": rightmost(map(abs,C)),
        "Z": rightmost(map(abs,Z)), "P_jump": rightmost(jumps),
        "vmd_maximum_C": rightmost(row[0] for row in vmd_sites),
        "vmd_maximum_speed": rightmost(row[1] for row in vmd_sites),
        "delta": delta,
        "secant_center_distance": rightmost(abs(y-v) for y,v in zip(secant,center_V)),
        "secant_vmd_distance": rightmost(
            abs(y-velocity_potential_reference(v,c)) for y,v,c in zip(secant,center_V,center_C)),

    }


def bootstrap(initial):
    assert initial["P_jump"] < CAP_INITIAL_JUMP, {k: float(v.a) for k,v in initial.items()}
    xi = iv.exp(I(-16))/initial["delta"]
    jump = CAP_INITIAL_JUMP+DURATION*(CAP_V_MAX-CAP_V_MIN)
    exponential_jump = iv.exp(jump)
    spatial_P = 2*xi**2*jump

    upper_V = initial["maximum_V"]+DURATION*(CAP_C**2+spatial_P)
    lower_V = initial["minimum_V"]-DURATION*(CAP_B**2+spatial_P)
    forcing_B = xi*((1+exponential_jump)*CAP_Z+2*jump*CAP_C)
    e_V = iv.exp(DURATION*CAP_V_MAX)
    upper_B = e_V*initial["B"]+forcing_B*(e_V-1)/CAP_V_MAX

    forcing_Z = xi*(1+exponential_jump)*CAP_B
    Z_growth = 1-CAP_V_MIN
    e_Z = iv.exp(DURATION*Z_growth)
    upper_Z = e_Z*initial["Z"]+forcing_Z*(e_Z-1)/Z_growth

    # u=1-V and C obey the VMD system plus the two stated forcing terms.
    # Bound its Jacobian along the segment between the VMD and full states.
    abs_u = max((1-CAP_V_MIN).b,(CAP_V_MAX-1).b,initial["vmd_maximum_speed"])
    abs_C = max(CAP_C.b,initial["vmd_maximum_C"])
    comparison_L = max((2*abs_C).b,(abs_u+abs_C).b)
    forcing_C = xi*exponential_jump*CAP_B
    forcing = max((CAP_B**2+spatial_P).b,forcing_C.b)
    vmd_error = forcing*(iv.exp(DURATION*comparison_L)-1)/comparison_L
    upper_C = initial["vmd_maximum_C"]+vmd_error

    assert lower_V.a > CAP_V_MIN and upper_V.b < CAP_V_MAX, (lower_V,upper_V)
    secant_center_mismatch = initial["secant_vmd_distance"]+forcing*(
        (iv.exp(DURATION*comparison_L)-1)/comparison_L-DURATION
    )/(DURATION*comparison_L)

    assert upper_B.b < CAP_B and upper_C.b < CAP_C and upper_Z.b < CAP_Z, (
        upper_B,upper_C,upper_Z)
    return {"xi":xi, "jump":jump, "exp_jump":exponential_jump, "spatial_P":spatial_P,
            "lower_V":lower_V, "upper_V":upper_V, "upper_B":upper_B,
            "upper_C":upper_C, "upper_Z":upper_Z, "vmd_error":vmd_error,
            "center_observation_mismatch_bound":secant_center_mismatch}


def comparison_matrix(bound):
    xi, E, jump = (bound[key] for key in ("xi","exp_jump","jump"))
    C, B, V = (bound[key].b for key in ("upper_C","upper_B","upper_V"))
    a = 2*max(C**2,B**2)+4*xi**2
    b = 2*B
    c = 2*xi*C*E
    z_p = 2*xi*C
    z_y = (2+2*E*(1+jump))*xi**2

    return [[I(0),I(1),I(0),I(0)], [a,I(0),c,b],
            [I(0),I(0),I(0),I(1)], [z_p,b,z_y,V]]


def comparison_flow(M,start):
    term = start[:]
    total = start[:]
    order = 24
    for k in range(1,order+1):
        term = [DURATION*sum(m*v for m,v in zip(row,term))/k for row in M]
        total = [x+y for x,y in zip(total,term)]
    scaled_norm = DURATION*rightmost(sum(row) for row in M)
    tail = rightmost(start)*scaled_norm**(order+1)/math.factorial(order+1)/(1-scaled_norm/(order+2))
    return total[0]+tail,tail,scaled_norm


def tangent_bound(bound):
    terminal,tail,scaled_norm = comparison_flow(comparison_matrix(bound),[I(0),I(1),I(0),I(0)])
    defect = terminal/DURATION-1
    assert defect.b < I("0.3"), defect
    return {"matrix_row_norm": upward(scaled_norm/DURATION),
            "exponential_series_tail": upward(tail),
            "observation_defect_upper_bound": upward(defect),
            "inverse_norm_upper_bound": upward(1/(1-defect))}


def certificate(n, radius):
    path = INPUTS/"bb-cbe3b"/f"cbe3b_slices_N{n}.npz"
    initial = validate_initial_data(path,n,radius)
    enclosure = bootstrap(initial)
    tangent = tangent_bound(enclosure)
    defect = I(tangent["observation_defect_upper_bound"])
    self_map_radius = (1-defect)*radius-enclosure["center_observation_mismatch_bound"]
    assert initial["secant_center_distance"] < radius
    assert self_map_radius.a > 0
    initial_report = {key:(downward(value) if key=="minimum_V" else upward(value))
                      for key,value in initial.items()}
    enclosure_report = {key:(downward(value) if key=="lower_V" else upward(value))
                        for key,value in enclosure.items()}
    return {"initial_velocity_radius": "0.06",
            "certified_additional_observation_error_radius": downward(self_map_radius),
            "velocity_distance_to_hidden_center_bound": upward(
                enclosure["center_observation_mismatch_bound"]/(1-defect)),
            "source_file": path.name, "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "initial_bounds": initial_report, "bootstrap_bounds": enclosure_report,
            "tangent_certificate":tangent}


if __name__ == "__main__":
    results = {str(n):certificate(n,I("0.06")) for n in (1024,2048,4096)}
    print(json.dumps({"platform":platform.platform(),"mpmath":mp.__version__,
                      "interval_digits":iv.dps,"results":results},indent=2,allow_nan=False))
