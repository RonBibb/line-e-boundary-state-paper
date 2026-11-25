"""Check the two observation sidebands when neither velocity is supplied."""
import json
import math
import platform

import numpy as np
from scipy.special import jv,yv

from check_nonlinear_sampling import Collision,theta_grid,series
from check_null_transport import tangent_at,angle_formula


def scalar_gain(a,c,k,order):
    return math.pi*k*math.sqrt(a*c)/2*(
        jv(order,k*a)*yv(order,k*c)-yv(order,k*a)*jv(order,k*c))


def check():
    base = Collision(a=1.6,mode=4,amplitude=0.12)
    delays = (math.pi/12,math.pi/28)
    L = math.lcm(12,28)
    results = {}
    for probe in (2*L,4*L,8*L):
        n = 4*probe+144
        angular = theta_grid(n)
        seeds = np.zeros((2,2,n))
        seeds[0,0] = probe*np.cos(probe*angular)
        seeds[1,1] = probe*np.cos(probe*angular)*math.exp(-base.exact_fields(
            base.a,np.array([0.]),1.)[0,0])
        response = tangent_at(base,delays,seeds,7e-11,include_Q_velocity=True)
        projected_norms = []
        for initial_direction,direction in ((0,1),(1,0)):
            projected = []
            predicted = []
            for sign in (-1,1):
                side = probe+sign*base.mode
                phases = np.array([scalar_gain(base.a,base.a+d,side,0. if direction==0 else 1.)
                                   for d in delays])
                observed = []
                expected = []
                for slot,delay in enumerate(delays):
                    t = base.a+delay
                    components = response[slot,[0,2],initial_direction].copy()
                    components[1] *= np.exp(base.exact_fields(t,angular,1.)[0])
                    observed.append(side/probe*math.sqrt(t/base.a)*series(components[direction])[side])
                    phase,_ = base.radial_phase(t)
                    sgn = 1 if initial_direction==0 else -1
                    expected.append(sgn*base.amplitude*base.q*phase*math.cos(probe*delay)
                                    *side/(4j*probe))
                observed = np.array(observed)
                expected = np.array(expected)
                projected.append(float(np.linalg.norm(observed-phases*np.dot(phases,observed)/np.dot(phases,phases))))
                predicted.append(float(np.linalg.norm(expected-phases*np.dot(phases,expected)/np.dot(phases,phases))))
            projected_norms.append({"input_frame_component":initial_direction,
                                    "observed_projected_sidebands":projected,
                                    "predicted_projected_sidebands":predicted})
        results[str(probe)] = projected_norms
    angular = theta_grid(1024)
    angles = []
    for epsilon in (0.03,0.06,0.12,0.24):
        model = Collision(a=base.a,mode=base.mode,amplitude=epsilon)
        L_one = 16
        delay = math.pi/L_one
        angles_per_x = -2*sum(angle_formula(model,model.a+delay,angular-2*j*delay)
                              for j in range(L_one))
        phase,_ = model.radial_phase(model.a+delay)
        angles.append({"amplitude":epsilon,
                       "mean_angle_divided_by_amplitude_squared":float(np.mean(angles_per_x))/epsilon**2,
                       "max_angular_variation":float(np.ptp(angles_per_x)),
                       "predicted_ratio":model.q*model.g*L_one*phase**2/2})
    return {"period_lcm":L,"two_delays":list(delays),"two_time_sidebands":results,
            "single_time_holonomy":angles}


if __name__=="__main__":
    full = {"platform":platform.platform(),**check()}
    print(json.dumps(full,indent=2,allow_nan=False))
