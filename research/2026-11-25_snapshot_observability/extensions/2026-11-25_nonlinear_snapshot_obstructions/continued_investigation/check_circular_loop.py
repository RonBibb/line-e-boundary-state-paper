"""Nonlinear circular-loop nodes, detuning, and small subcritical speed."""
import json
import math
import platform

import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.optimize import brentq
from scipy.special import j0,y0


RATIO = math.e
def linear_radius(nu,r):
    return math.pi/2*(j0(nu)*y0(nu*r)-y0(nu)*j0(nu*r))


def nonlinear_placement(nu,epsilon,tolerance):
    def rhs(r,state):
        radius,derivative = state
        return derivative,-derivative/r-nu**2*math.sinh(radius)*math.cosh(radius)
    result = solve_ivp(rhs,(1.,RATIO),(0.,epsilon),method="DOP853",
                        rtol=tolerance,atol=tolerance*epsilon/20)
    assert result.success
    return float(result.y[0,-1])


def asymptotic_speed(nu,epsilon,tolerance):
    def rhs(x,state):
        radius,velocity = state
        return velocity,-nu**2*math.exp(-2*x)*math.sinh(radius)*math.cosh(radius)
    result = solve_ivp(rhs,(0.,42.),(0.,-epsilon),method="DOP853",
                        rtol=tolerance,atol=tolerance*epsilon/20)
    assert result.success
    return float(abs(result.y[1,-1]))


def experiment(tolerance):
    nu0 = brentq(lambda nu:linear_radius(nu,RATIO),1.5,2.2,xtol=4e-15)
    moments = [quad(lambda r:r*linear_radius(nu0,r)**exponent,1.,RATIO,epsabs=1e-12)[0]
               for exponent in (2,4)]
    cubic_shift = -(nu0/3)*moments[1]/moments[0]
    rows = {}
    for epsilon in (0.035,0.07,0.14,0.28,0.45):
        nu = brentq(lambda trial:nonlinear_placement(trial,epsilon,tolerance),
                     1.35,2.22,xtol=5e-15)
        a = nu/4
        theta = 2*np.pi*np.arange(512)/512
        r = nonlinear_placement(nu,epsilon,tolerance)
        circular_plus = np.log(np.cosh(r)+np.sinh(r)*np.cos(4*theta))
        circular_minus = np.log(np.cosh(r)-np.sinh(r)*np.cos(4*theta))
        speed = asymptotic_speed(nu,epsilon,tolerance)
        assert 0 < speed < epsilon < 1
        rows[str(epsilon)] = {
            "nu":nu,"a":a,"b":a*RATIO,
            "terminal_orbit_P_max_difference":float(np.max(abs(circular_plus-circular_minus))),
            "placement_shift_over_epsilon_squared":(nu-nu0)/epsilon**2,
            "unretimed_terminal_radius_over_epsilon_cubed":
                 nonlinear_placement(nu0,epsilon,tolerance)/epsilon**3,
            "limiting_speed":speed,"speed_over_epsilon":speed/epsilon}
    return {"linear_placement":nu0,"predicted_placement_shift":cubic_shift,
            "predicted_speed_over_epsilon":float(abs(j0(nu0))),"amplitudes":rows}


if __name__=="__main__":
    print(json.dumps({"platform":platform.platform(),"ratio":RATIO,
                      "checks":{str(tol):experiment(tol) for tol in (2e-11,5e-13)}},
                     indent=2,allow_nan=False))
