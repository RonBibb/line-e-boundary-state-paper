"""Check mean-metric identities and the exact asymmetric rank-two collision."""
import json
import math
import platform

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import j0,y0

from check_production_modulus import from_initial_data
from explore_modulus_collision import roots,candidate_ratio


def scalar_virial():
    modes = from_initial_data(None)
    q = 0.4**2+0.8**2+0.5*sum(mode.limiting_velocity**2 for mode in modes)
    intercept = [mode.regular+(2*mode.singular/math.pi)*(
                 math.log(mode.k/2)+np.euler_gamma) for mode in modes]
    product = 0.5*sum(p*mode.limiting_velocity for p,mode in zip(intercept,modes))
    initial_mean, initial_energy = -0.36, 1.17
    virial = initial_mean+initial_energy-q+product

    def energy_balance(tau):
        return 0.4**2+0.8**2-q+0.5*sum(
            mode.velocity(tau)**2+(mode.k*math.exp(-tau)*mode.value(tau))**2
            for mode in modes)

    direct = initial_mean+quad(energy_balance,0.,28.,epsabs=1e-13)[0]
    assert abs(direct-virial)<1e-12
    return {"initial_mean":initial_mean,"initial_energy":initial_energy,
            "limiting_mean_slope":q,"mean_times_intercept":product,
            "mean_modulus_from_exact_boundary_identity":virial,
            "mean_modulus_from_independent_quadrature":direct}


def asymmetric_pair():
    log_R = brentq(lambda x:candidate_ratio(math.exp(x))-6,math.log(64),math.log(1e10),xtol=1e-13)
    R = math.exp(log_R)
    b = roots(R)[0]
    a = b/R
    phi = math.pi*math.sqrt(2)/2
    eta, amplitude = 0.035, 0.025
    g,q = 0.55,0.37
    sigma = math.hypot(g,q)
    D = math.log((1+g/sigma)/(1-g/sigma))

    def profile(x):
        return math.cos(x)+eta*math.cos(6*x+phi)

    zeros = [brentq(profile,x-.15,x+.15)
             for x in (math.pi/2,3*math.pi/2)]
    positive_length = 2*math.pi-zeros[1]+zeros[0]
    sign_mean = positive_length/math.pi-1
    exact_mean_difference = -2*D*sign_mean

    def phase(tau,x):
        return -amplitude*sum(coef*math.cos(k*x+shift)*radial(k,tau)/j0(k*a)
                              for k,coef,shift in ((1,1.,0.),(6,eta,phi)))

    def radial(k,tau):
        # Asymptotic Bessel evaluation used only in a regime where the
        # first omitted term is below double roundoff.
        if tau > 28:
            return j0(k*a)*(math.log(k/2)+np.euler_gamma-tau)-math.pi/2*y0(k*a)
        z = k*math.exp(-tau)
        return math.pi/2*(j0(k*a)*y0(z)-y0(k*a)*j0(z))

    def metric_difference(tau,x):
        w = phase(tau,x)
        first = np.logaddexp(math.log((1+g/sigma)/2)+sigma*w,
                             math.log((1-g/sigma)/2)-sigma*w)
        second = np.logaddexp(math.log((1-g/sigma)/2)+sigma*w,
                              math.log((1+g/sigma)/2)-sigma*w)
        return float(-2*(first-second))

    means = {}
    for tau in (32.,128.,512.,2048.,8192.,32768.):
        breaks = sorted([0.,2*math.pi]+[root for root in zeros]+
                        [root+step/(amplitude*tau) for root in zeros for step in (-2,2)
                         if 0<root+step/(amplitude*tau)<2*math.pi])
        integral = sum(quad(lambda x:metric_difference(tau,x),left,right,epsabs=1e-11)[0]
                       for left,right in zip(breaks,breaks[1:]))
        means[str(tau)] = integral/(2*math.pi)
    assert abs(means["32768.0"]-exact_mean_difference)<1e-4
    return {"a":a,"b":b,"ratio":R,"even_harmonic_relative_size":eta,
            "zeros_of_limiting_seed_velocity":zeros,
            "positive_sign_fraction":positive_length/(2*math.pi),
            "exact_limiting_mean_metric_difference":exact_mean_difference,
            "observed_mean_metric_differences":means,
            "common_velocity_range_lower_bound":1-sigma*amplitude*(1+eta),
            "common_velocity_range_upper_bound":1.}


if __name__=="__main__":
    print(json.dumps({"platform":platform.platform(),"scalar_virial":scalar_virial(),
                      "asymmetric_pair":asymmetric_pair()},indent=2,allow_nan=False))
