"""Independently compute the limiting modulus of the historical smooth holdout."""
from dataclasses import dataclass
import importlib.util
import json
import math
from pathlib import Path
import platform

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import j0, j1, y0, y1

INPUTS = Path(__file__).resolve().parents[5] / "reproduction" / "inputs"
SOURCE = INPUTS/"bb-cbe2i"/"run_cbe2i.py"
MODES = (("P",1,"cos",0.6,0.), ("P",2,"cos",0.,0.3),
         ("phi",1,"sin",0.5,0.), ("phi",2,"cos",0.,-0.2))
DEPTHS = (8,12,16,24)


@dataclass(frozen=True)
class BesselMode:
    field: str
    k: int
    parity: str
    radial_frequency: float
    spatial_frequency: float
    regular: float
    singular: float

    def value(self, tau):
        arg = self.radial_frequency*math.exp(-tau)
        return self.regular*j0(arg)+self.singular*y0(arg)

    def velocity(self, tau):
        arg = self.radial_frequency*math.exp(-tau)
        return arg*(self.regular*j1(arg)+self.singular*y1(arg))

    @property
    def limiting_velocity(self):
        return -2*self.singular/math.pi


def mode_frequencies(k, mesh_size):
    if mesh_size is None:
        return float(k), float(k)
    dx = 2*math.pi/mesh_size
    return 2*math.sin(k*dx/2)/dx, math.sin(k*dx)/dx


def from_initial_data(mesh_size):
    modes = []
    for field,k,parity,value,velocity in MODES:
        radial, spatial = mode_frequencies(k,mesh_size)
        matrix = [[j0(radial),y0(radial)],[radial*j1(radial),radial*y1(radial)]]
        regular, singular = np.linalg.solve(matrix,[value,velocity])
        modes.append(BesselMode(field,k,parity,radial,spatial,float(regular),float(singular)))
    return modes


def from_sample_pair(first, second, start, stop, mesh_size, use_stencil):
    modes = []
    means = tuple(float(np.mean(second[i]-first[i])/(stop-start)) for i in (0,2))
    series = {name:[2*np.fft.rfft(tape[i])/mesh_size for tape in (first,second)]
              for name,i in (("P",0),("phi",2))}
    for field,k,parity,_,_ in MODES:
        radial, spatial = mode_frequencies(k,mesh_size if use_stencil else None)
        left, right = (radial*math.exp(-tau) for tau in (start,stop))
        matrix = [[j0(left),y0(left)],[j0(right),y0(right)]]
        values = [x[k].real if parity=="cos" else -x[k].imag for x in series[field]]
        regular, singular = np.linalg.solve(matrix,values)
        modes.append(BesselMode(field,k,parity,radial,spatial,float(regular),float(singular)))
    return modes,means


def modulus(modes, means, tau, anchor):
    P1, P2, F1, F2 = modes
    vp = means[0]+P1.limiting_velocity+P2.limiting_velocity
    vf = means[1]+F2.limiting_velocity
    q = vp*vp+vf*vf

    def remainder(time):
        p = means[0]+P1.velocity(time)+P2.velocity(time)
        f = means[1]+F2.velocity(time)
        weighted_gradient = math.exp(-time)*F1.spatial_frequency*F1.value(time)
        return (p-vp)*(p+vp)+(f-vf)*(f+vf)+weighted_gradient**2

    stop = max(tau+8,32.)
    correction, quadrature_error = quad(remainder,tau,stop,epsabs=3e-14)
    return {"modulus":anchor-tau*q+correction,"q_at_origin":q,
            "tail_from_anchor_used_in_estimator":correction,
            "double_quadrature_reported_error":quadrature_error}


def reference_high_precision(decimal_digits):
    with mp.workdps(decimal_digits):
        k1, k2 = mp.mpf(1), mp.mpf(2)
        j = mp.besselj
        y = mp.bessely
        coefficients = []
        for (_,k,_,value,derivative),k_mp in zip(MODES,(k1,k2,k1,k2)):
            matrix = mp.matrix([[j(0,k_mp),y(0,k_mp)],[k_mp*j(1,k_mp),k_mp*y(1,k_mp)]])
            A,B = mp.lu_solve(matrix,mp.matrix([mp.mpf(str(value)),mp.mpf(str(derivative))]))
            coefficients.append((k_mp,A,B))

        def values(mode,tau):
            k,A,B = mode
            z = k*mp.exp(-tau)
            return A*j(0,z)+B*y(0,z),z*(A*j(1,z)+B*y(1,z))

        asymptotic = [-2*B/mp.pi for _,_,B in coefficients]
        pv = mp.mpf("0.4")+asymptotic[0]+asymptotic[1]
        fv = mp.mpf("0.8")+asymptotic[3]

        def integrand(time):
            p = mp.mpf("0.4")+values(coefficients[0],time)[1]+values(coefficients[1],time)[1]
            f = mp.mpf("0.8")+values(coefficients[3],time)[1]
            phi_mode = values(coefficients[2],time)[0]
            return (p-pv)*(p+pv)+(f-fv)*(f+fv)+mp.exp(-2*time)*phi_mode**2

        T = mp.mpf(34)
        integral = mp.quad(integrand,[0,2,6,12,20,T])
        # For tau=T+x, field modes are bounded by M*(1+x) using
        # a Volterra estimate of the Bessel equation in tau.
        decayed = mp.exp(-2*T)
        M = [max(map(abs,values(mode,T)))/(1-3*mode[0]**2*decayed/4)
             for mode in coefficients]
        p_bound = (M[0]+4*M[1])*decayed*mp.mpf("0.75")
        f_bound = 4*M[3]*decayed*mp.mpf("0.75")
        tail_bound = (mp.mpf("0.75")*(2*abs(pv)*p_bound+2*abs(fv)*f_bound)
                      +mp.mpf("1.25")*(p_bound**2+f_bound**2+M[2]**2*decayed))
        return {"modulus":mp.nstr(integral,45),"cutoff_tau":float(T),
                "analytic_integral_tail_bound":mp.nstr(tail_bound,12),
                "q_at_origin":mp.nstr(pv**2+fv**2,45)}


def reload_model():
    spec = importlib.util.spec_from_file_location("historical_cbe2i",SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archived_estimates(model, n):
    saved = tuple(sorted({depth-offset for depth in DEPTHS for offset in (0,2,4)}))
    end,_,series = model.evolve(model.production_data(n),0.,float(max(DEPTHS)),
                               n,store_taus=saved)
    series[max(DEPTHS)] = end
    discrete_ref = modulus(from_initial_data(n),(0.4,0.8),0.,0.)
    data = {}
    for depth in DEPTHS:
        t0,t1 = depth-4,depth-2
        first,second = series[t0],series[t1]
        estimates = np.array([(second[i][0]-first[i][0])/2 for i in (0,2)])
        q = float(estimates@estimates)
        prior_anchor = float(second[4][0])
        end_anchor = float(series[depth][4][0])
        data[str(depth)] = {
            "historical_endpoint_constant":end_anchor-depth*q,
            "earlier_anchor_constant":prior_anchor-t1*q,
            "q_zero_remainder_estimate":q,
            "tail_corrected_matched_stencil":modulus(
                *from_sample_pair(first,second,t0,t1,n,True),t1,prior_anchor),
            "tail_corrected_continuum_bessel":modulus(
                *from_sample_pair(first,second,t0,t1,n,False),t1,prior_anchor)}
    return {"discrete_exact_reference":discrete_ref,"estimates":data}


if __name__=="__main__":
    comparison = [reference_high_precision(48),reference_high_precision(68)]
    assert abs(float(comparison[0]["modulus"])-float(comparison[1]["modulus"])) < 1e-13
    reference = float(comparison[-1]["modulus"])
    model = reload_model()
    results = {str(n):archived_estimates(model,n) for n in (256,512,1024,2048)}
    errors = {}
    for depth in DEPTHS:
        errors[str(depth)] = {
            key:{n:row["estimates"][str(depth)][key]-reference for n,row in results.items()}
            for key in ("historical_endpoint_constant","earlier_anchor_constant")}
    print(json.dumps({"platform":platform.platform(),"independent_continuum_reference":comparison,
                      "grids":results,"continuum_errors":errors},indent=2,allow_nan=False))
