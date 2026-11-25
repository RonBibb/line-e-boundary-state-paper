"""Independent coordinate linearization checks of null transport and snapshot gains."""
from dataclasses import replace
import json
import math
import platform

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from check_nonlinear_sampling import Collision, series, theta_grid


def tuned_collision():
    mode = 4
    gap = 15*math.pi/(16*mode)
    seed = Collision(mode=mode, amplitude=0.24)
    find_a = lambda a: replace(seed, a=a).radial_phase(a+gap)[0]
    a = brentq(find_a, 1e-7, 1., xtol=1e-15)
    return replace(seed,a=a), gap


def spectral_derivatives(values, modes, n):
    transformed = np.fft.rfft(values,axis=-1)
    first = np.fft.irfft(1j*modes*transformed,n=n,axis=-1)
    second = np.fft.irfft(-modes**2*transformed,n=n,axis=-1)
    return first, second


def tangent_at(collision, delays, seeds, tolerance, include_Q_velocity=False):
    n = seeds.shape[-1]
    theta = theta_grid(n)
    modes = np.arange(n//2+1)
    if include_Q_velocity:
        starting = np.zeros((4,*seeds[0].shape))
        starting[1],starting[3] = seeds
    else:
        starting = np.zeros((4,*seeds.shape))
        starting[1] = seeds

    def rhs(t, flattened):
        P, Pt, Q, Qt = collision.exact_fields(t,theta,1.)
        (Px,Qx), _ = spectral_derivatives(np.stack([P,Q]),modes,n)
        p, pt, q, qt = flattened.reshape(starting.shape)
        (px,qx), (pxx,qxx) = spectral_derivatives(np.stack([p,q]),modes,n)
        nonlinearity_p = 2*np.exp(2*P)*((Qt**2-Qx**2)*p+Qt*qt-Qx*qx)
        nonlinearity_q = -2*(Qt*pt+Pt*qt-Qx*px-Px*qx)
        return np.stack([pt,pxx-pt/t+nonlinearity_p,qt,qxx-qt/t+nonlinearity_q]).ravel()

    times, order = np.unique(collision.a+np.array(delays),return_inverse=True)
    result = solve_ivp(rhs,(collision.a,times[-1]),starting.ravel(),method="DOP853",
                       rtol=tolerance,atol=tolerance*0.01,t_eval=times)
    if not result.success:
        raise RuntimeError(result.message)
    return result.y.T[order].reshape((len(delays),*starting.shape))


def seed_function(collision, w):
    sigma = math.hypot(collision.g,collision.q)
    return collision.P0+np.log(np.cosh(sigma*w)+collision.g*np.sinh(sigma*w)/sigma)


def angle_formula(collision, t, theta):
    u, _ = collision.radial_phase(t)
    w = collision.amplitude*u*np.cos(collision.mode*theta)
    nodes, weights = np.polynomial.legendre.leggauss(12)
    points = w[:,None]*(nodes+1)/2
    return collision.current*w/2*np.sum(weights*np.exp(-seed_function(collision,points)),axis=-1)


def independent_ray_integral(collision, delay, x, sign):
    nodes, weights = np.polynomial.legendre.leggauss(192)
    product = np.zeros_like(x)
    for parameter, weight in zip(delay*(nodes+1)/2, weights):
        time = collision.a+parameter
        ray_x = x+sign*(parameter-delay)
        _, ut = collision.radial_phase(time)
        base = collision.exact_fields(time,ray_x,1.)
        Qtheta = collision.current*time*collision.amplitude*ut*np.cos(collision.mode*ray_x)
        product += weight*np.exp(base[0])*(base[3]+sign*Qtheta)
    return -delay*product/2


def l2(fields):
    return float(np.sqrt(np.sum(np.mean(fields**2,axis=-1))))


def high_mode(collision, nodal_delay, delay, k, n, tolerance):
    theta = theta_grid(n)
    tangent = tangent_at(collision,(delay,nodal_delay),(k*np.cos(k*theta))[None,:],tolerance)
    results = []
    for offset, tangent_fields in zip((delay,nodal_delay),tangent):
        t = collision.a+offset
        rho = angle_formula(collision,t,theta)
        incoming, outgoing = np.sin(k*(theta+offset)),np.sin(k*(theta-offset))
        principal = 0.5*np.stack([np.cos(rho)*(incoming-outgoing),
                                  np.sin(rho)*(incoming+outgoing)])
        base = collision.exact_fields(t,theta,1.)
        exact = math.sqrt(t/collision.a)*np.stack([tangent_fields[0,0],
                                                  np.exp(base[0])*tangent_fields[2,0]])
        results.append({"delay": offset,
                        "full_gain_for_sine_primitive": l2(exact)/l2(np.sin(k*theta)),
                        "principal_gain_for_sine_primitive": l2(principal)/l2(np.sin(k*theta)),
                        "position_channel_L2": l2(exact[0]),
                        "transverse_channel_L2": l2(exact[1]),
                        "principal_error_L2": l2(exact-principal),
                        "k_times_principal_error_L2": k*l2(exact-principal)})
    return results


def weighted_output(values, exponent, cutoff):
    coeffs = series(values)[...,:cutoff+1]
    ks = np.arange(cutoff+1)
    weights = np.power(1+ks**2,exponent/2)
    return np.concatenate([coeffs[...,:1].real,
                           math.sqrt(2)*weights[1:]*coeffs[...,1:].real,
                           math.sqrt(2)*weights[1:]*coeffs[...,1:].imag],axis=-1)


def gains_by_configuration(collision, delays, cutoff, n, tolerance):
    theta = theta_grid(n)
    nonzero = np.arange(1,cutoff+1)
    phases = nonzero[:,None]*theta[None,:]
    scales = np.sqrt(2/(1+nonzero**2))[:, None]
    seeds = np.concatenate([np.ones((1,n)), scales*np.cos(phases), scales*np.sin(phases)])
    branches = [[] for _ in delays]
    for start in range(0,len(seeds),12):
        evolving = tangent_at(collision,delays,seeds[start:start+12],tolerance)
        for i,(d,values) in enumerate(zip(delays,evolving)):
            t = collision.a+d
            P = collision.exact_fields(t,theta,1.)[0]
            positions = weighted_output(values[0],2,min(n//3,cutoff+2*collision.mode))
            transverse = weighted_output(np.exp(P)*values[2],2,min(n//3,cutoff+2*collision.mode))
            branches[i].append(np.concatenate([positions,transverse],axis=-1))
    nodal, terminal, helper = [np.concatenate(block,axis=0).T for block in branches]
    P_size = helper.shape[0]//2
    layouts = {"nodal_full": nodal, "separated_full": terminal,
               "separated_position_only": terminal[:P_size],
               "separated_full_plus_earlier_position": np.concatenate([terminal,helper[:P_size]],axis=0),
               "separated_full_plus_nodal_position": np.concatenate([terminal,nodal[:P_size]],axis=0)}
    results = {}
    for name, matrix in layouts.items():
        values = np.linalg.svd(matrix,compute_uv=False)
        results[name] = {"minimum": float(values[-1]), "maximum": float(values[0]),
                         "smallest_five": values[-5:].tolist()}
    return results


if __name__ == "__main__":
    collision, b = tuned_collision()
    c = math.pi/16
    d1 = 0.0131
    angular = theta_grid(1024)
    rho = angle_formula(collision,collision.a+c,angular)
    ray_error = max(
        float(np.max(np.abs(independent_ray_integral(collision,c,angular,sign)+sign*rho)))
        for sign in (1.,-1.))
    cover = float(np.min(np.sin(rho)**2+
                         np.sin(angle_formula(collision,collision.a+c,angular-2*c))**2))
    assert ray_error < 1e-10 and cover > 0, (collision.a, ray_error, cover)
    modes = {}
    for k in (64,128,256,512,1024):
        n = max(384, 4*k+128)
        modes[str(k)] = high_mode(collision,b,c,k,n,2e-10)
    matrices = {}
    for K,n in ((24,256),(48,384),(96,512),(96,768)):
        matrices[f"{K}_{n}"] = gains_by_configuration(collision,(b,c,d1),K,n,8e-11)
    print(json.dumps({"platform": platform.platform(), "a": collision.a,
                      "delays":{"nodal":b,"separated":c,"helper":d1},
                      "seed_mode":collision.mode,"amplitude":collision.amplitude,
                      "maximum_ray_integral_disagreement":ray_error,
                      "minimum_sum_squared_cover_sample":cover,
                      "high_mode_results":modes, "finite_band_singular_values":matrices},
                     indent=2,allow_nan=False))
