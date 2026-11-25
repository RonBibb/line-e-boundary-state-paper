"""Decode perturbed rank-two data using one additional full fixed slice."""
import json
import math
import platform

import numpy as np

from check_nonlinear_sampling import Collision,evolve,theta_grid,series,to_grid,h1_norm,spectral_operator
from check_null_transport import weighted_output
from check_open_two_snapshot import tangent_matrix


def basis(grid,K):
    theta = theta_grid(grid)
    k = np.arange(1,K+1)[:,None]
    phases = k*theta[None,:]
    scale = np.sqrt(2/(1+k*k))
    return np.concatenate((np.ones((1,grid)),scale*np.cos(phases),scale*np.sin(phases)))


def basis_coordinates(fourier):
    modes = np.arange(1,len(fourier))
    scale = np.sqrt(2*(1+modes*modes))
    return np.concatenate((np.array([fourier[0].real]),
                           scale*fourier[1:].real,-scale*fourier[1:].imag))


def metric_gradient(state,t):
    P,Pt,Q,Qt = state
    spatial,_ = spectral_operator(state.shape[-1])
    (Px,Qx),_ = spatial(state[[0,2]])
    return -2*t*(Pt*Px+np.exp(2*P)*Qt*Qx)


def fixed_map(samples, terminal_P, cutoff):
    p,q,lam_x = samples
    position = weighted_output(p,2,cutoff)
    orbit = weighted_output(np.exp(terminal_P)*q,2,cutoff)
    metric = weighted_output(lam_x,1,cutoff)
    return {"only_P":position,
            "full_orbit":np.concatenate([position,orbit]),
            "P_and_metric":np.concatenate([position,metric])}


def run(model,K,reference_n,inverse_n,noise):
    a,d = model.a,math.pi/16
    c = a+d
    tolerance = 5e-12
    reference = model.exact_fields(a,theta_grid(reference_n),1.)
    rng = np.random.default_rng(121)
    sites = np.array([1,3,16,32,48,K])
    sites = np.unique(sites[sites<=K])
    transfer = series(reference[1])[:K+1].copy()
    for k in sites:
        transfer[k] += (0.007/(1+k))*(rng.normal()+1j*rng.normal())
    reference[1] = to_grid(transfer,reference_n)
    reference_terminal = evolve(reference,a,np.array([d]),"fourier",tolerance)[0]
    p_data,q_data,lam_x_data = (series(field) for field in (*reference_terminal[[0,2]],metric_gradient(reference_terminal,c)))
    initial = model.exact_fields(a,theta_grid(inverse_n),0.)
    p_grid,q_grid,lam_x_grid = (to_grid(samples[:inverse_n//3+1],inverse_n)
                     for samples in (p_data,q_data,lam_x_data))

    trial_basis = basis(inverse_n,K)
    secant = (p_data[:K+1]-series(initial[0])[:K+1])/d
    guess = basis_coordinates(secant)
    initial[1] = guess@trial_basis
    operators,terminal_P,cutoff = tangent_matrix(initial,a,c,K,
                                                  return_operators=True)
    exact_record = fixed_map((p_grid,q_grid,lam_x_grid),terminal_P,cutoff)
    first_resonant = min(4*model.mode,K)
    truth = basis_coordinates(transfer)
    matrices = {name:operators[name] for name in exact_record}
    results = {}
    for arm,matrix in matrices.items():
        measurement_noise = np.zeros_like(exact_record[arm])
        measurement_noise[first_resonant] = noise
        inverse = np.linalg.pinv(matrix,rcond=1e-10)
        y = exact_record[arm]+measurement_noise
        u = guess.copy()
        history = []
        for iteration in range(7):
            initial[1] = u@trial_basis
            candidate = evolve(initial,a,np.array([d]),"fourier",tolerance)[0]
            computed = fixed_map((*candidate[[0,2]],metric_gradient(candidate,c)),terminal_P,cutoff)[arm]
            residual = y-computed
            history.append({"iteration":iteration,"record_residual":float(np.linalg.norm(residual)),
                            "velocity_H1_error":float(np.linalg.norm(u-truth))})
            step = inverse@residual
            if np.linalg.norm(step)>0.75 or np.linalg.norm(residual)<1e-10:
                break
            u += step
        singular = np.linalg.svd(matrix,compute_uv=False)
        results[arm] = {"smallest_gain_at_secant":float(singular[-1]),
                        "linear_noise_amplification_in_injected_direction":
                         float(np.linalg.norm(inverse@measurement_noise)/noise) if noise else None,
                        "iterations":history}
    return {"bandwidth":K,"reference_grid":reference_n,"inverse_grid":inverse_n,
            "noise_norm":noise,"delay":d,"initial_velocity_norm":h1_norm(transfer),
            "results":results}


if __name__=="__main__":
    pair = Collision(a=1.6,mode=4,amplitude=0.24)
    results = [run(pair,K,reference,grid,noise)
               for K,reference,grid,noise in
               ((48,768,384,0.),(96,1536,512,0.),(96,1536,768,0.),
                (96,1536,768,3e-5))]
    print(json.dumps({"platform":platform.platform(),"runs":results},indent=2,allow_nan=False))
