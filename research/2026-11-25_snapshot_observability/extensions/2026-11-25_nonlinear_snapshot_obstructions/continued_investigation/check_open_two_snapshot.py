"""Finite-grid stress tests of two-slice recovery with orbit and metric data."""
import json
import math
import platform

import numpy as np
from scipy.integrate import solve_ivp

from check_nonlinear_sampling import spectral_operator, theta_grid
from check_null_transport import tuned_collision, weighted_output

DELAY = math.pi/16


def initial_fields(n, normal_velocity, polarized=False):
    theta = theta_grid(n)
    w = 0.9+0.24*np.cos(4*theta)
    transport = np.tanh(w)
    transverse = 1/np.cosh(w)
    speed = 0.16
    if polarized:
        return np.stack([w,np.full(n,speed),np.zeros(n),np.zeros(n)])
    normal = normal_velocity*(0.65+0.35*np.sin(3*theta))
    return np.stack([np.log(np.cosh(w)),speed*transport+normal*transverse,
                     transport,speed*transverse**2-normal*transverse*transport])


def forward_rhs(spatial, products, fields, t):
    P, Pt, Q, Qt = fields
    (Px,Qx),(Pxx,Qxx) = spatial(np.stack([P,Q]))
    e2P = np.exp(2*P)
    primary, other = products(np.stack([e2P*(Qt**2-Qx**2),-2*(Pt*Qt-Px*Qx)]))
    return np.stack([Pt,Pxx-Pt/t+primary,Qt,Qxx-Qt/t+other])


def background(start, a, c, tolerance, spatial, products):
    result = solve_ivp(lambda t,y:forward_rhs(spatial,products,y.reshape(start.shape),t).ravel(),
                       (a,c),start.ravel(),method="DOP853",dense_output=True,
                       rtol=tolerance,atol=tolerance*0.02)
    if not result.success:
        raise RuntimeError(result.message)
    return lambda t:result.sol(t).reshape(start.shape)


def tangent_matrix(start, a, c, K, tolerance=8e-11, return_operators=False):
    n = start.shape[-1]
    spatial, products = spectral_operator(n)
    base = background(start,a,c,tolerance,spatial,products)
    theta = theta_grid(n)
    modes = np.arange(1,K+1)
    phases = modes[:,None]*theta[None,:]
    amp = np.sqrt(2/(1+modes**2))[:,None]
    seeds = np.vstack((np.ones((1,n)),amp*np.cos(phases),amp*np.sin(phases)))

    def solve(seed_group):
        tangent = np.zeros((4,*seed_group.shape))
        tangent[1] = seed_group

        def rhs(t, flat):
            P, Pt, Q, Qt = base(t)
            (Px,Qx), _ = spatial(np.stack([P,Q]))
            p, pt, q, qt = flat.reshape(tangent.shape)
            (px,qx),(pxx,qxx) = spatial(np.stack([p,q]))
            rp, rq = products(np.stack([
                2*np.exp(2*P)*((Qt**2-Qx**2)*p+Qt*qt-Qx*qx),
                -2*(Qt*pt+Pt*qt-Qx*px-Px*qx)]))
            return np.stack([pt,pxx-pt/t+rp,qt,qxx-qt/t+rq]).ravel()

        run = solve_ivp(rhs,(a,c),tangent.ravel(),method="DOP853",rtol=tolerance,
                        atol=tolerance*0.01,t_eval=[c])
        if not run.success:
            raise RuntimeError(run.message)
        return run.y[:,-1].reshape(tangent.shape)

    P, Pt, Q, Qt = base(c)
    (Px,Qx),_ = spatial(np.stack([P,Q]))
    cutoff = min(n//3,K+48)
    p_columns, orbit_columns, metric_columns, partial_metric_columns = [], [], [], []
    for offset in range(0,len(seeds),12):
        p,pt,q,qt = solve(seeds[offset:offset+12])
        (px,qx),_ = spatial(np.stack([p,q]))
        lam_x = -2*c*(Px*pt+Pt*px+np.exp(2*P)*(Qx*qt+Qt*qx+2*Qx*Qt*p))
        p_record = weighted_output(p,2,cutoff)
        physical_q_record = weighted_output(np.exp(P)*q,2,cutoff)
        lam_record = weighted_output(lam_x,1,cutoff)
        orbit = np.column_stack([p_record,physical_q_record])
        p_columns.append(p_record)
        orbit_columns.append(orbit)
        metric_columns.append(np.column_stack([orbit,lam_record]))
        partial_metric_columns.append(np.column_stack([p_record,lam_record]))
    matrices = {name:np.concatenate(columns,axis=0).T for name,columns in
                (("only_P",p_columns),("full_orbit",orbit_columns),("full_metric",metric_columns),
                 ("P_and_metric",partial_metric_columns))}
    if return_operators:
        return matrices,P,cutoff
    result = {}
    for name,matrix in matrices.items():
        values = np.linalg.svd(matrix,compute_uv=False)
        result[name] = {"minimum_gain":float(values[-1]),"maximum_gain":float(values[0]),
                        "smallest_three":values[-3:].tolist()}
    (initial_Px,initial_Qx),_ = spatial(start[[0,2]])
    result["initial_rank_witness"] = float(np.max(np.abs(start[1]*initial_Qx-start[3]*initial_Px)))
    result["initial_momentum_average"] = float(np.mean(initial_Px*start[1]+
                                                      np.exp(2*start[0])*initial_Qx*start[3]))
    return result


def analytic_path(n, blend):
    collision, _ = tuned_collision()
    theta = theta_grid(n)
    a = collision.a
    start = collision.exact_fields(a,theta,1.)
    P, Pt, Q, Qt = start
    Qmean = float(np.mean(Q))
    R = max(3.2*np.max(np.abs(Q-Qmean)),math.exp(-np.mean(P))/math.sqrt(1-0.43**2))
    center = Qmean-0.43*R
    geodesic_P = -0.5*np.log(R**2-(Q-center)**2)
    P = (1-blend)*P+blend*geodesic_P
    Pt, Qt = (1-blend)*Pt, (1-blend)*Qt
    spatial, _ = spectral_operator(n)
    (Px,Qx),_ = spatial(np.stack([P,Q]))
    momentum = np.mean(Px*Pt+np.exp(2*P)*Qx*Qt)
    Qt -= momentum*Qx/np.mean(np.exp(2*P)*Qx**2)
    return a,np.stack([P,Pt,Q,Qt])


if __name__ == "__main__":
    results = {}
    for K,n in ((24,256),(48,384),(96,512),(96,768),(160,1024)):
        configurations = (("polarized",0.,True),("geodesic_image",0.,False),
                          ("rank_two",0.07,False))
        results[f"{K}_{n}"] = {name:tangent_matrix(initial_fields(n,normal,polarized),0.9,0.9+DELAY,K)
                              for name,normal,polarized in configurations}
    path = {}
    for blend in (0.,0.1,0.25,0.5,0.75,1.):
        a, initial = analytic_path(512,blend)
        path[str(blend)] = tangent_matrix(initial,a,a+math.pi/32,64)
    print(json.dumps({"platform":platform.platform(),"short_delay":DELAY,
                      "regular_examples":results,"analytic_path":path},indent=2,allow_nan=False))
