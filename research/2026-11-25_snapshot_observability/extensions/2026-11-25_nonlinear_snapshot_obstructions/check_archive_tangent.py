"""Measure the mixed-observation Jacobian on the archived Gowdy trajectory."""
import json
from pathlib import Path
import platform
from math import exp, factorial

import numpy as np
from scipy.integrate import solve_ivp

INPUTS = Path(__file__).resolve().parents[4] / "reproduction" / "inputs"
GAP = 2.0


def d1(f,dx):
    return (np.roll(f,-1)-np.roll(f,1))/(2*dx)


def d2(f,dx):
    return (np.roll(f,-1)-2*f+np.roll(f,1))/dx**2


def field_rhs(t,fields,dx):
    P,V,Q,R=fields
    a=np.exp(-2*t)
    x,y=d1(P,dx),d1(Q,dx)
    return np.stack([V,a*d2(P,dx)+np.exp(2*P)*(R**2-a*y**2),
                     R,a*d2(Q,dx)-2*(V*R-a*x*y)])


def coefficients(t,fields,dx):
    P,V,Q,R=fields
    tau_factor=np.exp(-t)
    e=np.exp(P)
    A=tau_factor**2/dx**2
    B=e*R
    C=tau_factor*e*d1(Q,dx)
    q_scale=tau_factor*C/dx
    lean=tau_factor**2*d1(P,dx)/dx
    minus=np.exp(P-np.roll(P,1))
    plus=np.exp(P-np.roll(P,-1))
    return V,A,B,C,q_scale,lean,minus,plus


def scaled_jacobian(coeffs,tangent):
    V,A,B,C,q_scale,lean,minus,plus=coeffs
    p,v,y,z=tangent
    A=float(A)
    V,B,C,q_scale,lean,minus,plus=(g[:,None] for g in (V,B,C,q_scale,lean,minus,plus))
    y_m,y_p=np.roll(y,1,axis=0),np.roll(y,-1,axis=0)
    p_m,p_p=np.roll(p,1,axis=0),np.roll(p,-1,axis=0)
    return np.stack([
        v,
        (2*(B**2-C**2)-2*A)*p+A*(p_m+p_p)+2*B*z+q_scale*(minus*y_m-plus*y_p),
        V*y+z,
        -V*z-2*B*v-2*A*y+(A-lean)*minus*y_m+(A+lean)*plus*y_p+q_scale*(p_p-p_m),
    ])


def generator_row_sums(coeffs):
    V,A,B,C,q_scale,lean,minus,plus=coeffs
    local=max(1.,float(np.max(np.abs(2*(B**2-C**2)-2*A)+2*np.abs(B))),
              float(np.max(np.abs(V)+1.)),float(np.max(np.abs(V)+2*np.abs(B)+2*A)))
    adjacent=max(float(np.max(2*A+np.abs(q_scale)*(plus+minus))),
                 float(np.max(np.abs(A-lean)*minus+np.abs(A+lean)*plus+2*np.abs(q_scale))))
    return local,adjacent


def projection(tape,dx,seeds,step):
    n,colors=seeds.shape
    initial=np.zeros((4,n,colors))
    initial[1]=seeds
    def rhs(t,y):
        f=y[:4*n].reshape(4,n)
        v=y[4*n:].reshape(4,n,colors)
        return np.concatenate((field_rhs(t,f,dx).ravel(),
                              scaled_jacobian(coefficients(t,f,dx),v).ravel()))
    sol=solve_ivp(rhs,(18.,16.),np.concatenate((tape.ravel(),initial.ravel())),
                  method="DOP853",rtol=1e-12,atol=3e-15,max_step=step,
                  t_eval=np.linspace(18.,16.,129))
    assert sol.success
    max_diagonal=max_adjacent=0.
    for t,reg in zip(sol.t,sol.y[:4*n,:].T):
        d,e=generator_row_sums(coefficients(t,reg.reshape(4,n),dx))
        max_diagonal=max(max_diagonal,d)
        max_adjacent=max(max_adjacent,e)
    endpoint=sol.y[:4*n,-1].reshape(4,n)
    output=-sol.y[4*n:,-1].reshape(4,n,colors)[0]/GAP
    return endpoint,output,max_diagonal,max_adjacent


def analyze(tape,earlier,dx,colors,step):
    n=tape.shape[1]
    assert n%colors==0
    idx=np.arange(n)
    seeds=np.eye(colors)[idx%colors]
    end,sampled,d,e=projection(tape,dx,seeds,step)
    radius=colors//2-1
    offsets=np.arange(-radius,radius+1)
    entries=sampled[idx[:,None],(idx[:,None]+offsets[None,:])%colors].copy()
    diag=entries[:,radius].copy()
    entries[:,radius]-=1.
    row_sums=np.sum(np.abs(entries),axis=1)
    # For exact uniform d,e in the scaled generator, this expression
    # bounds the tail in the row norm via the nearest-neighbor Dyson series.
    # Here d,e are sampled maxima, so the reported value is diagnostic.
    x=GAP*e
    tail=exp(GAP*d+x)*x**(radius+1)/factorial(radius+1)/GAP
    return {"colors":colors,"max_step":step,
            "maximum_row_defect":float(np.max(row_sums)),
            "largest_off_diagonal_sum":float(np.max(row_sums-np.abs(diag-1.))),
            "diagonal_range":[float(np.min(diag)),float(np.max(diag))],
            "maximum_row_defect_index":int(np.argmax(row_sums)),
            "sampled_local_generator_norm":d,"sampled_adjacent_generator_norm":e,
            "locality_tail_scale_from_sampled_bounds":tail,
            "archive_P_at_16_absolute_difference":float(np.max(np.abs(end[0]-earlier[0])))}


def directional_check(tape,dx):
    n=tape.shape[1]
    x=np.arange(n)*dx
    input_v=(0.7*np.cos(3*x)+0.3*np.cos(63*x))
    _,got,_,_=projection(tape,dx,input_v[:,None],0.04)
    eps=1e-5
    results=[]
    for sign in (1,-1):
        start=tape.copy()
        start[1]+=sign*eps*input_v
        sol=solve_ivp(lambda t,y: field_rhs(t,y.reshape(4,n),dx).ravel(),
                      (18.,16.),start.ravel(),method="DOP853",
                      rtol=1e-12,atol=3e-15,max_step=0.04)
        assert sol.success
        results.append(sol.y[:n,-1])
    fd=-(results[0]-results[1])/(2*eps*GAP)
    return float(np.max(np.abs(fd-got[:,0])))


if __name__=="__main__":
    rows={}
    for n in (1024,2048,4096):
        with np.load(INPUTS/"bb-cbe3b"/f"cbe3b_slices_N{n}.npz") as archive:
            tape,earlier=(np.stack([archive[f"tau{t}_{i}"] for i in range(4)]) for t in (18.0,16.0))
        dx=2*np.pi/n
        checks=[analyze(tape,earlier,dx,8,0.04)]
        if n==4096:
            checks.append(analyze(tape,earlier,dx,16,0.02))
        rows[str(n)]={"tangent_runs":checks,"FD_direction_absolute_gap":directional_check(tape,dx)}
    print(json.dumps({"platform":platform.platform(),"archive":"bb-cbe3b","results":rows},
                     indent=2,allow_nan=False))
