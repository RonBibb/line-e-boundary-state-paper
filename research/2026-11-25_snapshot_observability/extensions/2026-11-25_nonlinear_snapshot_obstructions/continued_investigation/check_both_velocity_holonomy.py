"""Test fixed-time recovery of both velocity components in a circular loop."""
import json
import math
import platform

import numpy as np
from scipy.integrate import quad,solve_ivp

J = np.array([[0.,-1.],[1.,0.]])


def radial_background(a,c,winding,radius,tolerance):
    def rhs(t,y):
        R,v = y
        return v,-v/t-winding**2*math.sinh(R)*math.cosh(R)
    result = solve_ivp(rhs,(a,c),(radius,0.),method="DOP853",dense_output=True,
                       rtol=tolerance,atol=tolerance*0.02)
    assert result.success
    return result.sol


def vector_transfer(a,c,k,winding,radial,tolerance,return_terminal=False):
    start = np.zeros((2,len(k),2,2),dtype=np.complex128)
    start[1,:,:,:] = np.eye(2)

    def rhs(t,flat):
        z,zt = flat.reshape(start.shape)
        R,v = radial(t)
        omega = winding*math.cosh(R)
        potential = np.array([(winding*math.sinh(R))**2,-v*v])
        ztt = (-zt/t-(k[:,None,None]**2+omega*omega)*z
               +2j*k[:,None,None]*omega*(J@z)-potential[None,:,None]*z)
        return np.stack((zt,ztt)).ravel()

    result = solve_ivp(rhs,(a,c),start.ravel(),method="DOP853",rtol=tolerance,
                       atol=tolerance*0.005,t_eval=[c])
    assert result.success
    terminal, rate = result.y[:,-1].reshape(start.shape)
    if return_terminal:
        return terminal, rate
    scaled = (math.sqrt(c/a)*np.sqrt(1+k**2))[:,None,None]*terminal
    return np.linalg.svd(scaled,compute_uv=False)


def experiment(a,d,winding,radius,tolerance):
    radial = radial_background(a,a+d,winding,radius,tolerance)
    omega = quad(lambda x:winding*math.cosh(radial(x)[0]),a,a+d,epsabs=1e-13)[0]
    def principal(k):
        return np.minimum(abs(np.sin(k*d+omega)),abs(np.sin(k*d-omega)))
    sampled = np.arange(129)
    low = vector_transfer(a,a+d,sampled,winding,radial,tolerance)
    # Include high modes near the predicted most poorly observed phases.
    windows = ((128,256),(256,512),(512,1024),(1024,2048),(2048,4096))
    selected = []
    for left,right in windows:
        k = np.arange(left,right+1)
        selected.append(int(k[np.argmin(principal(k))]))
    high = vector_transfer(a,a+d,np.array(selected),winding,radial,tolerance)
    expected = principal(np.array(selected))
    return {"radius":radius,"delay":d,"winding":winding,"integrated_connection":omega,
            "connection_excess_over_center":omega-winding*d,
            "minimum_exact_gain_through_128":float(np.min(low[:,1])),
            "high_frequency_probes":[{"k":k,"principal_gain":float(e),
                                       "exact_gain":float(s[1])}
                                      for k,e,s in zip(selected,expected,high)]}


if __name__=="__main__":
    a,winding = 1.6,4
    rational = math.pi/16
    irrational = math.pi*math.sqrt(2)/24
    result = {"platform":platform.platform(),"a":a,
              "rational":[experiment(a,rational,winding,radius,2e-10)
                           for radius in (0.,0.08,0.16,0.32)],
              "irrational":experiment(a,irrational,winding,0.32,2e-10)}
    print(json.dumps(result,indent=2,allow_nan=False))
