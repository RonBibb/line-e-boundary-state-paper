"""Find a pair of commensurate radial eigenvalues with asymmetric phase."""
import math
import json
import numpy as np
from scipy.special import j0,y0
from scipy.optimize import brentq


def roots(ratio, count=5):
    x = np.linspace(0.01,(count+1)*np.pi/(1-1/ratio),4096)
    f = lambda z:j0(z/ratio)*y0(z)-y0(z/ratio)*j0(z)
    values = f(x)
    return [brentq(f,a,b,xtol=1e-14) for a,b,fa,fb in
            zip(x,x[1:],values,values[1:]) if fa*fb<0][:count]


def candidate_ratio(ratio):
    values = roots(ratio)
    return values[-1]/values[0]


if __name__=="__main__":
    ratios = (1.5,2.,4.,8.,16.,32.,64.,128.,1024.,1000000.)
    root = brentq(lambda log_R:candidate_ratio(math.exp(log_R))-6,math.log(2),25.,xtol=1e-13)
    R = math.exp(root)
    k = roots(R)
    a,b = k[0]/R,k[0]
    print(json.dumps({"bracket_values":{str(r):candidate_ratio(r) for r in ratios},
                      "endpoints_roots":{str(x):roots(x) for x in (64.,1e10)}, "R":R,"a":a,"b":b,"roots":k,"J0_inner":[j0(a),j0(6*a)]},
                     indent=2,allow_nan=False))
