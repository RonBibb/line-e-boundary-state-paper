"""Interval Bessel brackets for the fundamental and sixth-mode collision."""
import json
import math
import platform

import mpmath as mp

mp.iv.dps = 55
iv = mp.iv
I = iv.mpf


def symmetric_interval(radius):
    return I([-radius.b,radius.b])


def bessel_pair(x):
    y = x*x/4
    term = I(1)
    harmonic = I(0)
    regular = term
    harmonic_sum = I(0)
    for j in range(1,101):
        term *= -y/(j*j)
        harmonic += I(1)/j
        regular += term
        harmonic_sum -= harmonic*term
    next_term = abs(term*y/(101**2))
    series_ratio = y/(102**2)
    regular += symmetric_interval(next_term/(1-series_ratio))
    harmonic_sum += symmetric_interval(
        (harmonic+I(1)/101)*next_term/(1-2*series_ratio))
    singular = 2/iv.pi*((iv.ln(x/2)+iv.euler)*regular+harmonic_sum)
    return regular,singular


def determinant(z,r):
    j_outer,y_outer = bessel_pair(z)
    j_inner,y_inner = bessel_pair(z/r)
    return j_inner*y_outer-y_inner*j_outer


def report_interval(x):
    return [math.nextafter(float(x.a),-math.inf),
            math.nextafter(float(x.b),math.inf)]


def certify():
    brackets = {"64":(("2.850","2.852"),("15.737","15.739")),
                "10000000000":(("2.473","2.475"),("15.007","15.009"))}
    result = {}
    for r,(first,fifth) in brackets.items():
        signs = []
        for below,above in (first,fifth):
            below_value = determinant(I(below),I(r))
            above_value = determinant(I(above),I(r))
            signs.append(report_interval(below_value)+report_interval(above_value))
            assert (below_value.b<0<above_value.a or above_value.b<0<below_value.a)
        ratio = I([fifth[0],fifth[1]])/I([first[0],first[1]])
        result[r] = {"root_brackets":[first,fifth],"determinant_sign_intervals":signs,
                     "fifth_to_first_ratio_interval":report_interval(ratio)}
        if r == "64":
            assert ratio.b < 6
        else:
            assert ratio.a > 6
    return result


if __name__=="__main__":
    print(json.dumps({"platform":platform.platform(),"interval_digits":iv.dps,
                      "results":certify()},indent=2,allow_nan=False))
