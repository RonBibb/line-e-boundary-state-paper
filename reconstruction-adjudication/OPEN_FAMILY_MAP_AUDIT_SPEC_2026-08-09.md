# E1-RA3 Open-Family Map Audit Specification

**Date:** 2026-08-09  
**Test ID:** `E1-RA3`  
**Status:** frozen before execution

## 1. Question

Within a bounded smooth polarized family, is the map

\[
\mathcal O\circ\mathcal R:
\mathscr B_{\rm full}^{\rm adm}\to\mathscr V
\]

locally sensitive to each admitted boundary-data direction, with only the
symmetry/null directions predicted by `E1-F1`?

This is a local numerical rank and conditioning test, not a proof of global
injectivity.

## 2. Admissible Fourier family

Use

\[
v_P=a_0+a_2\cos2\theta,
\qquad
P_\infty=b_1\cos\theta,
\]

\[
v_\phi=c_0+c_2\cos2\theta,
\qquad
\phi_\infty=d_0+d_1\sin\theta.
\]

Parity makes

\[
\oint(v_PP_{\infty,\theta}+v_\phi\phi_{\infty,\theta})d\theta=0,
\]

so every parameter perturbation remains on the admissible constraint
manifold.

Include the global metric mode \(c_\lambda\).

The parameter vector is

\[
\mathbf p=(a_0,a_2,b_1,c_0,c_2,d_1,d_0,c_\lambda).
\]

Two base points are frozen:

\[
\mathbf p_A=(0.40,0.30,0.60,0.80,-0.20,0.50,0.10,0.50),
\]

\[
\mathbf p_B=(0.55,-0.18,0.35,0.68,0.14,-0.42,-0.20,-0.75).
\]

## 3. Reconstruction and observation

Construct the two-term VTD state at \(\tau_0=16\), including

\[
\lambda=q\tau_0+\lambda_\infty,
\qquad
q=v_P^2+v_\phi^2,
\]

and evolve to \(\tau=2\) using the frozen CBE-2I equations.

Run

\[
N\in\{256,512\}.
\]

Evaluate every 8th grid point at \(N=512\) and every 4th at \(N=256\), giving
the same 64 physical sample locations.

## 4. Perturbations

Use centered differences with

\[
\delta p_j=10^{-4}
\]

for every parameter. Verify linearity by repeating the base-A, \(N=256\)
Jacobian at \(\delta p_j=5\times10^{-5}\). Corresponding normalized Jacobian
columns must agree to 2% for active directions.

## 5. Observable layers

### Shape layer

\[
\mathcal O_{\rm shape}=(\Omega,\Sigma^2,s_P)
\]

at the 64 sample points.

### Full layer

Append

\[
\mathcal O_{\rm abs}=(\Theta\lambda_A,\log\lambda_\rho).
\]

Each output component is normalized by its baseline sup scale with a floor
of \(10^{-3}\). Each parameter column is expressed per unit fractional
change using \(\max(|p_j|,0.1)\).

## 6. Preregistered null and active directions

- \(d_0\) is a scalar-shift symmetry and must be null in both layers.
- \(c_\lambda\) must be null in the shape layer.
- \(c_\lambda\) must be active in the full layer.
- The remaining six shape directions must be active in both layers.

A direction is null if its normalized column norm is below \(10^{-8}\) and
active if above \(10^{-5}\).

## 7. Rank and conditioning gates

After deleting preregistered null columns:

- shape Jacobian expected rank: 6;
- full Jacobian expected rank: 7.

Numerical rank uses relative singular-value threshold \(10^{-7}\). The
smallest retained singular value and condition number must agree within a
factor of 2 between \(N=256\) and \(N=512\) for each base family.

Failure does not invalidate RA1 or RA2. It identifies local redundancy or
ill-conditioning in the chosen observable representation.

