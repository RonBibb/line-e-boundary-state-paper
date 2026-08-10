# E1-RA2 Global Metric-Mode Classification Specification

**Date:** 2026-08-09  
**Test ID:** `E1-RA2`  
**Status:** frozen before numerical verification

## Question

Classify the additive Gowdy mode

\[
\lambda(\tau,\theta)\mapsto\lambda(\tau,\theta)+c.
\]

Determine whether it is removed by residual gauge, derived from the four
field asymptotic functions and constraints, or constitutes an additional
global metric/modulus datum when the \(T^3\) identifications are fixed.

## Analytic tasks

1. Derive constant coordinate transformations preserving the polarized Gowdy
   areal form.
2. Track the periods of \(\theta,x,y\).
3. Determine the transformation of lapse, directional expansion, scalar
   density, normalized observables, clocks, and proper cycle lengths.
4. Integrate the Hamiltonian and momentum quadratures and identify all
   undetermined zero modes.
5. Translate the result into a geometric singular-data interpretation.

## Bounded verification

Use the frozen CBE-2I evolution at \(N=512\), sampled at \(\tau=8\). Apply

\[
c\in\{-2,1,3\}
\]

without further evolution and verify:

\[
N_c/N_0=e^{c/4},
\quad
H_{i,c}/H_{i,0}=e^{-c/4},
\quad
\rho_c/\rho_0=e^{-c/2},
\]

\[
\Omega_c=\Omega_0,
\quad
\Sigma_c^2=\Sigma_0^2,
\quad
(s_P)_c=(s_P)_0,
\]

\[
(\Theta\lambda_A)_c/(\Theta\lambda_A)_0=e^{-c/4},
\quad
(\lambda_\rho)_c/(\lambda_\rho)_0=e^{c/4}.
\]

The momentum and Hamiltonian derivative equations must remain unchanged.

Tolerance: maximum relative discrepancy \(10^{-11}\), with an absolute floor
of \(10^{-13}\) for quantities crossing zero.

## Frozen readings

- **Gauge:** a period-preserving diffeomorphism removes \(c\).
- **Derived:** boundary constraints determine \(c\) from the four functions.
- **Independent modulus:** field equations determine only derivatives of
  \(\lambda\), and removing \(c\) changes one or more fixed torus periods.
- **Nonclosure:** analytic and numerical transformation laws disagree.

