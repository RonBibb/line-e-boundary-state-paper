# E1-F1 Boundary-State and Map Formalization

**Date:** 2026-08-09  
**Formalization ID:** `E1-F1`  
**Status:** controlling formal model for Paper E1 development  
**Supersedes:** the undifferentiated four-function completeness language in
the original CBE-2I result and Paper-E1 transfer

## 1. Purpose and scope

This document defines the mathematical objects actually tested by Paper E1.
It separates:

1. finite-time Cauchy data;
2. constrained asymptotic boundary data;
3. extraction and reconstruction maps;
4. geometric comparison observables;
5. time-parity conventions; and
6. normalized shape equivalence from full metric-state equivalence.

It does not prescribe a physical passage through the singular boundary. A
physical scattering or continuation law remains a separate map-selection
problem.

## 2. Declared polarized Gowdy system

The spacetime metric is

\[
ds^2=e^{(\lambda+\tau)/2}
\left(-e^{-2\tau}d\tau^2+d\theta^2\right)
+e^{-\tau}\left(e^Pdx^2+e^{-P}dy^2\right),
\]

on a marked \(T^3\), with fixed coordinate periods

\[
(L_\theta,L_x,L_y).
\]

The reduced Einstein–massless-scalar equations are

\[
P_{,\tau\tau}=e^{-2\tau}P_{,\theta\theta},
\qquad
\phi_{,\tau\tau}=e^{-2\tau}\phi_{,\theta\theta},
\]

\[
\lambda_{,\tau}
=P_{,\tau}^2+e^{-2\tau}P_{,\theta}^2
+\phi_{,\tau}^2+e^{-2\tau}\phi_{,\theta}^2,
\]

\[
\lambda_{,\theta}
=2(P_{,\tau}P_{,\theta}+\phi_{,\tau}\phi_{,\theta}).
\]

The normalization remains \(8\pi G=1\), with the scalar rescaled to enter the
reduced equations symmetrically with \(P\). A manuscript-level action
derivation remains required.

## 3. Finite-time Cauchy-data space

At a regular areal slice \(\tau=\tau_0\), define

\[
\mathscr C(\tau_0)
=\left{
(P,\Pi_P,\phi,\Pi_\phi,\lambda;
L_\theta,L_x,L_y)
\right},
\]

where

\[
\Pi_P=P_{,\tau},
\qquad
\Pi_\phi=\phi_{,\tau}.
\]

Admissible Cauchy data satisfy the momentum constraint

\[
\lambda_{,\theta}
=2(\Pi_PP_{,\theta}+\Pi_\phi\phi_{,\theta})
\]

and the global periodicity condition

\[
\oint
(\Pi_PP_{,\theta}+\Pi_\phi\phi_{,\theta}),d\theta=0.
\]

The Hamiltonian quadrature determines \(\lambda\)'s evolution after its
spatially constant mode is supplied.

We distinguish contracting and expanding copies

\[
\mathscr C_-(\tau_0),
\qquad
\mathscr C_+(\tau_0),
\]

by their physical future-normal orientation, not merely by the direction in
which the numerical integrator traverses \(\tau\).

## 4. Asymptotic boundary-data space

### 4.1 Shape data

For the declared smooth polarized AVTD form,

\[
P(\tau,\theta)
=v_P(\theta)\tau+P_\infty(\theta)+r_P,
\]

\[
\phi(\tau,\theta)
=v_\phi(\theta)\tau+\phi_\infty(\theta)+r_\phi,
\]

where the remainders decay in the admitted smooth sector.

Define

\[
\mathfrak B_{\rm shape}
=\left(v_P,v_\phi,P_\infty,\phi_\infty\right).
\]

These functions determine

\[
q(\theta)=v_P^2+v_\phi^2,
\]

and the leading metric asymptotics

\[
\lambda(\tau,\theta)
=q(\theta)\tau+\lambda_\infty(\theta)+r_\lambda.
\]

The asymptotic momentum constraint gives

\[
\partial_\theta\lambda_\infty
=2\left(v_P\partial_\theta P_\infty
+v_\phi\partial_\theta\phi_\infty\right).
\]

### 4.2 Admissibility condition

Periodic \(\lambda_\infty\) requires

\[
\boxed{
\oint
\left(v_P\partial_\theta P_\infty
+v_\phi\partial_\theta\phi_\infty\right)d\theta=0.
}
\]

The admissible shape space \(\mathscr B_{\rm shape}^{\rm adm}\) is the subset
of periodic function quadruples satisfying this condition and the smooth
AVTD-domain conditions required by the declared evolution.

### 4.3 Global metric datum

The constraint determines \(\lambda_\infty(\theta)\) only up to

\[
c_\lambda=\lambda_\infty(0).
\]

`E1-RA2` establishes that, with the marked torus periods fixed,
\(c_\lambda\) is an independent global metric/modulus datum rather than a
derived function of \(\mathfrak B_{\rm shape}\).

Define the full data object

\[
\mathfrak B_{\rm full}
=\left(
v_P,v_\phi,P_\infty,\phi_\infty,
c_\lambda;
L_\theta,L_x,L_y
\right).
\]

Only invariant combinations of \(c_\lambda\) and the torus periods are
physical under changes of marking. Paper E1 fixes the marking and treats
\(c_\lambda\) as the remaining relative base/orbit modulus.

The complete admissible boundary space is denoted

\[
\mathscr B_{\rm full}^{\rm adm}.
\]

## 5. Residual equivalences and symmetries

### 5.1 Period-preserving diffeomorphisms

Ordinary periodic translations and orientation reversals of the marked cycles
act on the data but do not alter their physical equivalence class when the
field profiles are transformed accordingly.

An arbitrary continuous rescaling of a periodic coordinate is not a
period-preserving degree-one diffeomorphism of the fixed marked circle.

### 5.2 Constant \(\lambda\) shift

Locally, \(\lambda\to\lambda+c\) is associated with

\[
\tau'=\tau+c/3,
\quad
\theta'=e^{-c/3}\theta,
\quad
x'=e^{c/6}x,
\quad
y'=e^{c/6}y.
\]

This changes the marked periods and therefore moves between different global
modulus data in the frozen model.

### 5.3 Scalar shift

The massless scalar admits

\[
\phi\mapsto\phi+\phi_c.
\]

Observables built only from derivatives cannot detect the constant scalar
mode. Paper E1 must state whether \(\phi_\infty\)'s zero mode is retained as a
matter boundary datum or quotiented by this global symmetry. Until a physical
coupling fixes an absolute scalar origin, shape scoring is shift blind.

### 5.4 Killing-direction relabeling

Interchanging \(x\leftrightarrow y\) sends

\[
P\mapsto-P,
\qquad
s_P\mapsto-s_P.
\]

Either the marked Killing directions must remain fixed or comparisons must
minimize over this discrete relabeling. The completed CBE-2I campaign used a
fixed marking.

## 6. Extraction map

At finite depth \(\tau_f\) and window length \(\Delta>0\), define the original
two-surface extraction

\[
\mathcal E_{\tau_f,\Delta}:
\mathscr C_-	o\mathscr B_{\rm full}^{\rm adm}
\]

by

\[
v_f^{(\tau_f)}
=\frac{f(\tau_f)-f(\tau_f-\Delta)}{\Delta},
\qquad
f_\infty^{(\tau_f)}
=f(\tau_f)-v_f^{(\tau_f)}\tau_f,
\]

for \(f\in\{P,\phi\}\), together with the extracted metric zero mode.

For predictive holdout, define

\[
\mathcal E^{\rm hold}_{\tau_f,\Delta}
\]

using only \([\tau_f-2\Delta,\tau_f-\Delta]\) and extrapolating to
\(\tau_f\). `E1-RA1` shows that this holdout map reproduces the tested smooth
polarized endpoint and retains the original second-order comparison residuals.

The extraction map is a finite-depth estimator of asymptotic data. It is not
identical to the limiting Fuchsian data map until

\[
\tau_f\to\infty,
\qquad
h\to0,
\]

with an appropriate remainder and conditioning estimate.

## 7. Reconstruction map

Given admissible boundary data and a regular initialization depth \(\tau_0\),
define

\[
\mathcal R_{\tau_0}:
\mathscr B_{\rm full}^{\rm adm}
\to\mathscr C_+(\tau_0)
\]

through

\[
P_0=v_P\tau_0+P_\infty,
\qquad
\Pi_{P0}=v_P,
\]

\[
\phi_0=v_\phi\tau_0+\phi_\infty,
\qquad
\Pi_{\phi0}=v_\phi,
\]

\[
\lambda_0=q\tau_0+\lambda_\infty,
\]

with \(\lambda_\infty\) obtained from its derivative constraint and
\(c_\lambda\). The regular equations are then evolved away from the boundary.

At finite \(\tau_0\), omitted asymptotic remainders produce a reconstruction
error. The present implementation uses the two-term VTD truncation.

## 8. Time orientation and parity map

Let \(n_{(+\tau)}\) be the unit normal associated with increasing \(\tau\).
For the collapsing branch, the physical future normal is

\[
n_-^{\rm future}=n_{(+\tau)}.
\]

For the reconstructed branch evolved away from the boundary toward decreasing
\(\tau\), the physical future normal is

\[
n_+^{\rm future}=-n_{(+\tau)}.
\]

The parity map \(\mathcal T\) preserves intrinsic data and reverses all data
linear in the future normal:

\[
\mathcal T:
(g_{ij},K_{ij},\phi,n\phi)
\mapsto
(g_{ij},-K_{ij},\phi,-n\phi).
\]

For the scored observables:

| Quantity | Parity |
|---|---|
| \(\rho\) | even |
| \(\Omega\) | even |
| \(\Sigma^2\) | even |
| \(s_P=(H_x-H_y)/\Theta\) | even |
| \(\Theta\) | odd |
| \(\Theta\lambda_A\) | odd |
| \(\lambda_A,\lambda_\rho\) | even |
| \(c_\lambda\) and intrinsic torus moduli | even |

`E1-RA1` verifies that the explicit future-normal reversal followed by
\(\mathcal T\) agrees with the original common-normal numerical comparison.

## 9. Observable map

Define

\[
\mathcal O:
\mathscr C_\pm
\to\mathscr V
\]

with two layers.

### 9.1 Normalized shape observables

\[
\mathcal O_{\rm shape}
=\left(\Omega,\Sigma^2,s_P,\ldots\right).
\]

These are insensitive to \(c_\lambda\) in the declared model.

### 9.2 Absolute metric/clock observables

\[
\mathcal O_{\rm abs}
=\left(
\lambda_A,
\lambda_\rho,
\Theta\lambda_A,
\text{proper cycle lengths},
\text{curvature scales}
\right).
\]

These require a common metric-modulus and unit normalization. In particular,
under \(c_\lambda\to c_\lambda+c\),

\[
\lambda_\rho\to e^{c/4}\lambda_\rho,
\qquad
\Theta\lambda_A\to e^{-c/4}\Theta\lambda_A.
\]

Constraint residuals are validation diagnostics and are not components of
either physical observable layer.

## 10. Equivalence relations

### 10.1 Shape equivalence

Two admitted branches are shape-equivalent at the boundary,

\[
\mathscr C_-\sim_{\rm shape}\mathscr C_+,
\]

when

\[
\lim_{\epsilon\to0}
d_{\rm shape}
\left(
\mathcal O_{\rm shape}[\mathscr C_-](\epsilon),
\mathcal T\mathcal O_{\rm shape}[\mathscr C_+](\epsilon)
\right)=0
\]

on the declared geometric clock surfaces and after the fixed discrete marking
conventions are applied.

This relation intentionally quotients the global metric scale detected by
`E1-RA1` and `E1-RA2`.

### 10.2 Full boundary-state equivalence

Two branches are fully boundary-state equivalent,

\[
\mathscr C_-\sim_{\rm full}\mathscr C_+,
\]

only if shape equivalence holds and the intrinsic metric/modulus data match:

\[
c_{\lambda,-}=c_{\lambda,+}
\]

under a common marking, or if an explicitly declared global isometry maps the
complete marked metric data. Absolute clocks and rates must then agree after
parity in common units.

### 10.3 Physical continuation

Neither relation selects

\[
\mathcal S:
\mathscr B_-	o\mathscr B_+.
\]

A physical continuation requires an independently derived scattering or
selection law. Equality under \(\sim_{\rm shape}\) or \(\sim_{\rm full}\) is
not such a law.

## 11. Composition tested by E1

The operational object is

\[
\mathcal P_{\tau_f,\Delta,\tau_0}
=\mathcal O\circ\mathcal R_{\tau_0}
\circ\mathcal E_{\tau_f,\Delta}.
\]

Paper E1 compares this with direct observables from the original contracting
branch after parity and clock matching.

For the tested smooth polarized realization:

- the holdout composition closes with second-order residual convergence for
  \(\mathcal O_{\rm shape}\);
- the full observable comparison closes only when \(c_\lambda\) is retained;
  and
- deleting \(c_\lambda\) leaves shape observables unchanged but destroys
  absolute clock/rate agreement.

## 12. Formal claims admitted after E1-F1

Paper E1 may state:

1. the reduced four-function data predictively reconstruct the normalized
   shape state for the declared smooth polarized realization;
2. a separate global metric mode is required for full state and clock
   reconstruction;
3. constraint admissibility restricts the function space;
4. explicit time parity is consistent with the implemented common-normal
   comparison; and
5. boundary-state equality remains distinct from a physical scattering law.

Paper E1 may not state:

1. that four functions completely specify the physical boundary state;
2. that the scale mode is gauge under fixed torus identifications;
3. that one tested realization proves injectivity or stability on an open
   family;
4. that the construction applies to exceptional unpolarized data; or
5. that a continuation through the singularity has been selected.

## 13. RA3 admission requirements

The next bounded test, `E1-RA3`, is admitted only as an open-family map audit.
It must:

1. define a preregistered family of smooth admissible Fourier data;
2. vary each shape-function sector independently while preserving periodic
   compatibility;
3. vary \(c_\lambda\) independently;
4. evaluate local sensitivity, conditioning, and distinguishability;
5. separate shape and absolute-observable Jacobians;
6. test whether any data directions are numerically null or redundant;
7. include more than one base realization; and
8. remain bounded to the smooth polarized model, with no ASC0/ASC1 execution.

Passing RA3 would support local injectivity and stability only over the tested
family. It would not establish a global theorem.

## 14. E1-F1 disposition

\[
\boxed{
\mathscr B_{\rm full}^{\rm adm}
=
\mathscr B_{\rm shape}^{\rm adm}
\times
\mathscr B_{\rm metric}
}
\]

with a constraint-qualified product rather than an unrestricted Cartesian
product, and with

\[
\boxed{
\sim_{\rm shape}
\;\text{strictly weaker than}\;
\sim_{\rm full}
\;\text{strictly weaker than a physical continuation law}.
}
\]

