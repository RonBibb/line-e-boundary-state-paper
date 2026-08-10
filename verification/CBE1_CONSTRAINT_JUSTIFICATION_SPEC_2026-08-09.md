# CBE1 Constraint-Justification Specification

**Date:** 2026-08-09  
**Test ID:** `E1-CJ1`  
**Status:** frozen before execution

## Question

For the implemented spatially flat homogeneous Einstein–massless-scalar
system, determine the exact propagation law of the Friedmann-constraint defect
and whether the Hamiltonian-constrained reduction is mathematically legitimate
without imposing the cross-branch comparison result.

## Definitions

With \(\pi=\dot\phi\), the free system is

\[
\dot a=aH,
\qquad
\dot H=-\frac12\pi^2,
\qquad
\dot\pi=-3H\pi.
\]

Define

\[
C=3H^2-\frac12\pi^2,
\qquad
\epsilon_C=\frac{|C|}{\rho},
\qquad
\rho=\frac12\pi^2.
\]

Analytically evaluate \(\dot C\) and the scaling of \(\epsilon_C\).

## Numerical arms

Integrate the expansion branch from

\[
t_0\in\{10^{-5},10^{-6},10^{-7}\}
\]

to \(t=1\), with the same asymptotic initialization and DOP853 tolerances
\(10^{-8}\) and \(10^{-11}\). Record:

- initial and final signed \(C\);
- maximum variation of \(C\);
- initial and final relative defect;
- fitted log-slope of \(\epsilon_C\) against \(a\); and
- dependence on the floating-point cancellation in the initial constraint.

The predicted continuum laws are

\[
\dot C=0,
\qquad
\rho\propto a^{-6},
\qquad
\epsilon_C\propto a^6
\]

for a nonzero constant numerical defect.

## Constrained reduction

For a chosen branch sign \(s=\pm1\), impose

\[
H=s\frac{|\pi|}{\sqrt6}
\]

and evolve only

\[
\dot a=aH,
\qquad
\dot\pi=-3H\pi.
\]

Verify algebraically that this implies

\[
\dot H=-\frac12\pi^2
\]

on the constraint surface. The reduction is legitimate if this identity
holds and if it contains no reference to the opposite branch or to the
comparison statistic \(D\).

## Frozen readings

- If \(C\) is constant and \(\epsilon_C\propto a^6\), retire the existing
  \(a^4\) curvature-error explanation.
- If \(C\propto a^{-2}\) and \(\epsilon_C\propto a^4\), retain the existing
  explanation.
- If neither closes, classify the formulation explanation as unresolved.

