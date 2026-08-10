# E1-FG1 Exact Fuchsian and Geometric Benchmark Specification

**Date:** 2026-08-09  
**Test ID:** `E1-FG1`  
**Status:** frozen before execution

## Purpose

Test the Paper-E1 extraction/reconstruction model against an exact multi-mode
polarized Gowdy–massless-scalar solution whose singular asymptotic functions
are known analytically. Include the constrained metric variable \(\lambda\)
and its global mode.

## Exact fields

For each nonzero Fourier mode, use

\[
f_k(\tau,\theta)
=\left[A_kJ_0(ke^{-\tau})+B_kY_0(ke^{-\tau})\right]T_k(\theta),
\]

plus a zero mode \(f_0+v_0\tau\). Since

\[
Y_0(z)=\frac{2}{\pi}\left(\log(z/2)+\gamma\right)+o(1),
\]

the exact singular data are

\[
v_k=-\frac{2B_k}{\pi},
\qquad
f_{\infty,k}=A_k+\frac{2B_k}{\pi}
\left(\log(k/2)+\gamma\right).
\]

Use two nonzero modes in each of \(P\) and \(\phi\), with cosine modes for
\(P\) and sine modes for \(\phi\), so periodic momentum compatibility holds
mode by mode.

## Geometric data

Construct

\[
q=v_P^2+v_\phi^2,
\qquad
\lambda_\infty'=2(v_PP_\infty'+v_\phi\phi_\infty'),
\]

with frozen \(c_\lambda=0.4\). Determine the exact finite-depth
\(\lambda(\tau,0)\) from

\[
\lambda(\tau,0)=q(0)\tau+c_\lambda
-\int_\tau^\infty[\lambda_{,s}(s,0)-q(0)]ds,
\]

then recover its spatial profile from the exact momentum constraint.

## Arms

1. `FUCHSIAN`: initialize at \(\tau_f\) from exact two-term singular data.
2. `WINDOW`: estimate \((v,f_\infty)\) from exact fields on
   \([\tau_f-2,\tau_f]\), retaining the exact finite-depth metric anchor.

Evolve each arm from \(\tau_f\) to \(\tau=2\) with the frozen CBE-2I solver.

## Grid

\[
\tau_f\in\{8,12,16\},
\qquad
N\in\{256,512\}.
\]

## Scores

- endpoint asymptotic-data error;
- field and momentum sup errors at \(\tau=2\);
- \(\lambda\) sup error;
- shape-observable distance;
- full-observable distance;
- momentum-constraint residual.

## Gates

1. Exact asymptotic-data errors must decrease from depth 8 to 12 to 16.
2. At depth 16, field and geometric errors must decrease under resolution
   doubling.
3. Both arms must achieve full-observable distance below \(10^{-3}\) at
   depth 16, \(N=512\).
4. Fine-grid momentum residual must be below \(10^{-3}\).
5. The exact asymptotic data must satisfy periodic compatibility below
   \(10^{-12}\).

Passing supports compatibility of the E1 finite-depth method with the exact
Fuchsian data structure in this declared model. It is not a general
well-posedness theorem.

