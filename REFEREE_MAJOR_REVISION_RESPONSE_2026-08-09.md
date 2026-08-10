# Paper E1 Major-Revision Response

**Date:** 2026-08-09  
**Manuscript:** `PAPER_E1_MANUSCRIPT_DRAFT_01_2026-08-09.md`, revised to Draft 0.2

## Governing disposition

Both assessments recommended major revision before typesetting. All identified
scientific blockers have either been resolved or converted into explicit
limitations.

## Actions completed

1. **Endpoint-copy vulnerability.** Predictive holdout is now the primary
   finite-window result. The endpoint-window estimator is explicitly a lossy
   baseline. Full three-resolution holdout values and endpoint-prediction
   errors are reported. The retained collapse-side metric anchor is disclosed.

2. **Exact/window discrepancy.** A preserved bounded extension isolates the
   discrepancy to the finite-slice metric quadrature. At (N=512), using
   window field data with the asymptotic metric route changes the full distance
   from (1.4077\times10^{-4}) to (4.4696\times10^{-5}), equal to the exact
   arm at displayed precision. This quantitative error budget is now in the
   manuscript.

3. **Third FG1 resolution.** (N=1024) was added without altering the frozen
   FG1 record. Exact-arm full distance, lambda error, and momentum residual all
   exhibit two successive observed orders essentially equal to two.

4. **Constraint and lambda reporting.** Both exact and window momentum
   residual sequences are reported. The larger lambda amplitude is explained
   as accumulated differentiation/quadrature error and is shown to converge at
   second order.

5. **Global metric modulus.** The (c_\lambda) ablation is now the conceptual
   center of the paper. The manuscript derives why the local constant-shift
   coordinate transformation changes fixed torus periods and is not a
   period-preserving degree-one diffeomorphism.

6. **Scoring definitions.** The scored vector, area and density clock scales,
   identity weighting, shape and full norms, sampled surfaces, spatial
   supremum, and density-clock diagnostic are defined explicitly.

7. **Numerical formulation.** Initial data, periodic domain, second-order
   spatial operators, RK4 integrator, timestep rule, Fourier antiderivative,
   and lambda anchoring are stated.

8. **Multimode geometry.** The manuscript states that lambda is reconstructed
   from total multimode fields with all constraint-quadrature products, rather
   than assembled modewise.

9. **Internal vocabulary.** The internal FG1 classification string was removed
   while the malformed-gate disclosure and preserved original output remain.

10. **Literature.** The unfinished Beyer--LeFloch entry was repaired. Rendall--
    Weaver, Ringström, and Chruściel--Isenberg--Moncrief were added, and the
    spike limitation is tied to that literature.

11. **Independent verification.** The result-chain ledger now includes the
    third-resolution and metric-route decomposition. It passes 12/12 checks.

## Deliberately deferred

The suggested cross-program theorem about global or scale data disappearing
from normalized observables is not imported into Paper E1. It may motivate a
separate synthesis, but combining Palatini, junction, cosmological-boundary,
and galaxy-model examples would expand E1 beyond its frozen scientific scope.

Typesetting remains premature until the normalized action derivation,
figures, final bibliography verification, and author/repository metadata are
complete.

## Follow-up moderate revision

The second referee round identified a remaining seam between the asymptotic
modulus and the imported production anchor. Draft 0.2 has been further revised
as follows:

1. The manuscript now distinguishes
   \(c_\lambda=\lambda_\infty(0)\), the finite-depth anchor
   \(a_f=\lambda(\tau_f,0)\), and the reconstructed-slice anchor. It gives the
   formal estimator
   \(\widehat c_\lambda=a_f-\widehat q(0)\tau_f-\widehat r_\lambda\) and states
   that the production campaign does not estimate the remainder.
2. The principal ablation claim is now necessity/incompleteness, not recovery
   or five-part source-independent completeness.
3. Clock language now says “comparison on the areal clock, with density-scale
   consistency,” matching the implementation.
4. The exact benchmark coefficients, holdout width, and observation surface
   are tabulated.
5. A compact table reports every principal acceptance threshold.
6. The Jacobian dimensions, prescribed nulls, numerical ranks, smallest
   singular values, condition numbers, resolution stability, and linearity
   result are reported.
7. The homogeneous material is explicitly compressed into a validation
   summary, with reproducible details assigned to Appendix A.
8. Exceptional-data language is prospective; no unreported control is claimed.
9. The equality of endpoint and holdout distances is explained by the
   holdout's \(10^{-11}\)--\(10^{-12}\) field error lying far below the metric
   floor.
10. A bounded depth-24 route test confirms the common bottleneck. With the
    asymptotic spatial metric route and the same finite-depth anchor, the
    \(N=1024\) maximum distance is \(2.3432\times10^{-6}\) at both depths 16
    and 24. The prior depth growth was a finite-slice metric-quadrature effect.

The independent ledger now passes 13/13 checks. The remaining pre-typesetting
work is Appendix A, the normalized action derivation, figures, final
bibliographic verification, and submission metadata.

## Final scientific-draft revision

The final referee round was acted on as follows:

1. The manuscript-defined asymptotic metric construction
   \(\widehat q\tau_f+\widehat\lambda_\infty\) is now the primary production
   protocol. Its full three-resolution sequence is reported at both depths;
   the finite-slice quadrature is an error-isolation arm.
2. The common-anchor shift is explicit:
   \(\widehat c_\lambda^{(0)}=a_f-\widehat q(0)\tau_f\), with the periodic
   antiderivative zero at \(\theta=0\), so the reconstructed slice retains
   \(a_f\) without importing additional information.
3. A known-answer modulus estimator was run in the exact Bessel model. Without
   remainder subtraction its errors are \(6.13\times10^{-8}\),
   \(3.33\times10^{-11}\), and \(1.45\times10^{-14}\) at depths 8, 12, and
   16. With the analytic remainder it recovers \(c_\lambda=0.4\) to roundoff.
4. The abstract again states second-order convergence and distinguishes exact
   five-part recovery from the qualified production anchor.
5. The Jacobian parameterization, base points, perturbation, outputs, and
   sampling rule are now explicit. Its conditioning difference is interpreted.
6. Score terminology and the remaining clock criterion were normalized.
7. Both reported equation typos were corrected.

The independent result-chain ledger now passes 14/14 checks. No further
physics simulation is indicated before appendices, figures, archive, and
typesetting.

## Production modulus closure

The final requested bounded check was completed:

- the analytic benchmark's zero-remainder estimator errors decay at rates
  consistent with the expected \(\tau_f e^{-2\tau_f}\) remainder;
- the manuscript now states explicitly that this analytic-input test is
  resolution-independent and validates the formula rather than the numerical
  extraction chain;
- the production estimates converge spatially with difference ratios 3.9996
  at both depths;
- their Richardson values are 0.11803247427 and 0.11803247398 at depths 16
  and 24, differing by \(2.85\times10^{-10}\);
- no known production truth is asserted, and the depth gap is not claimed to
  converge spatially.

The manuscript now describes production-side modulus recovery as a spatially
convergent, depth-consistent estimate. The independent ledger passes 15/15
checks. The greater-framework publication checkpoint has also been updated:
Line E is scientifically closed and in manuscript production, independent of
ASC0-B and the exceptional spike limit.
