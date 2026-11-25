# Corrections to the 2026 Gowdy boundary-state preprint

This note accompanies a proposed successor to *Boundary-State Sufficiency in Gowdy Cosmology: Smooth Closure Does Not Extend Through Spikes*. The public preprint and its supporting archive establish the historical baseline. [Historical preprint](https://doi.org/10.5281/zenodo.21959052) · [Historical code release](https://doi.org/10.5281/zenodo.21874202)

## Einstein–scalar normalization and geometric scores

The displayed positive quadratures for the evolved metric potential require the transverse coefficient $e^{(\tau-\lambda)/2}$. The preprint and smooth scoring code instead used $e^{(\tau+\lambda)/2}$. The metric displayed by Lim uses a plus sign with negative quadratures. The reduced scalar was also rescaled. With $8\pi G=1$ and canonical scalar $\varphi$, use $\phi=\sqrt{2}\varphi$ in the equal-coefficient quadratures. Its normal energy density is $(\phi_\tau^2+e^{-2\tau}\phi_\theta^2)/(4N^2)$, where $N=e^{(-\lambda-3\tau)/4}$. The field equations and the existence of a free additive constant are unaffected by these changes to the metric and scoring formulas. The successor derives the correct equations from $R_{ab}=\nabla_a\varphi\nabla_b\varphi$. [Lim Eqs. (1), (4)–(5)](https://arxiv.org/pdf/0710.0628) · [Positive-quadrature convention](https://arxiv.org/pdf/0811.1266)

The benchmark now estimates the additive constant from the earlier metric reading at $\tau_f-2$ using the two-point field estimate. The numerical evidence therefore supports a pre-endpoint reconstruction with one independent metric datum. The table compares corrected observables on fixed areal slices, with density-scale mismatch measured on those same slices.

| Depth | Corrected score at $N=1024$, pre-endpoint metric | Score after setting limiting metric constant to zero |
|---:|---:|---:|
| 8 | $4.759\times10^{-4}$ | $0.23265$ |
| 12 | $1.060\times10^{-5}$ | $0.23157$ |
| 16 | $9.786\times10^{-6}$ | $0.23157$ |
| 24 | $9.786\times10^{-6}$ | $0.23157$ |

At depth 8 the same-data numerical behavior is a finite-window floor. Its corrected normalized observables give a fine-grid pass at the prior threshold of $10^{-3}$. The preprint's earlier threshold failure used the erroneous score. At depths 16 and 24 the corrected comparison improves by a factor of almost four on each grid doubling. The shape score is insensitive to an additive metric shift. It is also insensitive to an additive shift of $P$, which changes the lengths of the marked orbit cycles. Score agreement alone therefore leaves geometric degrees of freedom untested. [Corrected results](geometry_checks.json)

The preprint's numerical Bessel WINDOW arm reads the field at the restart point. A distinct exact-data arm initializes from the analytic coefficients. The new benchmark includes a genuinely withheld endpoint for Bessel fields and reports both cases. [Bessel records in the geometry results](geometry_checks.json)

## Meaning of the spike comparison

At a declared finite depth $f$, the archived tuple encodes the pair $P(f-2),P(f)$ along with $Q(f),Q_\tau(f)$ in a reversible algebraic parametrization. Replacing the missing $P_\tau(f)$ by the two-point slope is a decoder choice. On the centered Lim family the previously reported successful fit used quantities already available in that tuple and knowledge of the family. Off-grid checks additionally supplied the true Fourier shift to the restart. Neither result establishes intrinsic information loss in the four finite-depth arrays.

The published production numbers approximately $0.82$ refer to a fixed-depth reduced wave-map score on a mask of detector components. At the finest resolution the contemporaneous complementary-mask score peaked near $3.034$ at $\tau=2$. The detector-width numbers measure coordinate radians. Component counts can include multiple edges of a single feature. Neither decrease of the entire far-field score nor a count of geometric spikes follows. [Original data and scores](../../reproduction/inputs/bb-cbe3/CBE3.json) · [Original exact controls](../../reproduction/inputs/bb-cbe3t/CBE3T_R1.json)

The amended true-spike coefficient $-Q_0w$ already appears in the later primary formula. The original centered exact-spike model is rejected on the specified production center and training interval in the archived follow-up. In the successor the exact family is a local known-answer control only. [Published coefficient, Eq. (A4)](https://arxiv.org/pdf/0904.1546) · [Archived family rejection](../../reproduction/inputs/bb-cbe-limfit1/LIMFIT1_RESULT.json)

## Updated position

The research contribution is now an analysis of the map from finite measurements to a specified target, with direct benchmarks for the corrected Einstein–scalar benchmark and the nonlinear finite-window inverse. The sharp result concerns polarized Gowdy on a marked periodic torus. The independent review gives the full historical reference audit, including the originally omitted bidirectional Gowdy solver of Amorim, Bernardi and LeFloch. The new literature comparisons include the loaded string, the wave snapshot problem, singular scattering from full Cauchy data, and dynamical sampling. [Historical claim and citation audit](../../reviews/2026-11-25_novelty_validity_strategy/CITATION_MATRIX.md) · [Successor source](snapshot_observability.tex)

A new reconstruction matches the original observation window by solving for the missing derivative with the reduced PDE. On the original depth-18, $N=4096$ case, the same finite data and a fourth-order restart give a mask score of $0.00298$, compared with $0.82215$ for the secant and $0.00302$ for the fourth-order restart supplied with the endpoint velocity. An independent exact-template fit now estimates the subcell center from the samples by fitting a quadratic to $e^{-P}$. These remedies establish finite-window and known-family benchmarks at the stated settings. [New primary runs](nonlinear_checks.json)

## Follow-up to the production modulus and finite-window checks

The November continuation evaluated the historical depth-16 versus depth-24 check and compared each estimate to a separately computed continuum reference. The depth-gap prediction did not hold under refinement because the two depth estimates share the leading spatial discretization error. The earlier-anchor estimator with the Bessel tail converges at second order against the independent reference. The continuation also gives a rigorous finite-grid reconstruction basin for the archived velocity and exact nonlinear counterexamples at different observation placements. [Follow-up and evidence](extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/READOUT.md)
