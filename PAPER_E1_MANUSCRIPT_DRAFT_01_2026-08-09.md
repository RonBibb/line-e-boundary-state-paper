# Testing Singular Boundary-State Sufficiency in Gowdy Cosmology: Smooth Closure and Spike-Localized Failure

**Manuscript status:** Submission candidate  
**Target:** *Classical and Quantum Gravity*  
**Article class:** Regular research article  
**Author:** Ronald Bibb  
**Affiliation:** Independent Researcher, Lilburn, GA 30047, USA  
**Email:** ronbibb@gmail.com  
**ORCID:** 0009-0004-1153-2464

## Abstract

Can asymptotic data inferred from ordinary Cauchy evolution be shown
sufficient for a declared set of singular geometric observables? We formulate
a falsifiable extraction--reconstruction diagnostic: data are estimated on a
finite-depth window, the intervening history is discarded, a second branch is
initialized independently, and invariant observables are compared after a
declared time-parity map.

The answer depends sharply on the dynamical sector. In smooth polarized Gowdy
\(T^3\) Einstein--scalar data, predictive field reconstruction and exact
multimode Bessel tests converge at second order; four field functions recover
normalized shape, while full marked geometry requires the familiar additive
metric modulus. In nonlinear vacuum unpolarized Gowdy data, however, the same
logic gives a resolution-stable localized failure: the reconstruction distance
near two resolved spike candidates is
\((0.8201,0.8169,0.8222)\) for
\(N=(1024,2048,4096)\). Exact true-spike controls explain the failure.
Ordinary four-function velocity-dominated data lose the spike, reaching a
known-answer distance of \(0.5715\), whereas an enriched centered-spike state
restores reconstruction to \(2.38\times10^{-4}\). Analytic differentiation
verifies the exact control equations to \(3.11\times10^{-15}\), and unseen
shifted controls converge at second order and reproduce their centered
counterparts exactly under discrete periodic translation. The exceptional
status of Gowdy spikes is known analytically; the new result is a quantitative,
observable-relative measure of reconstruction loss and a controlled
identification of the enrichment that repairs it. Thus the method both closes in a weakly stressed
smooth sector and rejects an insufficient representation in a nonlinear
exceptional sector. Generic almost-everywhere sufficiency remains unresolved;
the result is the localized information-loss diagnosis and the identification
of a tested enriched state that repairs it on the exact centered family.

## 1. Introduction

Spacelike cosmological singularities are commonly studied by replacing
regular Cauchy data with asymptotic data defined at the singular boundary.
Fuchsian methods establish existence and, in suitable settings, uniqueness of
solutions possessing prescribed asymptotics [1--4]. Numerical Fuchsian
algorithms provide constructive approximations from such prescribed data
[5,6]. Geometric singular initial-value formulations and singularity
scattering maps further clarify which data may be assigned at a singular
hypersurface and how two sides may be related [7--9]. These results establish
that singular data can parameterize important classes of solutions. They do
not, by themselves, test whether a *particular reduced representation inferred
from finite-time Cauchy evolution* has retained all information needed to
recover a chosen geometric state.

We address that complementary validation problem. Starting from an ordinary
contracting Cauchy evolution, we extract a candidate boundary state, discard
the bulk history, independently initialize a branch directed away from the
boundary, and compare the two using geometric observables. We call this the
observable-relative extraction--reconstruction diagnostic. Its purpose is
diagnostic: it asks whether the declared representation is sufficient for the
declared class of observables and models. Passing is not a claim of injectivity
on the full geometric state space.

Three issues make the test nontrivial. First, formal time reversibility alone
cannot establish observable-relative completeness: an implementation may silently retain endpoint
fields, constraint data, or a global metric anchor. Second, equality of
coordinate-dependent fields can conceal disagreement in physical geometry.
We therefore score normalized geometric observables on geometric clocks and
state the parity induced by the future-normal reversal explicitly.
We consequently use matched positive, negative, ablation, convergence, and
constraint tests. The sharpest outcome is not the analytic existence of the
familiar additive metric mode. It is the sector dependence revealed by the
same diagnostic: smooth polarized data close, whereas a four-function
velocity-dominated representation fails locally in nonlinear unpolarized
spike data and is repaired only after the exact-family state is enlarged.

The contribution is the combined validation protocol

\[
\mathscr C_-
\xrightarrow{\;\mathcal E\;}
\mathscr B_{\rm adm}
\xrightarrow{\;\mathcal R\;}
\mathscr C_+
\xrightarrow{\;\mathcal O,\mathcal T\;}
\text{observable-relative residual},
\]

where \(\mathcal E\) extracts constrained boundary data from a regular Cauchy
solution, \(\mathcal R\) reconstructs a separately initialized branch,
\(\mathcal O\) produces comparison observables, and \(\mathcal T\) is the
declared time-parity map. The method differs from a numerical Fuchsian solver
in its inferential direction and validation target: the singular data are not
initially prescribed as truth but estimated from one Cauchy evolution and then
tested for their ability to reproduce another.

Accordingly, the novelty is not a claim to new Gowdy asymptotics or discovery
of the additive metric mode. It is a controlled inference about a proposed
boundary representation: the calculation distinguishes shape from marked
geometry, rejects the representation at resolved nonlinear spike candidates,
and identifies on exact controls which additional spike data restore the lost
solution. These outcomes are not consequences of formal reversibility alone.

The finite-depth logic has a close numerical-relativity analogue in waveform
extraction: measurements at finite radius must be separated from asymptotic
data, tested across extraction surfaces, and checked against an independent
evolution or extrapolation [15--17]. Here the limiting surface is singular
rather than null. We apply that logic to homogeneous controls, smooth
polarized Einstein--scalar cosmologies, exact polarized Bessel data, and a
nonlinear unpolarized vacuum calculation with exact true- and false-spike
controls.

The polarized evidence has three distinct independence levels. First, predictive
holdout extraction reconstructs the field sector without reading endpoint
field values. Second, the production metric reconstruction retains the
finite-depth collapse-side anchor \(a_f\) and is not a source-independent
recovery of the limiting metric modulus. Third, the exact benchmark supplies
the limiting five-part boundary object analytically and validates the full
forward reconstruction without importing a production endpoint anchor.

The paper is organized as follows. Section 2 defines the boundary spaces,
maps, parity, observables, and equivalence relations. Section 3 gives the exact
homogeneous controls. Section 4 introduces the polarized and unpolarized Gowdy
systems and their asymptotic data. Section 5 describes the numerical protocol
and acceptance criteria. Section 6 reports the smooth reconstruction,
geometric ablation, and spike-localized failure. Section 7
interprets the result relative to prior singular initial-value methods.
Section 8 states the limitations and conclusions.

## 2. Observable-relative boundary-state framework

### 2.1 Cauchy and boundary data

Let \(\mathscr C_-(\tau_0)\) and \(\mathscr C_+(\tau_0)\) denote regular
Cauchy-data spaces distinguished by physical future-normal orientation. In
the polarized Gowdy model considered below, a Cauchy state is

\[
\mathscr C(\tau_0)=
\{P,\Pi_P,\phi,\Pi_\phi,\lambda;
L_\theta,L_x,L_y\},
\qquad
\Pi_P=P_{,\tau},\quad \Pi_\phi=\phi_{,\tau}.
\]

The asymptotic shape data are

\[
\mathfrak B_{\rm shape}
=\left(v_P,v_\phi,P_\infty,\phi_\infty\right),
\]

with asymptotic expansions

\[
P=v_P\tau+P_\infty+r_P,
\qquad
\phi=v_\phi\tau+\phi_\infty+r_\phi.
\]

They determine \(q=v_P^2+v_\phi^2\) and

\[
\lambda=q\tau+\lambda_\infty+r_\lambda,
\qquad
\partial_\theta\lambda_\infty
=2(v_P\partial_\theta P_\infty
+v_\phi\partial_\theta\phi_\infty).
\]

Periodic admissibility requires

\[
\oint
(v_P\partial_\theta P_\infty
+v_\phi\partial_\theta\phi_\infty)\,d\theta=0.
\]

The constraint fixes \(\lambda_\infty\) only up to one constant
\(c_\lambda=\lambda_\infty(0)\). With the torus marking fixed, the full data
object is therefore

\[
\mathfrak B_{\rm full}=
(v_P,v_\phi,P_\infty,\phi_\infty,c_\lambda;
L_\theta,L_x,L_y)\in\mathscr B_{\rm full}^{\rm adm}.
\]

The constant scalar mode may be quotiented by the shift symmetry
\(\phi\mapsto\phi+\phi_c\) when only derivative observables are scored. The
present calculations retain a fixed scalar origin and a fixed marking of all
three torus cycles.

Here a marking means an ordered choice of primitive generators of
\(\pi_1(T^3)\), represented by the \(\theta\), \(x\), and \(y\) cycles together
with their declared periods \((L_\theta,L_x,L_y)\). Marked states are quotiented
only by diffeomorphisms that preserve those generators and periods (up to the
fixed orientation convention). A constant coordinate dilation changes the
identification lattice and therefore relates different marked states rather
than removing a degree of freedom within one state. If the marking and its
periods are relaxed, the associated overall scale direction may instead be
quotiented; the claim that \(c_\lambda\) is a modulus is explicitly relative to
the fixed-marking equivalence relation.

Three metric quantities must be distinguished. The singular modulus is

\[
c_\lambda=\lambda_\infty(0).
\]

At finite depth the evolved solution instead supplies the anchor
\(a_f=\lambda(\tau_f,0)\). An estimator of the limiting modulus would require

\[
\widehat c_\lambda(\tau_f)=a_f-\widehat q(0)\tau_f
-\widehat r_\lambda(\tau_f,0).
\]

The smooth calculation neglects the remainder
\(\widehat r_\lambda\) and uses
\(\widehat c_\lambda^{(0)}=a_f-\widehat q(0)\tau_f\). Its validity is assessed
by spatial convergence and agreement between extraction depths, since no
production truth value is known. In the exact benchmark, by contrast,
\(c_\lambda\) and the remainder are known analytically.

### 2.2 Extraction and reconstruction

The singular limit in the areal coordinate used below is
\(\tau\rightarrow+\infty\). Reconstruction is initialized at finite
\(\tau_f\) and evolved toward smaller \(\tau\).

Our primary finite-window estimator is predictive. It uses only the holdout
interval \([\tau_f-2\Delta,\tau_f-\Delta]\):

\[
v_f^{\rm hold}=
\frac{f(\tau_f-\Delta)-f(\tau_f-2\Delta)}{\Delta},
\qquad
f_\infty^{\rm hold}=f(\tau_f-\Delta)
-v_f^{\rm hold}(\tau_f-\Delta),
\]

and predicts \(f(\tau_f)=v_f^{\rm hold}\tau_f+f_\infty^{\rm hold}\).
It does not read the endpoint field. The calculation retains the collapse-side
finite-depth anchor \(a_f=\lambda(\tau_f,0)\); it is therefore
independent in the field-extraction sector, but not an independent recovery of
every full-metric datum.

For comparison, the elementary endpoint-window baseline is

\[
v_f^{(\tau_f)}=
\frac{f(\tau_f)-f(\tau_f-\Delta)}{\Delta},
\qquad
f_\infty^{(\tau_f)}=f(\tau_f)-v_f^{(\tau_f)}\tau_f,
\]

for \(f\in\{P,\phi\}\). Because this formula reproduces the endpoint field
algebraically, it is only a lossy endpoint-compression baseline and is not
principal evidence for field extraction. The exact benchmark of Sec. 6.3 uses
analytic limiting Fuchsian data and contains neither endpoint copy nor an
imported finite-depth metric anchor.

Given admissible boundary data, the two-term reconstruction at depth
\(\tau_0\) is

\[
P_0=v_P\tau_0+P_\infty,
\quad \Pi_{P0}=v_P,
\qquad
\phi_0=v_\phi\tau_0+\phi_\infty,
\quad \Pi_{\phi0}=v_\phi,
\]

\[
\lambda_0=q\tau_0+\lambda_\infty.
\]

The regular field equations are then evolved away from the boundary. At finite
depth, the omitted remainders contribute asymptotic truncation error.

### 2.3 Time parity and comparison observables

Reversal of the future-directed unit normal acts as

\[
\mathcal T:(g_{ij},K_{ij},\phi,n\phi)
\mapsto(g_{ij},-K_{ij},\phi,-n\phi).
\]

Thus \(\rho\), \(\Omega\), \(\Sigma^2\), and the normalized polarization
shape \(s_P=(H_x-H_y)/\Theta\) are even; \(\Theta\) is odd. Products such as
\(\Theta\lambda_A\), where \(\lambda_A\) is a boundary-distance clock, are
compared only after the corresponding orientation convention is applied.

We distinguish a normalized shape layer

\[
\mathcal O_{\rm shape}=(\Omega,\Sigma^2,s_P,\ldots)
\]

from an absolute metric/clock layer containing \(\lambda_A\),
\(\lambda_\rho\), \(\Theta\lambda_A\), proper cycle lengths, and curvature
scales. Shape observables are insensitive to \(c_\lambda\); absolute
observables are not. Constraint residuals are numerical-health diagnostics,
not physical state-vector components.

### 2.4 Observable-relative completeness criterion

In the Gowdy calculation the scored vector and two clock scales are

\[
\mathfrak O_s=(\Omega,\Sigma^2,s_P,\Theta\lambda_A),
\qquad \lambda_A=e^{-\tau/2},
\qquad \lambda_\rho=\rho^{-1/2}.
\]

The first clock is the symmetry-orbit area scale \(\lambda_A\), equivalently
\(\tau=-2\log\lambda_A\) after normalization of the marked Killing
coordinates. The independently constructed density scale \(\lambda_\rho\) is
used as a clock-consistency diagnostic. On each sampled areal surface
\(\tau\in\{2,4,8,12\}\), define

\[
D(\tau)=\sup_\theta
\left\|
\mathfrak O_{s,-}(\tau,\theta)
-\mathcal T\mathfrak O_{s,+}(\tau,\theta)
\right\|_2,
\]

with the identity weight matrix: the first three entries are dimensionless,
and \(\Theta\lambda_A\) is dimensionless under the fixed normalization. The
reported distance is the maximum over all periodic grid points and all four
surfaces. We separately report

\[
D_\rho=\max_{\tau,\theta}
\frac{|\lambda_{\rho,-}-\lambda_{\rho,+}|}{\lambda_{\rho,-}}.
\]

For the exact benchmark, the shape distance uses
\((\Omega,\Sigma^2,s_P)\); the full distance adds
\(\Theta\lambda_A\) and
\(\log(\lambda_{\rho,-}/\lambda_{\rho,+})\) in the same pointwise Euclidean
norm. Thus no unspecified weight enters a principal result. A candidate
representation passes the observable-relative numerical test in a stated family if
the residual decreases consistently under refinement, remains below the
specified criterion, agrees on the areal clock with density-clock consistency,
and passes compatible, incompatible, constraint, and data-ablation controls.
Passing is a sufficient validation result only for the declared model family,
marking, reconstruction map, observable vector, clocks, parity convention,
norm, and tolerances. It neither proves that the representation is injective on
the full geometric state space nor excludes information loss in an unmeasured
observable.

## 3. Homogeneous controls

The homogeneous suite serves to validate the comparison machinery rather than
to carry the novelty claim. A stiff, spatially flat FLRW collapse paired with
its time-reversed expansion is the exact positive control. The comparison
closes on both the areal and density clocks, with fine residuals
\(D_{\lambda_\rho}=4.98\times10^{-7}\) and
\(D_{\lambda_R}=3.10\times10^{-7}\).

A Schwarzschild/Kantowski--Sachs interior compared with stiff FLRW supplies an
exact adverse control. It terminates at the invariant-mismatch stage because
the vacuum-versus-matter Ricci structure, normalized shear, and Weyl fraction
cannot agree. This demonstrates that the protocol is not a generic
time-reversal matcher.

For homogeneous Einstein--massless-scalar solutions with anisotropy, closure
is stratified by the conserved normalized shear
\(\Sigma^2=\sigma^2/(6H^2)\). States in different shear strata are not
equivalent even when their volume behavior is similar. The constrained
homogeneous evolution is an exact reduction of the Hamiltonian constraint
after the expansion/contraction sign is chosen. These controls are not
principal evidence; their state variables, clock interpolation, mismatch
vector, and analytic derivations are collected in Appendix A.

## 4. Gowdy models and asymptotic states

### 4.1 Polarized Einstein--scalar sector

We use the areal-time metric

\[
ds^2=e^{(\lambda+\tau)/2}
(-e^{-2\tau}d\tau^2+d\theta^2)
+e^{-\tau}(e^Pdx^2+e^{-P}dy^2)
\]

on a marked \(T^3\), in units \(8\pi G=1\), with the scalar normalized to
enter symmetrically with \(P\). The reduced equations are

\[
P_{,\tau\tau}=e^{-2\tau}P_{,\theta\theta},
\qquad
\phi_{,\tau\tau}=e^{-2\tau}\phi_{,\theta\theta},
\]

\[
\lambda_{,\tau}=P_{,\tau}^2+e^{-2\tau}P_{,\theta}^2
+\phi_{,\tau}^2+e^{-2\tau}\phi_{,\theta}^2,
\]

\[
\lambda_{,\theta}
=2(P_{,\tau}P_{,\theta}+\phi_{,\tau}\phi_{,\theta}).
\]

The evolution of \(P\) and \(\phi\) is free, while \(\lambda\) is recovered
by constrained quadrature, including its global mode. This use of constrained
evolution avoids treating drift in a redundant free subsystem as physical
disagreement.

For a Fourier mode with \(z=ke^{-\tau}\), the exact solution is

\[
f_k(\tau,\theta)=
[A_kJ_0(z)+B_kY_0(z)]T_k(\theta),
\]

and its two-term singular data are

\[
v_k=-\frac{2B_k}{\pi},
\qquad
f_{\infty,k}=A_k+\frac{2B_k}{\pi}
[\log(k/2)+\gamma].
\]

This identity supplies a known-answer bridge between our operational boundary
state and standard Fuchsian asymptotic data.

### 4.2 Unpolarized vacuum sector and spike data

The nonlinear extension uses vacuum unpolarized Gowdy \(T^3\), the setting in
which numerical studies first exposed the small-scale structures and
high-velocity spike dynamics relevant here [19,20]. In the same areal
coordinate convention, the dynamical fields \((P,Q)\) form a wave map
into the hyperbolic plane, with equations

\[
P_{,\tau\tau}=e^{-2\tau}P_{,\theta\theta}
+e^{2P}\left(Q_{,\tau}^2-e^{-2\tau}Q_{,\theta}^2\right),
\]

\[
Q_{,\tau\tau}=e^{-2\tau}Q_{,\theta\theta}
-2\left(P_{,\tau}Q_{,\tau}-e^{-2\tau}P_{,\theta}Q_{,\theta}\right),
\]

\[
\lambda_{,\tau}=P_{,\tau}^2+e^{-2\tau}P_{,\theta}^2
+e^{2P}\left(Q_{,\tau}^2+e^{-2\tau}Q_{,\theta}^2\right),
\]

\[
\lambda_{,\theta}=2\left(P_{,\tau}P_{,\theta}
+e^{2P}Q_{,\tau}Q_{,\theta}\right).
\]

Away from exceptional points, the tested velocity-dominated representation is

\[
\mathfrak B_{\rm VTD}=(v,P_\infty,Q_\infty,\pi_Q),
\qquad \pi_Q=e^{2P}Q_{,\tau},
\]

with

\[
P\sim v\tau+P_\infty,
\qquad
Q\sim Q_\infty-\frac{e^{-2P}\pi_Q}{2v}.
\]

This four-function state is not assumed to remain sufficient at a spike. The
exact controls use the true- and false-spike families generated by Gowdy
solution transformations [11,18]. Their role is twofold: they validate the
nonlinear evolution and detector against known answers, and they test whether
the ordinary VTD compression preserves the localized feature. The enriched
centered-spike representation retains the exact-family parameters
\((w,Q_0,Q_2,x_0)\) in addition to the background data. It is tested as a
family-specific repair, not asserted as a universal spike state.

## 5. Numerical protocol and validation design

Spatial derivatives of evolved fields use centered second-order periodic
differences. Time evolution uses classical fourth-order Runge--Kutta with
segmentwise steps satisfying \(|\Delta\tau|\leq0.2\Delta\theta\) and landing
exactly on every stored surface. The spatially varying part of \(\lambda\) is
reconstructed by a Fourier antiderivative after removal of the zero Fourier
mode; its additive anchor is supplied separately. All fields are periodic on
\(0\leq\theta<2\pi\).

The production initial data at \(\tau=0\) are

\[
P=0.6\cos\theta,\quad P_{,\tau}=0.4+0.3\cos2\theta,
\]

\[
\phi=0.5\sin\theta,\quad
\phi_{,\tau}=0.8-0.2\cos2\theta,
\]

with \(\lambda\) obtained from the periodic momentum quadrature and anchored
at zero.

The predictive holdout width is \(\Delta=2\): data are fitted on
\([\tau_f-4,\tau_f-2]\) and extrapolated to \(\tau_f\).

The exact benchmark parameters are:

**Table 1. Exact multimode Bessel benchmark parameters.**

| Field | Homogeneous constant | Homogeneous velocity | \(k\) | Parity | \(A_k\) | \(B_k\) |
|---|---:|---:|---:|---|---:|---:|
| \(P\) | 0.2 | 0.35 | 2 | cosine | 0.30 | 0.12 |
| \(P\) | 0.2 | 0.35 | 3 | cosine | -0.15 | 0.05 |
| \(\phi\) | -0.1 | 0.65 | 1 | sine | 0.25 | -0.08 |
| \(\phi\) | -0.1 | 0.65 | 4 | sine | 0.10 | 0.04 |

Here \(c_\lambda=0.4\), the extraction depths are \(8,12,16\), and the
observation surface is \(\tau=2\).

The principal numerical acceptance criteria are:

**Table 2. Principal numerical acceptance and adverse-control criteria.**

| Check | Decision rule |
|---|---|
| Observable-relative distance | fine-grid maximum \(<10^{-3}\) |
| Momentum constraint | fine-grid normalized residual \(<10^{-3}\) |
| Density-clock consistency | maximum relative gap \(<10^{-2}\) |
| Refinement | distance decreases at each resolution doubling |
| Periodic compatibility | absolute integral defect \(<10^{-12}\) |
| Parity consistency | transformed and common-normal scores agree within \(10^{-12}\) |
| Ablation failure | maximum distance \(>10^{-3}\) or nondecreasing under refinement |
| Jacobian directions | null norm \(<10^{-8}\), active norm \(>10^{-5}\) |
| Jacobian rank | relative singular-value threshold \(10^{-7}\) |
| Jacobian linearity | normalized columns agree within 2% under halved perturbation |

The momentum threshold applies to the designated fine grid; the exact test's
coarse \(N=256\) residual, \(1.1162\times10^{-3}\), is therefore a reported
coarse-grid value rather than a failed acceptance criterion.

The original exact-data refinement test required strict decrease. The
exact-input calculation returned \((0,0,0)\), so its implementation was
corrected to accept either exact closure or strict decrease. No tolerance,
datum, numerical parameter, reconstruction, or reported value changed; both
machine records are preserved in the archive.

The smooth calculation uses two extraction depths,
\(\tau_f=16,24\), and three spatial resolutions,
\(N=256,512,1024\). The exact Fuchsian/geometric benchmark uses depths
\(8,12,16\), initially at \(N=256,512\), and compares reconstructed fields
and geometry at \(\tau=2\). The acceptance criteria were fixed before result
interpretation. A bounded extension adds \(N=1024\) and decomposes the two
metric-initialization routes without altering the original result.

The unpolarized calculation uses

\[
P(0)=0,\quad P_{,\tau}(0)=2+1.5\cos\theta,
\quad Q(0)=\cos\theta,\quad Q_{,\tau}(0)=0,
\]

at \(N=(1024,2048,4096)\), with extraction depths \(\tau_f=12,18\) and
comparison surfaces \(\tau=2,4,8\). The four-function VTD state is extracted
over a width-2 window and evolved backward independently. The pointwise score
combines hyperbolic target-space distance, hyperbolic-velocity error, and
normalized \(\lambda_{,\tau}\) error. Localized components must span at least
four cells; a resolution-stable component score above \(10^{-2}\) is adverse.
Exact transformed-unpolarized, false-spike, homogeneous, and momentum-
constraint convergence tests must pass before the generic calculation is
interpreted.

The exact-spike tests use \(N=(1024,2048,4096)\). Analytically differentiated
field equations require residual below \(10^{-12}\), with the literal printed
formula alternatives retained as negative controls. Ordinary VTD loss is
compared with an enriched centered-family reconstruction at \(\tau=2,4\).
Unseen translated cases are preserved as independent tests of periodic-domain
robustness. Translation is the discrete circle action: the centered periodic
representative and every exact comparison slice are shifted by the same
integer number of cells. Centered and shifted errors must agree, remain finite,
converge at second order, and fall below \(10^{-2}\) on the finest grid.

The validation hierarchy is:

1. exact compatible and incompatible homogeneous controls;
2. manufactured-solution and constraint-convergence controls;
3. spatial refinement at fixed extraction depth;
4. extraction-depth comparison;
5. predictive holdout extraction;
6. ablation of the global metric datum;
7. local rank and null-direction tests of the boundary-to-observable map;
8. exact analytic Fuchsian and geometric reconstruction;
9. nonlinear unpolarized reconstruction near resolved localized features; and
10. exact true/false-spike information-loss and enriched-state controls.

The preserved results show compatible and incompatible controls behaving as
designed, second-order spatial refinement in the smooth and exact benchmarks,
predictive holdout below the declared distance criterion, constraint
convergence, failure under metric-anchor ablation, stable local ranks, and
closure of the metric-route decomposition. The nonlinear extension and its
footprint-aware successor are independently recomputed from archived
machine-readable results. Code, specifications, machine-readable results, and
correction records are publicly available.

## 6. Results

The three independence levels govern the interpretation of every result below.
The predictive holdout is endpoint-independent only in the field sector. The
production metric leg imports \(a_f\), a finite-depth collapse-side datum. Only
the exact benchmark supplies \(c_\lambda\) and the complete five-part limiting
boundary object without importing that production anchor.

### 6.1 Predictive finite-window reconstruction

The predictive holdout, not the endpoint-copy baseline, is the primary
finite-window result. Following the reconstruction defined in Sec. 2.2, its
primary metric slice is

\[
\lambda_{\rm rec}=\widehat q\tau_f+\widehat\lambda_\infty,
\qquad
\widehat c_\lambda^{(0)}=a_f-\widehat q(0)\tau_f,
\]

where the periodic antiderivative defining
\(\widehat\lambda_\infty-\widehat c_\lambda^{(0)}\) is zero at \(\theta=0\).
Consequently \(\lambda_{\rm rec}(\tau_f,0)=a_f\): the asymptotic spatial route
and the finite-depth reconstruction use exactly the same anchor.

At depth 16 the maximum smooth-calculation distances for
\(N=(256,512,1024)\) are

\[
(3.7477961\times10^{-5},
9.3715322\times10^{-6},
2.3431847\times10^{-6}).
\]

The absolute depth-16/depth-24 differences at the same resolutions are

\[
(1.5283\times10^{-11},
1.0007\times10^{-11},
2.2968\times10^{-11}),
\]

so the curves agree below \(2.30\times10^{-11}\), rather than merely after
rounding to five significant figures. Both refinement ratios are 4.000 to the
reported precision. The two-term truncation remainder is thus below the
spatial discretization floor at these depths.

The largest endpoint field-prediction
errors are \(6.4\times10^{-11}\) at depth 16 and
\(2.8\times10^{-12}\) at depth 24. The smooth field sector is therefore only
weakly stressed in this family.

A bounded depth-and-norm test moves the same smooth family
into a less asymptotic regime. It preserves the predictive estimator,
asymptotic metric route, and imported finite-depth anchor, but uses extraction
depths \(8,12,16,24\) and scores every case on the common observation surfaces
\(\tau=2,4\). At \(N=1024\), the identity-weight distances are

\[
2.58457\times10^{-3},\quad
1.76861\times10^{-6},\quad
8.69854\times10^{-7},\quad
8.69878\times10^{-7},
\]

respectively. The depth-8 case fails the \(10^{-3}\) criterion, with
\(\|\delta\Omega\|_\infty=2.43070\times10^{-3}\) and an endpoint \(P\)
prediction error of \(2.20225\times10^{-4}\). The diagnostic therefore detects
measurable two-term truncation rather than closing automatically in this
family. Depth 12 marks the transition to the spatial-error floor; depths 16
and 24 are consequently described as cross-depth consistent, not as an
empirical measurement of a depth-convergence rate.

The same calculation reports componentwise maxima and three alternative norms. At
depths 16 and 24, the identity, component-maximum, shape-emphasized weights
\((2,2,2,1/2)\), and clock-emphasized weights \((1/2,1/2,1/2,2)\) all remain
below \(1.74\times10^{-6}\) on the fine grid and converge with refinement
ratios within \(0.1\%\) of four. At depth 12 the fine-grid component maxima in
the order \((\Omega,\Sigma^2,s_P,\Theta\lambda_A)\) are

\[
(1.55630\times10^{-6},\ 8.63388\times10^{-7},\
1.87371\times10^{-7},\ 8.44130\times10^{-7}).
\]

Thus the deep-production pass is not an artifact of identity weighting, while
the shallow failure remains visible under the identity, shape-emphasized, and
component-maximum scores.

As an error-isolation test, we also reconstruct \(\lambda\) directly from
differentiated finite-depth endpoint fields and a periodic quadrature, while
retaining the same \(a_f\). This finite-slice route gives
\((2.0089\times10^{-4},5.0213\times10^{-5},1.2553\times10^{-5})\) at depth
16 and
\((2.9665\times10^{-4},7.4144\times10^{-5},1.8535\times10^{-5})\) at
extraction depth 24. Its \(\lambda\) error at depth 16 decreases from
\(1.64\times10^{-3}\) to \(1.03\times10^{-4}\). The larger magnitude and
depth dependence are therefore properties of this alternate metric
quadrature, not of the declared reconstruction or the extracted field data.
The endpoint-copy field estimator is retained separately as a lossy baseline.

### 6.2 Shape versus full marked geometry: the global modulus

In the exact benchmark, the five-part data object
\((v_P,v_\phi,P_\infty,\phi_\infty,c_\lambda)\) reproduces both normalized
shape and absolute metric observables. The corresponding production object
uses \(a_f\) in place of the independently recovered limiting modulus.
Removing the global datum leaves the shape
layer unchanged but produces order-unity-to-large nonclosure in the absolute
layer: maximum residuals lie near 29 at depth 16 and 106 at depth 24 and do
not decrease with resolution. Direct transformations
\(c_\lambda\mapsto c_\lambda+c\) leave normalized shape observables invariant
to numerical precision while rescaling absolute clock and metric quantities
according to the analytic laws. With fixed marked periods, no nonzero shift is
generated by a period-preserving diffeomorphism. Thus \(c_\lambda\) is an
independent global metric/modulus datum in the declared model.

This is the principal *observable-relative incompleteness* finding: the operational protocol
detects that the four field asymptotic functions can reconstruct the
normalized shape layer while failing for the absolute marked geometry if the
global metric datum is omitted. The ablation establishes necessity of that
datum. The exact benchmark identifies it with \(c_\lambda\); the production
holdout estimates its limiting value from \(a_f\) and \(\widehat q(0)\), with
depth consistency replacing a known-answer comparison. The periodic momentum
constraint determines only \(\partial_\theta\lambda_\infty\), so it cannot
determine \(c_\lambda\). Locally, a shift \(\lambda\mapsto\lambda+c\) can be
accompanied by

\[
\tau'=\tau+c/3,\quad \theta'=e^{-c/3}\theta,
\quad x'=e^{c/6}x,\quad y'=e^{c/6}y.
\]

For \(c\ne0\), this changes at least one fixed coordinate period. It is not a
degree-one, period-preserving diffeomorphism of the declared marked torus.
Thus \(c_\lambda\) may be quotiented as a scale direction only if rescaling the
torus marking is admitted; with \((L_\theta,L_x,L_y)\) fixed, it is a modulus
of the declared state.

The fixed-marking dependence can be read directly from proper cycle lengths
on a constant-\((\tau,\theta)\) orbit:

\[
\ell_\theta=L_\theta e^{(\lambda+\tau)/4},\qquad
\ell_x=L_xe^{(-\tau+P)/2},\qquad
\ell_y=L_ye^{(-\tau-P)/2}.
\]

At fixed areal clock and fixed field data, \(\lambda\mapsto\lambda+c\)
multiplies \(\ell_\theta\) by \(e^{c/4}\). The polarization-independent marked
ratio

\[
\frac{\ell_\theta^2}{\ell_x\ell_y}
=\frac{L_\theta^2}{L_xL_y}e^{(\lambda+3\tau)/2}
\]

therefore changes by \(e^{c/2}\). This quantity isolates the modulus within
the fixed identification lattice. If the lattice scale is not part of the
state, the same change can be absorbed into a redefinition of the marking and
is not asserted to be an unmarked observable.

For the local boundary-to-observable Jacobian we use

\[
v_P=a_0+a_2\cos2\theta,\quad P_\infty=b_1\cos\theta,
\]
\[
v_\phi=c_0+c_2\cos2\theta,\quad
\phi_\infty=d_0+d_1\sin\theta,
\]

and the eight-entry vector
\(\mathbf p=(a_0,a_2,b_1,c_0,c_2,d_1,d_0,c_\lambda)\). The two base points
are

\[
\mathbf p_A=(0.40,0.30,0.60,0.80,-0.20,0.50,0.10,0.50),
\]
\[
\mathbf p_B=(0.55,-0.18,0.35,0.68,0.14,-0.42,-0.20,-0.75).
\]

Centered perturbations use \(\delta p_j=10^{-4}\). At \(N=256\) every fourth
grid point is sampled, and at \(N=512\) every eighth, giving the same 64
physical locations. The shape outputs are
\((\Omega,\Sigma^2,s_P)\); the full layer adds
\((\Theta\lambda_A,\log\lambda_\rho)\).

As an independent local sufficiency check, the boundary/observable Jacobian
changes rank from 6 for normalized-shape outputs to 7 when absolute
marked-geometric outputs are included; the added active direction is precisely
\(c_\lambda\). After the scalar-shift and normalized-scale null directions are
accounted for, the rank,
direction, resolution-stability, and linearity criteria all pass. Specifically,
the parameter vector has eight entries. The shape output contains three fields
at 64 sample locations (192 components), with the scalar offset \(d_0\) and
\(c_\lambda\) as prescribed null directions; its measured rank is 6. The full
output adds two fields (320 components total), activates \(c_\lambda\), and has
rank 7, with only \(d_0\) null. At the two base points, the smallest retained
singular values are \(0.1712\) and \(0.0345\) for shape and \(0.5079\) and
\(0.6804\) for full data at \(N=512\). Corresponding condition numbers are
121.8, 564.6, 44.9, and 30.7; each changes by less than 0.03% between
\(N=256\) and 512. Halving the parameter perturbation changes normalized
Jacobian columns by \(2.5\times10^{-8}\), well below the 2% criterion.
The shape layer is more poorly conditioned because its three normalized
outputs leave several active field directions nearly degenerate; the two
absolute outputs separate those directions and lower the full-layer condition
number. This independently reflects the distinction between normalized-shape
and marked-geometric information.

### 6.3 Exact Fuchsian and geometric benchmark

The benchmark uses genuine multimode sums: two nonzero modes occur separately
in each of \(P\) and \(\phi\), and \(\lambda\) is reconstructed from the total
fields. All products generated by the constraint quadrature are retained; the
geometry is not assembled modewise.

For exact analytic singular data at depth 16, refinement through
\(N=256,512,1024\) gives observed orders between 1.9988 and 2.0000 for the
principal errors.
The errors in \(P\), \(P_{,\tau}\), \(\phi\), and \(\lambda\) decrease from
\(4.32\times10^{-6}\), \(8.16\times10^{-6}\),
\(3.29\times10^{-6}\), and \(4.96\times10^{-4}\) to
\(1.08\times10^{-6}\), \(2.04\times10^{-6}\),
\(8.22\times10^{-7}\), and \(1.24\times10^{-4}\), respectively, and then to
\(2.70\times10^{-7}\), \(5.10\times10^{-7}\),
\(2.05\times10^{-7}\), and \(3.10\times10^{-5}\). The full distance decreases
from \(1.79\times10^{-4}\) through \(4.47\times10^{-5}\) to
\(1.12\times10^{-5}\). The momentum residual falls from
\(1.1162\times10^{-3}\) through \(2.7925\times10^{-4}\) to
\(6.9824\times10^{-5}\). The larger absolute \(\lambda\) and constraint errors
arise from differentiation and accumulated spatial quadrature, but their
second-order convergence distinguishes them from a failure of geometric
reconstruction.

For window-extracted data, the maximum asymptotic-data error falls from
\(2.31\times10^{-5}\) at depth 8 to \(1.88\times10^{-8}\) at depth 12 and
\(1.16\times10^{-11}\) at depth 16. At depth 16 and \(N=512\), the shape and
full distances are \(1.94\times10^{-5}\) and \(1.41\times10^{-4}\). Periodic
compatibility closes between exact zero and \(2.84\times10^{-16}\). The window
arm momentum residual falls from \(3.3979\times10^{-3}\) through
\(8.5006\times10^{-4}\) to \(2.1259\times10^{-4}\).

The same known-answer solution tests recovery of the global modulus. With
\(a_f=\lambda(\tau_f,0)\), the finite-depth estimator that neglects the
remainder,

\[
\widehat c_\lambda^{(0)}=a_f-q(0)\tau_f,
\]

approaches the known \(c_\lambda=0.4\) with absolute errors
\(6.1277\times10^{-8}\), \(3.3291\times10^{-11}\), and
\(1.4544\times10^{-14}\) at depths 8, 12, and 16. The successive error ratios,
1841 and 2289, agree with the ratios 1987 and 2236 predicted by a
\(\tau_f e^{-2\tau_f}\) remainder over the same depth increments. This is a
consistency check on the assumed two-term expansion. Subtracting the analytic
remainder recovers \(c_\lambda\) to roundoff.

This estimator test uses analytic values for \(q(0)\), the anchor, and the
remainder; its rows are therefore resolution-independent. It validates the
formula and the asymptotic remainder rate, not the numerical extraction chain.
The latter is tested separately in the production holdout. There,
\(\widehat c_\lambda^{(0)}\) at \(N=(256,512,1024)\) is

\[
(0.11804140434,\ 0.11803470693,\ 0.11803303243)
\]

at depth 16 and

\[
(0.11804140406,\ 0.11803470667,\ 0.11803303215)
\]

at depth 24. Successive spatial-difference ratios are 3.9996 at both depths.
Second-order Richardson estimates are 0.11803247427 and 0.11803247398, a gap
of \(2.85\times10^{-10}\). Because the production truth is unknown, this is a
depth-consistent recovery estimate rather than an absolute-error measurement;
the depth gap itself is not claimed to converge spatially.

The exact/window full-distance gap is accounted for by the two discrete metric
routes. At \(N=512\), retaining the window-extracted field data but replacing
only the finite-slice metric quadrature by \(q\tau_f+\lambda_\infty\) changes
the full distance from \(1.4077\times10^{-4}\) to
\(4.4696\times10^{-5}\), equal to the exact-data calculation at displayed precision.
At \(N=1024\), the corresponding values are \(3.5197\times10^{-5}\) and
\(1.1175\times10^{-5}\). The four-function data error is therefore consistent
with the exact calculation; the excess is a separate second-order metric-quadrature
contribution amplified by \(\tau_f\partial_\theta v\) terms in endpoint
gradients.

The exact-closure implementation correction is stated with the acceptance
criteria in Sec. 5; both machine records remain archived.

### 6.4 Nonlinear unpolarized extension: localized insufficiency

The unpolarized extension asks a question for which the smooth polarized
answer is not assumed: does the ordinary four-function VTD state preserve a
developing localized feature after the finite history has been discarded?
All launch controls pass. An exact \(SL(2,\mathbb R)\)-transformed solution
converges at orders 2.0002 and 2.0001. The exact false-spike control converges
at orders 1.9976 and 1.9994 and is classified as transient rather than true.
The homogeneous conserved quantity \(\pi_Q\) drifts by
\(1.51\times10^{-12}\), and the production momentum residual falls from
\(7.27\times10^{-3}\) through \(1.84\times10^{-3}\) to
\(4.61\times10^{-4}\), with orders 1.98 and 2.00.

At \(\tau_f=18\), two symmetric localized components are found. Their widths
are 5, 10, and 19 cells on \(N=(1024,2048,4096)\), corresponding to physical
angular widths \((0.03068,0.03068,0.02915)\) radians. They are therefore
resolved rather than single-grid artifacts. On the common \(\tau=2\) surface,
the four-function reconstruction distances are

\[
D_{\rm loc}=(0.82012,0.81689,0.82215).
\]

The spread is below 1% and the values exceed the adverse criterion by more
than a factor of 80. This is not a truncation-floor result: the proposed VTD
state is insufficient in the localized region. The original ordinary-region
supremum is not interpreted because the fixed mask is contaminated by the
wider feature footprint; generic almost-everywhere sufficiency therefore
remains unresolved.

Exact true-spike calculations identify the missing-information mechanism.
Analytic differentiation of the amended closed forms gives a maximum field-
equation residual of \(3.11\times10^{-15}\); literal printed alternatives
give residuals as large as \(3.59\) and serve as negative controls. Numerical
evolution of the true and false exact families converges at second order. If
the exact true spike is compressed into ordinary four-function VTD data, its
known-answer center distance reaches \(0.57150\). Retaining the enriched
centered-family parameters \((w,Q_0,Q_2,x_0)\) reduces the maximum distance to
\(2.3812\times10^{-4}\) at \(\tau=2\) and
\(2.2421\times10^{-5}\) at \(\tau=4\). Thus the diagnostic does more than
flag a mismatch: on a known exact family it distinguishes an insufficient
state from a tested state that preserves the spike.

The translated-control qualification is now closed. The earlier implementation
re-evaluated a nonperiodic local Lim formula as \(x-x_0\), thereby changing its
artificial seam jump instead of translating a periodic grid state. In the
corrected test, periodic translation is an integer cyclic shift applied to the
centered representative and to every exact comparison slice. For the lower-
velocity unseen family, the finest \((\tau=2,4)\) errors are
\((2.6369\times10^{-3},2.2809\times10^{-4})\); for the higher-velocity family
they are \((2.1388\times10^{-6},1.9881\times10^{-7})\). Both sequences converge
at second order over \(N=(1024,2048,4096)\), remain finite, and have zero
centered--shifted error gap at each resolution. The generic production feature is therefore called a resolved spike
candidate, not asserted to be a member of the exact family. A subsequent
dual-scale footprint calculation localized the exceptional region but could
not establish zero-measure decay within its sealed window; deeper runs stopped
when the narrowing feature fell below the constraint-reconstruction
resolution criterion. Those stopped outcomes are retained as numerical
limits, not converted into evidence for equivalence.

## 7. Relation to Fuchsian and geometric singular-data methods

The analytical literature already establishes that spikes are exceptional
points for ordinary velocity-dominated asymptotics. Rendall and Weaver
construct spike solutions, Lim supplies exact spike families, and
Ringström's strong-cosmic-censorship analysis treats the exceptional set
separately [11,18,13]. We therefore do not claim that failure of ordinary VTD
data at a spike is new. The contribution here is operational and quantitative:
the declared observable score measures the reconstruction loss, demonstrates
its resolution stability in a nonlinear evolution, and uses exact controls to
identify which enriched data repair that loss. The translated-control campaign
further verifies covariance under the discrete periodic translation used by
the numerical scheme.

The present construction does not compete with Fuchsian existence theory on
its strongest ground. Fuchsian theory begins with admissible singular data and
establishes a development with controlled remainders. A numerical Fuchsian
algorithm approximates that development. Our primary question runs in the
opposite inferential direction: starting with a finite-time Cauchy solution,
does a proposed compression into singular data preserve enough information to
recover the selected invariant state?

The numerical structure is analogous to finite-radius gravitational-wave
extraction. In that setting, near-zone contamination is diagnosed by varying
the extraction radius, extrapolating to infinity, or supplying worldtube data
to an independent characteristic evolution [15--17]. Here extraction depth
replaces radius and the candidate asymptotic state replaces worldtube data.
The analogy is methodological rather than geometric, but it supplies a mature
precedent for treating finite-surface contamination, discretization error, and
independent reconstruction as separate uncertainties.

The exact benchmark makes this distinction precise. In the polarized linear
field sector, the four functions
\[
(v_P,v_\phi,P_\infty,\phi_\infty)
\]
are exactly the coefficients of the analytic two-term Fuchsian expansion. The momentum
constraint supplies the spatially varying part of \(\lambda_\infty\), while
\(c_\lambda\) supplies its independent global mode in the exact benchmark.
The smooth calculation estimates the modulus by combining the finite-depth
anchor \(a_f\) with extracted \(\widehat q(0)\) and neglecting the asymptotic
remainder; convergence and cross-depth agreement test that estimate. Hence the operational
boundary object is a constrained geometric data set rather than four
unrestricted functions. Reconstruction from exact data validates the forward
leg; convergence of finite-window estimates toward those data validates the
extraction leg.

The smooth-sector contribution beyond recognizing familiar Fuchsian field
coefficients is an operational test of the separation between field/shape
data and full marked-geometric data. The
constraint recovers the nonzero modes of \(\lambda_\infty\), but ablation shows
that a global mode must be retained independently. This necessity is familiar
in constrained metric reconstruction and is not claimed as a newly discovered
degree of freedom. The stronger result comes from the nonlinear extension:
the same extraction--reconstruction test returns a converged localized
failure, and exact controls trace that failure to information removed by the
ordinary VTD state. Enriched spike data repair the exact family, and translated
controls reproduce the centered errors exactly while converging at second
order. In the smooth
sector, the global mode's load-bearing role under predictive extraction,
reconstruction, ablation, and geometric scoring remains useful supporting
evidence. It is
invisible to the normalized shape score and decisive for absolute clocks and
lengths.

The method also differs from a singularity scattering law. We test the
identity-plus-parity relation as a diagnostic hypothesis. We do not derive or
select a physical law mapping one side of a singularity to another. In
particular, agreement of boundary states does not imply geodesic extension, a
bounce, topology change, or transmission of matter through the singularity.

**Table 3. Relation of the present diagnostic to neighboring singular-data frameworks.**

| Feature | Fuchsian theory/algorithm | Scattering-map framework | Present diagnostic |
|---|---|---|---|
| Input data | Prescribed singular data | Data on one or both singular sides | Extracted from regular Cauchy evolution |
| Primary question | Existence/development | Classification or prescription of a two-sided map | Sufficiency for a declared observable set |
| Independent reconstruction | Forward construction from singular data | Depends on selected map | Required validation leg |
| Scoring | Solution/remainder error | Geometric map conditions | Pointwise residual on the areal clock plus density-clock consistency |
| Failure controls | Method-dependent | Not the present numerical protocol | Positive, negative, ablation, constraint, holdout |
| Physical continuation selected | No | Potentially prescribed | No |

## 8. Limitations and conclusion

The evidence remains model-scoped. It does not establish a theorem over the
full polarized or unpolarized function spaces, and the numerical schemes are
second order. In the generic unpolarized calculation, the resolved localized
feature has not been proved to belong to the exact Lim family. The original
ordinary-region mask is footprint-contaminated, and deeper production calculations
stop as the narrowing feature becomes underresolved. Consequently this paper
does not claim almost-everywhere sufficiency, zero-measure exceptional support,
or a universal spike boundary state. BKL-oscillatory regimes remain untested.
The diagnostic also selects no physical continuation, causal or thermodynamic
arrow, topology change, or transmission law through the singularity.

Within these limits, the results support four conclusions. First,
boundary-state sufficiency can be posed as a falsifiable numerical property of
an extraction--reconstruction chain rather than inferred from time
reversibility. Second, normalized shape sufficiency and full marked-geometric
sufficiency are distinct in the smooth polarized model: the latter requires a
global metric datum relative to the fixed marking. Third, the smooth field
holdout and exact-Fuchsian reconstruction converge consistently, but this
positive result is weakly stressed because the analytic remainder is already
below the spatial error at the deepest surfaces. Fourth, and most importantly,
ordinary four-function VTD data fail by a resolution-stable margin near
nonlinear unpolarized spike candidates. Exact controls reproduce the same
information loss and show that enriched centered-spike data restore the tested
family.

The central result is therefore sector dependent: smooth closure does not
license exceptional closure. A boundary representation must be tested against
the observables and localized structures it is intended to preserve.

## Appendix A. Normalization and homogeneous controls

### A.1 Action and conventions

We use signature \((-+++)\) and units \(8\pi G=1\). For a massless scalar,
the Einstein--scalar action is

\[
S[g,\varphi]=\int d^4x\,\sqrt{-g}\left[
\frac12 R-\frac12 g^{\mu\nu}\nabla_\mu\varphi\nabla_\nu\varphi
\right].
\]

Its Euler--Lagrange equations are \(G_{\mu\nu}=T_{\mu\nu}\) and
\(\nabla^\mu\nabla_\mu\varphi=0\), with

\[
T_{\mu\nu}=\nabla_\mu\varphi\nabla_\nu\varphi
-\frac12g_{\mu\nu}(\nabla\varphi)^2.
\]

For the polarized Gowdy ansatz of Sec. 4, integration over the two Killing
directions, removal of boundary terms, and a constant scalar rescaling
\(\varphi\mapsto\phi\) give, up to an overall constant torus-volume factor,

\[
S_{\mathrm{red}}=\frac12\int d\tau\,d\theta\,
\left[P_{,\tau}^2-e^{-2\tau}P_{,\theta}^2
+\phi_{,\tau}^2-e^{-2\tau}\phi_{,\theta}^2\right].
\]

The scalar rescaling is fixed here by the requirement that \(P\) and \(\phi\)
enter the metric constraints with equal coefficients. Variation gives the two
wave equations in Sec. 4. The remaining Einstein equations are the Hamiltonian
and momentum quadratures for \(\lambda\). Their compatibility follows directly:
differentiating \(\lambda_{,\tau}\) with respect to \(\theta\), differentiating
\(\lambda_{,\theta}\) with respect to \(\tau\), and using the wave equations
gives \(\partial_\theta\lambda_{,\tau}=\partial_\tau\lambda_{,\theta}\).
Consequently the local gradients determine \(\lambda\) only up to the additive
global mode tested in Secs. 6.2--6.3.

### A.2 Exact stiff-FLRW positive control

For the spatially flat metric \(ds^2=-dt^2+a(t)^2d\mathbf{x}^2\) and a
homogeneous scalar, the field equations reduce to

\[
3H^2=\rho=\frac12\dot\varphi^2,
\qquad
\dot H=-\frac12\dot\varphi^2,
\qquad
\ddot\varphi+3H\dot\varphi=0.
\]

After the sign \(\epsilon=\operatorname{sgn}H\) is chosen, constrained
evolution uses \(H=\epsilon|\dot\varphi|/\sqrt6\). With proper-time distance
\(s\) from the singular boundary, the exact stiff solution obeys

\[
a\propto s^{1/3},\qquad |H|=\frac1{3s},\qquad
\rho=\frac1{3s^2},\qquad \dot\varphi^2\propto a^{-6}.
\]

The curvature invariants used to construct the independent clocks are

\[
R_{\mu\nu}R^{\mu\nu}
=9(\dot H+H^2)^2+3(\dot H+3H^2)^2,
\]

\[
R_{\alpha\beta\gamma\delta}R^{\alpha\beta\gamma\delta}
=12\left[(\dot H+H^2)^2+H^4\right].
\]

Thus \(\lambda_\rho=\rho^{-1/2}\) and
\(\lambda_R=(R_{\alpha\beta\gamma\delta}R^{\alpha\beta\gamma\delta})^{-1/4}\)
are monotone on each branch. Collapse and expansion are compared at equal
clock value after applying \(H_-\mapsto-H_+\). The scored vector is

\[
\left(\Omega,\Sigma^2,\frac{R_{\mu\nu}R^{\mu\nu}}{36H^4},
w_{\mathrm{eff}},\lambda(H_-+H_+)\right),
\]

with unit weights and a pointwise Euclidean norm. The exact matched solution
has zero distance; the reported nonzero fine residuals measure interpolation
and integration error in the two-clock numerical control.

### A.3 Adverse and anisotropic controls

The Schwarzschild interior written as a Kantowski--Sachs cosmology approaches
the vacuum Kasner exponents \((-1/3,2/3,2/3)\). It therefore has nonzero
normalized shear and unit Weyl fraction, whereas spatially flat stiff FLRW has
zero shear and zero Weyl fraction. Its Ricci tensor is also vacuum rather than
scalar sourced. The comparison is rejected before any reconstruction score is
formed, providing an exact invariant adverse control.

For the homogeneous anisotropic Einstein--scalar control, let \(N=\ln a\),
\(\Sigma\) be the signed normalized shear amplitude, and
\(\Omega=1-\Sigma^2\). The autonomous equation is

\[
\Sigma'=(q-2)\Sigma,
\qquad
q=2\Sigma^2+\frac12(1+3w)\Omega,
\]

where a prime denotes \(d/dN\). For a massless scalar, \(w=1\), hence
\(q=2\) and \(\Sigma'=0\): normalized shear labels conserved boundary
strata. The control suite uses \(\Sigma_0^2=10^{-4},10^{-2},10^{-1}\) over
40 e-folds toward the boundary. A radiation comparison, \(w=1/3\), evolves
toward shear domination and supplies a matter-model contrast. These tests
justify treating \(\Sigma^2\) as boundary data rather than as a removable
time-orientation convention.

### A.4 Reproducibility role

The homogeneous computations evolve both a redundant free system and the
constraint-reduced system. The latter supplies the scored trajectory; the
former is retained only as a constraint-drift diagnostic. This separation
prevents numerical violation of the Friedmann constraint from being counted
as physical failure of collapse--expansion equivalence. All clock matches are
interpolated on their common invariant domains, and the same parity and norm
definitions are applied to the positive, adverse, and shear-stratified tests.

## Acknowledgments

This research was conducted independently, without institutional or grant
support. R.B. thanks Hazel for her tireless support and patience throughout
this research—and for saying yes on June 20th, in the middle of all of it.

## Data and code availability

The specifications, source code, machine-readable outputs, and independent
result-chain verifier supporting this study are publicly available at
https://github.com/RonBibb/line-e-boundary-state-paper together with the linked
computational repositories and machine-readable result records.

## References

1. L. Andersson and A. D. Rendall, “Quiescent cosmological singularities,”
   *Commun. Math. Phys.* **218**, 479--511 (2001), arXiv:gr-qc/0001047.
2. S. Kichenassamy and A. D. Rendall, “Analytic description of singularities
   in Gowdy spacetimes,” *Class. Quantum Grav.* **15**, 1339--1355 (1998).
3. E. Ames, F. Beyer, J. Isenberg, and P. G. LeFloch, “A class of solutions to
   the Einstein equations with AVTD behavior in generalized wave gauges,”
   *Ann. Henri Poincaré* **18**, 2295--2344 (2017), arXiv:1602.03018.
4. E. Ames, F. Beyer, and J. Isenberg, “Contracting asymptotics of the
   linearized lapse-scalar field sub-system of the Einstein-scalar field
   equations,” arXiv:1904.02854.
5. F. Beyer and P. G. LeFloch, “Second-order hyperbolic Fuchsian systems.
   Gowdy spacetimes and the Fuchsian numerical algorithm,”
   *Class. Quantum Grav.* **27**, 245012 (2010), arXiv:1006.2525.
6. F. Beyer and P. G. LeFloch, “Second-order hyperbolic Fuchsian systems and
   applications,” arXiv:1011.0270.
7. A. Franco-Grisales, “Developments of initial data on big bang singularities
   for the Einstein--nonlinear-scalar field equations,” arXiv:2409.17065.
8. B. Le Floch, P. G. LeFloch, and G. Veneziano, “Cyclic spacetimes through
   singularity scattering maps. The laws of quiescent bounces,”
   arXiv:2005.11324.
9. P. G. LeFloch, “Singularity scattering laws for bouncing cosmologies: a
   brief overview,” arXiv:2106.07958.
10. L. Andersson, H. van Elst, and C. Uggla, “Gowdy phenomenology in
    scale-invariant variables,” *Class. Quantum Grav.* **21**, S29--S57
    (2004), arXiv:gr-qc/0310127.
11. A. D. Rendall and M. Weaver, “Manufacture of Gowdy spacetimes with
    spikes,” *Class. Quantum Grav.* **18**, 2959--2975 (2001),
    arXiv:gr-qc/0103102.
12. H. Ringström, “Asymptotic expansions close to the singularity in Gowdy
    spacetimes,” *Class. Quantum Grav.* **21**, S305--S322 (2004),
    arXiv:gr-qc/0303051.
13. H. Ringström, “Strong cosmic censorship in \(T^3\)-Gowdy spacetimes,”
    *Ann. Math.* **170**, 1181--1240 (2009),
    doi:10.4007/annals.2009.170.1181.
14. P. T. Chruściel, J. Isenberg, and V. Moncrief, “Strong cosmic censorship
    in polarized Gowdy spacetimes,” *Class. Quantum Grav.* **7**, 1671--1680
    (1990).
15. M. Boyle and A. H. Mroué, “Extrapolating gravitational-wave data from
    numerical simulations,” *Phys. Rev. D* **80**, 124045 (2009),
    arXiv:0905.3177.
16. N. T. Bishop, R. Gómez, L. Lehner, and J. Winicour,
    “Cauchy-characteristic extraction in numerical relativity,”
    *Phys. Rev. D* **56**, 6298--6309 (1997), arXiv:gr-qc/9705033.
17. C. Reisswig, N. T. Bishop, D. Pollney, and B. Szilágyi,
    “Characteristic extraction in numerical relativity: binary black hole
    merger waveforms at null infinity,” *Class. Quantum Grav.* **27**, 075014
    (2010), arXiv:0912.1285.
18. W. C. Lim, “New explicit spike solution---non-local component of the
    billiard attractor,” *Class. Quantum Grav.* **25**, 045014 (2008),
    arXiv:0710.0628.
19. B. K. Berger and V. Moncrief, “Numerical investigation of cosmological
    singularities,” *Phys. Rev. D* **48**, 4676--4687 (1993),
    arXiv:gr-qc/9307032.
20. D. Garfinkle and M. Weaver, “High velocity spikes in Gowdy spacetimes,”
    *Phys. Rev. D* **67**, 124009 (2003), arXiv:gr-qc/0303017.

## Draft completion ledger

- [ ] Confirm author list, affiliations, corresponding author, and ORCID IDs.
- [x] Derive the reduced equations from the normalized action in an appendix.
- [x] Freeze the exact definitions of both geometric clocks and the scoring
      norm in the main text.
- [x] Insert production-family initial data and integrator details.
- [x] Generate the extraction--reconstruction flow figure.
- [x] Generate convergence, ablation, and exact-benchmark figures.
- [x] Add Appendix A with reproducible homogeneous-control definitions and
      analytic derivations.
- [ ] Verify every bibliography entry against the published version.
- [ ] Add a reproducibility manifest and archived repository DOI.
- [ ] Convert to the selected journal's LaTeX class after scientific review of
      this draft.
