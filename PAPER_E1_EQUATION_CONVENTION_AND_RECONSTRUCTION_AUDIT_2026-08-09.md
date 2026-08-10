# Paper E1 Equation, Convention, and Reconstruction Audit

**Date:** 2026-08-09  
**Status:** completed initial audit; manuscript admission conditional on defects below  
**Sources audited:** CBE-1 result and implementation; CBE-2I frozen specification,
result, and implementation; Paper-E1 transfer and scope freeze

## 1. Governing convention sheet

### 1.1 Einstein–scalar normalization

The declared units are

\[
8\pi G=1.
\]

The matter model is a real massless scalar with stress tensor

\[
T_{\mu\nu}
=\nabla_\mu\phi\nabla_\nu\phi
-\frac12 g_{\mu\nu}(\nabla\phi)^2,
\qquad
\Box_g\phi=0.
\]

The CBE-2I scalar is scaled so that its reduced Gowdy evolution and metric
quadratures enter symmetrically with the polarized gravitational variable.
The manuscript must derive this scaling from the action rather than rely on
the implementation convention alone.

### 1.2 Homogeneous definitions

For mean expansion rate \(H\), shear tensor \(\sigma_{\mu\nu}\), and scalar
energy density \(\rho\), the frozen normalized variables are

\[
\Omega=\frac{\rho}{3H^2},
\qquad
\Sigma^2=\frac{\sigma_{\mu\nu}\sigma^{\mu\nu}}{6H^2}.
\]

The homogeneous clocks are

\[
\lambda_\rho=\rho^{-1/2},
\qquad
\lambda_{\mathcal R}=K^{-1/4},
\]

where \(K=R_{\alpha\beta\gamma\delta}R^{\alpha\beta\gamma\delta}\) is the
Kretschmann scalar. The symbol \(\lambda_{\mathcal R}\) should be replaced in
the manuscript by \(\lambda_K\) to avoid confusion with the Ricci scalar.

For flat stiff FLRW,

\[
H=\pm\frac{1}{3|t|},
\qquad
\rho=3H^2,
\qquad
\frac{R_{\mu\nu}R^{\mu\nu}}{H^4}=36.
\]

The implementation formulas for the Ricci-squared and Kretschmann invariants
are consistent with the flat-FLRW curvature identities under the declared
sign-independent squared contractions.

### 1.3 Polarized Gowdy gauge

The declared metric is

\[
ds^2=e^{(\lambda+\tau)/2}
\left(-e^{-2\tau}d\tau^2+d\theta^2\right)
+e^{-\tau}\left(e^Pdx^2+e^{-P}dy^2\right),
\]

with periodic \(\theta\) and collapse toward decreasing orbit area as
\(\tau\to+\infty\).

The lapse for increasing \(\tau\) is

\[
N=e^{(\lambda+\tau)/4-\tau}.
\]

The reduced equations used in the implementation are

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
=2\left(P_{,\tau}P_{,\theta}
+\phi_{,\tau}\phi_{,\theta}\right).
\]

Direct differentiation of the three spatial scale factors gives

\[
H_\theta=\frac{\lambda_{,\tau}+1}{4N},
\qquad
H_x=\frac{P_{,\tau}-1}{2N},
\qquad
H_y=\frac{-P_{,\tau}-1}{2N}.
\]

These agree with the implementation for the normal orientation associated
with increasing \(\tau\). The mean expansion is

\[
H=\frac{H_\theta+H_x+H_y}{3}
=\frac{\lambda_{,\tau}-3}{12N}.
\]

For the production regime this orientation is normally contracting. Any
expanding interpretation therefore requires an explicit reversal of the
unit normal; it cannot be inferred merely from integrating the same equations
toward decreasing \(\tau\).

### 1.4 Polarized comparison quantities

The code evaluates

\[
\rho=\frac12 e^{-(\lambda+\tau)/2}
\left(e^{2\tau}\phi_{,\tau}^2+\phi_{,\theta}^2\right),
\]

\[
\Omega=\frac{\rho}{3H^2},
\qquad
\Sigma^2=\frac{\sum_i(H_i-H)^2}{6H^2},
\qquad
s_P=\frac{H_x-H_y}{\Theta},
\]

and

\[
\lambda_A=e^{-\tau/2},
\qquad
\lambda_\rho=\rho^{-1/2}.
\]

The formulas agree algebraically with the declared metric. Their status is
more accurately described as **geometrically defined relative to the Gowdy
symmetry and areal foliation**, not as unrestricted four-dimensional scalar
invariance:

- \(\rho\), \(H\), and \(\Sigma^2\) depend on the geometrically selected
  unit normal;
- \(\lambda_A\) depends on normalization of the symmetry-orbit coordinates;
- \(s_P\) changes sign under interchange of the two Killing directions; and
- \(\Theta\lambda_A\) is time odd and requires an explicit parity rule.

The manuscript should use “symmetry- and foliation-geometric observables”
unless the residual gauge group is fixed and invariance under it is proved.

## 2. Boundary data versus comparison observables

The paper must distinguish:

\[
\mathfrak B
=\text{constrained singular/asymptotic data},
\]

from

\[
\mathfrak O
=\text{derived observables used to score reconstruction}.
\]

For the polarized implementation, the advertised reduced data are

\[
\mathfrak B_{\rm advertised}
=\{v_P(\theta),v_\phi(\theta),P_\infty(\theta),\phi_\infty(\theta)\}.
\]

The scored observable vector is

\[
\mathfrak O
=\{\Omega,\Sigma^2,s_P,\Theta\lambda_A,\lambda_\rho\},
\]

with \(\lambda_\rho\) used as a clock-consistency diagnostic rather than as
a component of the Euclidean distance.

Constraint residuals are numerical-health diagnostics, not physical boundary
coordinates and not components of \(\mathfrak O\).

## 3. Reconstruction audit

### 3.1 Implemented extraction

At extraction depth \(\tau_f\), for \(f\in\{P,\phi\}\), the code defines

\[
v_f=\frac{f(\tau_f)-f(\tau_f-2)}{2},
\qquad
f_\infty=f(\tau_f)-v_f\tau_f.
\]

The reconstructed endpoint fields are then

\[
f_{\rm rec}(\tau_f)=v_f\tau_f+f_\infty=f(\tau_f),
\qquad
\partial_\tau f_{\rm rec}(\tau_f)=v_f.
\]

Thus the reconstruction copies the finite-depth field values \(P(\tau_f)\)
and \(\phi(\tau_f)\) exactly by algebraic construction, while replacing the
instantaneous endpoint velocities by two-surface averaged velocities. The
existing phrase “no point-copy of finite-\(\tau\) fields” is false as written
and must be retired.

This does not make the test empty: discarded decaying-mode and velocity
information can still prevent a successful return. It does, however, mean
that the test is presently a **lossy endpoint compression and backward
reconstruction test**, not reconstruction from an independently prescribed
singular state in the stronger Fuchsian sense.

### 3.2 Hidden metric anchor

After computing the reconstructed momentum integrand, the code obtains a
periodic antiderivative for \(\lambda\) and then executes

```text
lam0 += col_store[depth][4][0]
```

Therefore the reconstruction imports the collapse-side scalar anchor

\[
c_\lambda=\lambda(\tau_f,\theta=0).
\]

The actual reconstruction data are at least

\[
\mathfrak B_{\rm implemented}
=\{v_P,v_\phi,P_\infty,\phi_\infty,c_\lambda\}.
\]

Whether \(c_\lambda\) is removable by a declared global coordinate
normalization or represents a physical scale datum must be derived. Until
then, the claim that the four functions alone completely characterize the
scored state is not established.

### 3.3 Constraint and periodic compatibility

The implementation enforces the momentum quadrature for \(\lambda\) and the
initial production data satisfy the zero-mean integrability condition by
parity. The paper-level boundary space must nevertheless state explicitly:

\[
\oint d\theta\,
\left(P_{,\tau}P_{,\theta}
+\phi_{,\tau}\phi_{,\theta}\right)=0,
\]

together with the condition determining the global \(\lambda\) mode. The
four functions are not unrestricted elements of a Cartesian function space;
they lie on an admissible constraint manifold.

## 4. Parity and orientation audit

The homogeneous implementation explicitly compares

\[
H_{\rm col}+H_{\rm exp},
\]

which correctly encodes odd time parity after the two branches are placed on
a common positive boundary-distance clock.

The polarized implementation evolves the reconstructed state from large
\(\tau\) toward smaller \(\tau\) but computes \(H_i\) using the same
increasing-\(\tau\) normal formula on both stored branches. It then directly
subtracts \(\Theta\lambda_A\) rather than applying a separately documented
parity transformation.

This can be internally consistent if the reconstructed branch is interpreted
as the same orientation-normalized geometric solution sampled by reverse
integration. It is not yet a demonstrated comparison between physically
opposite future-directed normals. Before manuscript use, the authors must
freeze one of two interpretations:

1. **representation test:** compare both solutions in the same areal-normal
   convention and avoid calling the second normal physically expanding; or
2. **contracting/expanding parity test:** define future-directed normals on
   each branch separately and apply the explicit sign map to all time-odd
   quantities.

The current result record conflates these interpretations.

## 5. Specification-to-code discrepancies

### D-EQ-1 — constraint normalization differs

The frozen specification defines the normalized momentum residual with a
denominator based on \(\sup(|\lambda_{,\theta}|+\epsilon)\). The code divides
by \(\sup|2(P_{,\tau}P_{,\theta}+\phi_{,\tau}\phi_{,\theta})|+10^{-30}\).
These coincide only in the continuum constraint-satisfying limit. The machine
gate did not implement the frozen expression literally.

**Disposition:** recompute the reported residuals under both definitions from
the preserved data or rerun only the bounded verifier if full fields are not
stored. Preserve both values.

### D-EQ-2 — advertised four-function completeness omits \(c_\lambda\)

The code imports a collapse-side \(\lambda\) anchor.

**Disposition:** determine its gauge/scale status and either include it as a
boundary datum or fix it identically on both branches by a prospective,
source-independent convention.

### D-EQ-3 — finite-depth fields are copied exactly

The extraction formulas guarantee equality of \(P\) and \(\phi\) at
\(\tau_f\).

**Disposition:** replace the independence claim with an exact information-flow
statement. If a stronger singular-data reconstruction claim is desired, add a
separately designed test based on an asymptotic remainder model; do not
reinterpret the completed run retrospectively.

### D-EQ-4 — polarized time parity is not explicit

The code directly compares a time-odd clock-rate quantity under a common
normal convention.

**Disposition:** derive the branch orientation and parity map before deciding
whether “contracting and expanding” is justified in the headline.

### D-EQ-5 — invariant terminology is too broad

Several quantities are geometric only after the symmetry, areal foliation,
Killing-direction labels, and normal orientation are fixed.

**Disposition:** prove invariance under the residual gauge group or narrow the
terminology.

### D-EQ-6 — boundary-depth and discretization errors are not separated

The campaign shows second-order spatial convergence at each of two depths,
but the fine residual increases modestly from depth 16 to 24 while remaining
below the frozen tolerance. The current classification treats this as a
numerical floor without a derived asymptotic truncation/conditioning model.

**Disposition:** construct an error budget separating resolution, extraction
depth, finite-window velocity estimation, and roundoff/conditioning. Existing
records may suffice; no long simulation is automatically authorized.

## 6. Homogeneous constrained-evolution audit

The constrained FLRW implementation replaces \(H\) algebraically using

\[
H=\pm\frac{|\dot\phi|}{\sqrt6},
\]

then evolves \((a,\dot\phi)\). Within spatially flat homogeneous GR plus a
massless scalar, this is an exact reduction of the Hamiltonian constraint and
scalar equation after the expansion/contraction sign is chosen.

The homogeneous numerical result is therefore a formulation control rather
than strong independent evidence for a new physical equivalence. Its proper
roles in E1 are:

- exact positive control;
- demonstration of the two-clock scoring implementation;
- demonstration of structural drift in the unreduced free formulation; and
- illustration of conserved-shear strata.

The paper should not make the homogeneous numerical agreement carry the
novelty claim.

## 7. Effect of imported prior-art structure

The competing literature suggests a stronger formal organization:

\[
\mathscr B_{\rm adm}
=\{\text{singular data satisfying boundary constraints}\},
\]

\[
\mathcal E:\mathscr C_-\to\mathscr B_{\rm adm},
\qquad
\mathcal R:\mathscr B_{\rm adm}\to\mathscr C_+,
\]

followed by a scored observable map

\[
\mathcal O:\mathscr C_\pm\to\mathscr V.
\]

E1 currently tests a finite-depth approximation to this chain, but its
implemented \(\mathcal E\) retains endpoint field values and a metric anchor.
The manuscript must describe that map exactly before using the language of
singular-data completeness.

No scope expansion is adopted by this audit. A new Fuchsian reconstruction,
geometric scattering-map classification, or additional numerical campaign
would require a documented scope amendment.

## 8. Step-3 disposition

**Classification:** `CONVENTION-CLOSED / RECONSTRUCTION-CLAIM-HELD`.

The governing equations and most metric-derived formulas are internally
consistent under the declared gauge. The following claims remain held:

1. completeness of the four-function boundary representation;
2. absence of finite-depth point copying;
3. fully invariant status of the entire scored state vector; and
4. physical contracting-versus-expanding parity in the polarized campaign.

Paper development may continue with the completed result treated as a
finite-depth compression/reconstruction experiment. Before a manuscript
claim is frozen, D-EQ-1 through D-EQ-6 must be adjudicated analytically or by
bounded verification.
