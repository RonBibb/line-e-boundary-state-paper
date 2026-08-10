# E1 Reconstruction-Adjudication Test Specification

**Date:** 2026-08-09  
**Test ID:** `E1-RA1`  
**Status:** frozen before execution  
**Purpose:** adjudicate the held CBE-2I reconstruction and completeness claims
without modifying the original frozen repository or result

## 1. Questions

1. Does an independent implementation reproduce the original CBE-2I baseline?
2. Is the imported collapse-side \(\lambda\) anchor load bearing?
3. Does reconstruction close when the extraction window withholds the endpoint
   fields rather than reproducing them algebraically?
4. Does the frozen constraint gate pass under both the specified and
   implemented normalizations?
5. Does an explicit future-normal/time-parity construction reproduce the
   original direct comparison?

## 2. Unchanged model

Use the equations, production initial data, centered spatial differences,
RK4 integrator, time-step factor, comparison surfaces, and observable formulas
from the frozen CBE-2I implementation. Import those functions read-only.

No original source or result file may be changed.

## 3. Test grid

\[
N\in\{256,512,1024\},
\qquad
\tau_f\in\{16,24\}.
\]

Comparison surfaces remain

\[
\tau_c\in\{2,4,8,12\}.
\]

## 4. Extraction variants

### V0 — original endpoint-window extraction

\[
v_f=\frac{f(\tau_f)-f(\tau_f-2)}{2},
\qquad
f_\infty=f(\tau_f)-v_f\tau_f.
\]

This reproduces \(f(\tau_f)\) exactly.

### V1 — predictive holdout extraction

Withhold all field values at \(\tau_f\) from extraction and use

\[
v_f^{\rm hold}
=\frac{f(\tau_f-2)-f(\tau_f-4)}{2},
\qquad
f_\infty^{\rm hold}
=f(\tau_f-2)-v_f^{\rm hold}(\tau_f-2).
\]

Initialize the reconstruction at \(\tau_f\) by extrapolation:

\[
f_{\rm rec}(\tau_f)
=v_f^{\rm hold}\tau_f+f_\infty^{\rm hold}.
\]

The endpoint field error is recorded before backward evolution.

## 5. Anchor variants

### A0 — imported collapse anchor

Use the original prescription

\[
\lambda_{\rm rec}(\tau_f,0)=\lambda_{\rm col}(\tau_f,0).
\]

### A1 — source-independent fixed anchor

Set

\[
\lambda_{\rm rec}(\tau_f,0)=0
\]

without reading the collapse-side endpoint anchor. This is an ablation, not a
claim that zero is the physically correct boundary normalization.

The comparison tests whether the omitted global mode is load bearing for the
currently scored observables.

## 6. Arms

| Arm | Extraction | Anchor | Purpose |
|---|---|---|---|
| `V0-A0` | original endpoint window | collapse anchor | baseline reproduction |
| `V0-A1` | original endpoint window | fixed zero | anchor ablation |
| `V1-A0` | predictive holdout | collapse anchor | endpoint-copy ablation |
| `V1-A1` | predictive holdout | fixed zero | strictest combined ablation |

## 7. Scoring

For each arm, record:

- maximum scored distance over all comparison surfaces;
- component-wise sup residuals;
- pointwise relative \(\lambda_\rho\) gap;
- endpoint prediction errors in \(P,\phi,P_{,\tau},\phi_{,\tau}\);
- specified and implemented momentum-residual normalizations; and
- explicit parity residual.

The future-directed normal on the expansion branch is defined as the negative
of the increasing-\(\tau\) normal. Therefore its expansion scalar is

\[
\Theta_{+}^{\rm future}=-\Theta_{+}^{(+\tau)}.
\]

The time-reversal map reverses this time-odd quantity once more before
comparison. The parity-transformed result must agree with the original direct
common-normal comparison to roundoff.

## 8. Frozen readings

1. **Baseline reproduction:** `V0-A0` must reproduce the original headline
   distances to relative tolerance \(5\%\) or absolute tolerance \(10^{-8}\),
   whichever is larger.
2. **Anchor independent:** if `V0-A1` remains below \(10^{-3}\), the anchor is
   not load bearing for the frozen score on this test. Otherwise the
   four-function completeness claim fails unless the anchor is proved gauge.
3. **Predictive reconstruction:** `V1-A0` must decrease with resolution and
   remain below \(10^{-3}\) at \(N=1024\) for both depths to support a stronger
   non-endpoint-copy reconstruction claim.
4. **Strict four-function arm:** `V1-A1` passing the same gate supports, but
   does not prove, source-independent four-function completeness. Failure
   identifies the missing global mode and/or finite-depth remainder as load
   bearing.
5. **Constraint discrepancy:** both residual definitions are reported. A
   classification change under the specified definition is a formal frozen-gate
   failure.
6. **Parity closure:** explicit parity and common-normal scores must agree to
   \(10^{-12}\) absolute.

No failed ablation invalidates the original numerical evolution. It narrows
the interpretation of what information the reconstruction actually uses.

