# E1 depth and norm-sensitivity result

**Test ID:** `E1-DEPTH-NORM-SENSITIVITY`  
**Machine record:** `E1_DEPTH_NORM_SENSITIVITY_RESULT.json`  
**Status:** completed

## Question

Does the predictive production audit remain stable under alternative scoring
weights, and does a shallower extraction expose asymptotic truncation that is
invisible at depths 16 and 24?

## Frozen extension

The production initial data, holdout width, integrator, asymptotic metric
route, and imported finite-depth anchor were unchanged. Extraction depths
were (8,12,16,24), resolutions were (N=256,512,1024), and every case was
scored on the common surfaces (	au=2,4). The observable differences were
((\delta\Omega,\delta\Sigma^2,\delta s_P,\delta(\Theta\lambda_A))).

Four norms were recorded: the manuscript identity-weight Euclidean norm; a
shape-emphasized Euclidean norm with weights ((2,2,2,1/2)); a
clock-emphasized norm with weights ((1/2,1/2,1/2,2)); and the maximum
component norm.

## Results

At (N=1024), the identity-weight maxima at depths (8,12,16,24) were

\[
2.58457\times10^{-3},\quad
1.76861\times10^{-6},\quad
8.69854\times10^{-7},\quad
8.69878\times10^{-7}.
\]

The shallow depth-8 case fails the (10^{-3}) gate under the identity,
shape-emphasized, and component-maximum scores. Its largest component is
(\|\delta\Omega\|_\infty=2.43070\times10^{-3}), and its endpoint field
prediction errors reach (2.20225\times10^{-4}). The audit therefore detects
measurable two-term truncation rather than automatically closing.

At depth 12, the fine-grid component maxima are

\[
(1.55630\times10^{-6},\ 8.63388\times10^{-7},\
1.87371\times10^{-7},\ 8.44130\times10^{-7})
\]

in the order ((\Omega,\Sigma^2,s_P,\Theta\lambda_A)). At depths 16 and 24,
the clock-rate component dominates the fine-grid identity score at
(8.70\times10^{-7}), while the field-shape components are substantially
smaller.

At depths 16 and 24 all four norm choices converge with difference ratios
within (0.1\%\) of four. Their fine-grid values remain below (1.74\times
10^{-6}) even when the clock component is doubled. The pass/fail conclusion
at the production depths is therefore not an artifact of identity weighting.

## Interpretation

This is an observable-relative sensitivity audit, not source-independent
recovery of the limiting metric modulus: every production row retains the
finite-depth metric anchor. It adds a deliberately harder shallow extraction
and shows a transition from truncation-dominated failure at depth 8 to
spatial-error-dominated closure at depths 16 and 24.
