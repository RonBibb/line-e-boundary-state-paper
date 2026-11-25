# Nonlinear snapshot extension

Started 25 November 2026 at 09:38 Eastern under the remaining original four-hour allocation. The observation map is the four finite fields from the original restart, written as $P(a),Q(a),Q_s(a),P(b)$ for $s=\log(t/a)$ and $0<a<b$. It may also contain one metric anchor at $a$. Continue with the known homogeneous hyperbolic geodesic, then linearize the full vacuum wave map around it. Preserve the earlier submission and successor artifacts.

## Analytic targets

1. Solve the homogeneous inverse with a bounded hyperbolic speed. Seek an exact injectivity criterion, a collision, and a smallest new datum that resolves the collision. Verify the criterion both by the Jacobi transfer and the explicit homogeneous ODE.
2. Derive the frequency-resolved mixed observation map directly from the wave-map linearization. Check the case in which the two channels' order-$k^{-1}$ contributions cancel, and distinguish a rate loss from an exact linearized kernel. Check which additional channels restore the estimate.
3. At an exact linearized kernel, quantify the obstruction to stability of the underlying nonlinear map. Supplement the proof with a converged nonlinear wave-map check on smooth periodic perturbed data. Preserve all comparisons, including detuned data if used.

## Numerical scope

Run only on the configured Linux guest. Use explicit geodesic formulas, independent coordinate ODE integration, and high-precision Bessel functions for small multipliers. For the nonlinear test use Fourier differentiation with three resolutions and two temporal tolerances. Use equal $P(a),Q(a),Q_s(a)$ and compare signed perturbations of $P_s(a)$ at a declared fixed areal gap. Measure the stronger kernel criterion over the full terminal orbit state, with the extra metric change checked at linear order and the combined terminal metric difference compared nonlinearly. Numerical checks test the analytic formulas and a bounded finite-time witness. They do not amount to a continuum analysis near a true Gowdy spike.

## Production-history check

Use the already archived centered-grid production traces in the sibling CBE3b directory, on the same initial-value problem as the earlier failed secant restart. Starting at the archived depth-18 fields, differentiate the reduced semidiscrete map from endpoint $P_\tau$ to the two-point slope between depths 16 and 18. Check the row-sum distance from the identity at three resolutions using periodic coloring, temporal refinement at the finest resolution, and finite differences in a fixed test direction. Compute the locality tail scale from sampled generator bounds and report it as a numerical diagnostic. It is not an interval enclosure of the archived orbit or the continuum derivative.

## Exact inhomogeneous collision and detuning

For the exact collision at log gap 3, choose the reference time using a zero of the scalar two-snapshot Bessel determinant. Lift a single invisible scalar mode through both homogeneous target geodesics. Check equality at the two observations, different interior values, and the uniform subcritical asymptotic-speed bound. Also vary the center of the equal-weight Jacobi background to check the asymptotic distance of modes with an exact first-order kernel.

## Extension beyond the geodesic-image class

Use the common geodesic current to apply the periodic Gowdy-to-Ernst transformation to a matched pair. Confirm equality of both induced metrics and the measured velocity component on the nodal slices. Check that the target differential has rank two on an open set. For a two-mode example use the ratio of the first two radial Dirichlet eigenvalues to obtain simultaneous nodes, and audit the possible boundary isometries using the quotient circle and mean curvature. Test the exact formulas in Linux. Add an explicit nonresonant observation placement in the geodesic-image class and verify the predicted high-frequency sensitivity.
