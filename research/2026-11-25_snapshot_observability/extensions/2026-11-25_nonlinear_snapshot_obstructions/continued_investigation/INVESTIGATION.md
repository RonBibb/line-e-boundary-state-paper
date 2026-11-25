# Continued four-hour Gowdy investigation

The original allocation runs from 08:33 to 12:33 Eastern on 25 November 2026. The continuation began at 10:37 Eastern. Earlier proofs and computations are preserved in the parent directory. This pass tests whether finite-time observations can robustly recover endpoint velocity in fully inhomogeneous Gowdy evolution and whether doing so requires new timing.

## Investigation tracks

1. Prove a short-window nonlinear stability theorem for frequency-assigned displacement measurements. Derive a sampling schedule with a positive minimum delay, a finite-band tail estimate, a decoder and a period-count lower bound in the matching observation norm. Test against arbitrary inhomogeneous coefficients under the stated Sobolev bounds and against constant-background arithmetic obstructions.
2. Check the theorem and observation design on Linux. Use independent forward data and inverse grids or stencils. Use exact rank-two collisions with matched original data, then perturbed rank-two histories and archive-derived data. Compare actual competing observation placements, timing perturbations and noise. Record failed cases and quantitative evidence.
3. Analyze the possible stronger statement with a fixed finite set of full snapshots. Compute the high-frequency observation symbol around moving inhomogeneous backgrounds. Identify a condition that controls nearly blind modes or show why a uniform energy estimate still fails. Separate a linearized result from any nonlinear inference.
4. Compare the specific outcomes to prior observability, data-assimilation and dynamical-sampling work. Review proof details and consolidate the mature additions into the manuscript. Include remaining open questions with precise hypotheses.

Scientific runs use the configured Linux guest. The earlier source paper and release are outside this revision.

The first proof pass found the bandwidth-independent nonlinear Duhamel estimate and a fixed three-slice continuum estimate at the explicit rank-two backgrounds. The latter uses the two null-ray transports and a finite rotation cover to control high frequencies. A second terminal reading excludes a possible finite-dimensional kernel. Both results require independent checks and numerical stress tests.


## Completed results and checks

- Relative null transport gives stable continuum recovery near inhomogeneous backgrounds. The regular-time proofs cover one or two withheld velocities. The metric gradient yields necessary and sufficient conditions for the high-frequency estimate at both rational and irrational gaps. The circular-loop reverse estimate bounds all frequencies throughout an explicit interval. A separate open low-velocity class gives stability for five asymptotic fields.
- Exact common-node histories show full-geometry ambiguities with different limiting speeds, and with identical speeds but different limiting metric scales. The earlier nonlinear proof's metric determinant sign was corrected and its PDF rebuilt.
- The archive calculation bounds the exact centered-grid derivative over a ball and under simultaneous changes of recorded fields. The independent production calculation explains the historical depth-gap discrepancy and shows second-order error versus a continuum reference.
- All 23 checks in the final Linux replay returned valid evidence. A separate Linux check illustrates the longer-delay weak modes with both metric profiles supplied. The final manuscript was built and visually checked. The working result is ready for independent mathematical and priority review.

[Working manuscript](output/pdf/Gowdy_Finite_Time_Observability.pdf) · [Replay summary](REPLAY_SUMMARY.json) · [Correction ledger](CORRECTION_LEDGER.md)
