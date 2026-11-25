# A nonlinear extension for the Gowdy paper

**Position the successor around finite-time observability.** The original archive can support accurate reconstruction from its four observations. The mathematical question is where the same observation design loses information. The extension gives sharp criteria in tractable classes and exact indistinguishable histories in the full vacuum equations, including solutions with two independent target directions.

## Results and scope

| Result | What it establishes |
|---|---|
| Sharp homogeneous threshold | The four-field record \(P(a),Q(a),Q_t(a),P(b)\), with one metric anchor, determines the missing \(P_t(a)\) across the speed-bounded homogeneous class exactly when \(V\log(b/a)\le 2.399357\ldots\). Above it there are exact two-branch collisions. A terminal \(Q\) or \(\lambda\) resolves the homogeneous fold. |
| Exact nonlinear examples | At resonant times, subcritical inhomogeneous histories can share both full induced spatial metrics, \(Q_t\) at both ends, and \(\lambda\) at every intermediate time. The proof uses geodesic-image solutions. A second example lifts the partial-record fold with strictly positive subcritical asymptotic speeds everywhere. |
| Two independent field directions | A periodic Gowdy-to-Ernst lift carries the full-metric collision outside the geodesic-image class. Each target map has rank two on open sets. The measured \(Q\) is spatially nonconstant at both ends. One two-mode construction excludes equivalence by time-orientation-preserving, torus-equivariant boundary isometries even after changing the Killing basis. Within this family, the full-metric observation has exactly two branches when the inferred parameter magnitude is positive. One off-node scalar reading distinguishes them. Transformed examples can reach asymptotic speed one where the seed logarithmic velocity vanishes. |
| A constructive third snapshot | With constant initial target position in the geodesic-image class, full orbit snapshots are injective exactly when they have no common Dirichlet Fourier zero. An explicit rational areal gap gives one derivative of sharp loss in the energy-normalized velocity using the initial snapshot and the added one. |
| Moving-background stability | At speed \(1/\sqrt2\), opposite Jacobi phase corrections give a sharp three-derivative loss on a resonant subsequence in the energy-normalized velocity. Nearby backgrounds have exact linearized zeros. Adding \(Q(b)\) or \(\lambda(b)\) reduces the loss to one on the displayed background. At speed \(1/2\), even both induced metrics have an exact linearized kernel. Smooth signed pairs there separate cubically, excluding any pairwise inverse Hölder exponent above \(1/3\). |
| Historical test window | The centered semidiscrete tangent of the normalized two-snapshot observation is within about \(0.06201\) of the identity in maximum row norm over three sampled resolutions. Four-data reconstruction reduced the observed masked error from \(0.82215\) to \(0.00298\) on the same tape. The tangent check supports a local finite-grid explanation. A rigorous continuum neighborhood estimate remains open. |

The exact claims follow from formulas and proofs in the working draft. The Linux computations independently check the geodesic equations, leading Fourier coefficients, cubic response, and finite-grid archive behavior. The transformed rank-two witness has equal recorded fields within \(9\times10^{-16}\) and differs by \(0.1028\) in interior \(P\) at its sampled slice. Those values illustrate the exact proof. The finite-grid row bounds use sampled generator coefficients.

## What remains for priority review

- Two-time wave observability, the arithmetic of sample gaps, and the use of extra snapshots precede this work. Vest treats free and loaded strings. Christensen and coauthors prove uniqueness and compatibility results for two and three free-wave snapshots. The single scalar mode with nodal observations belongs to this tradition. [Vest](https://www.math.u-szeged.hu/ejqtde/p2411.pdf) · [Christensen et al.](https://arxiv.org/pdf/2308.12208)
- The exact Gowdy building blocks are also established. Berger and Moncrief give the pseudo-unpolarized Bessel construction in §§III–IV, including Eq. (3.19). Lim gives the homogeneous rotated Kasner solutions in §4.2. Rendall and Weaver give the Gowdy-to-Ernst transformation in §2 and discuss its periodicity obstruction. The contribution to check is how two different branches match the specified measurements, including the strengthened rank-two pair and its boundary-isometry test. [Berger and Moncrief](https://arxiv.org/pdf/gr-qc/9307032) · [Lim](https://arxiv.org/pdf/0710.0628) · [Rendall and Weaver](https://arxiv.org/pdf/gr-qc/0103102)
- Amorim, Bernardi and LeFloch already compute past and future Gowdy evolution, including limitations around steep spike data. This extension addresses what initial velocity can be inferred when it was omitted. The specific candidates for independent priority and proof review are the sharp speed–window fold, matching both induced metrics and the measured velocity component in a rank-two pair, and the cancellation of the two Jacobi channels with its sharp rate. [Bidirectional Gowdy computation](https://arxiv.org/pdf/0811.1266)

## Where to take the paper

Develop a standalone observability successor. Open with the rank-two counterexamples and the sharp homogeneous threshold. Present the original archive and corrected benchmark as the motivating example. Use the short third snapshot as a positive result, stating the geodesic-image hypothesis. Keep the polarized sharp target-slice and storage analysis in the companion manuscript. The two texts can be merged later if readers see one coherent paper.

The next review should pair a mathematical relativist with a PDE observability reader. The most consequential checks are the periodicity and mean-curvature argument for the two-mode lift, and the sign and cubic order of the mixed Jacobi coefficient. The next research question is whether a suitably timed third slice gives a quantitative inverse near a genuinely inhomogeneous rank-two background. The draft proves a positive estimate in the geodesic-image class and identifies a single off-node reading that separates the two explicit rank-two branches.

## Reproduction

Run the following in this directory on Linux with NumPy 2.2.6, SciPy 1.15.3, and mpmath 1.3.0. Each program prints JSON. The archive program uses the frozen `bb-cbe3b` traces in the repository’s `reproduction/inputs/` directory.

```sh
python -B check_obstruction.py > obstruction_checks.json
python -B check_exact_collision.py > exact_collision_checks.json
python -B check_rank_two.py > rank_two_checks.json
python -B check_cubic.py > cubic_checks.json
python -B check_archive_tangent.py > archive_tangent_checks.json
```

Build `nonlinear_snapshot_extension.tex` with pdfLaTeX from this directory. The proof sections it imports are stored alongside it.

[Working draft](output/pdf/Gowdy_Nonlinear_Snapshot_Obstructions.pdf) · [Earlier corrected manuscript](../../output/pdf/Gowdy_Two_Snapshot_Recovery.pdf) · [Corrections and impact](../../CORRECTIONS.md)

## Four-hour continuation

The continuation adds full-metric recovery at irrational gaps, a finite-grid certificate on the archived restart, and a production-side check of the historical smooth modulus. It also corrects the sign of the metric exponential in this working draft's boundary-isometry proof. [Continuation and proposed paper](continued_investigation/READOUT.md)
