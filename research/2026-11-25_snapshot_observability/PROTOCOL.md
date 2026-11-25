# Investigation scope

Started November 25, 2026 at 08:33 Eastern. Four-hour cap ends at 12:33 Eastern. The August submission and recorded experiments are source material. The proposed research article lives in this directory. Numerical work runs in the Linux guest, using a read-only import of the original solvers.

## Analytic question

For a marked periodic polarized Gowdy field, determine what a pair of finite-areal-time displacement snapshots stably determines. Test the conjectured exact dependence on the target slice and recoverability of the limiting velocity. Compare with Vest's loaded string, the wave snapshot literature, Li's singular scattering theorem, and classic Bessel cross-products. Require a proof and an explicit adverse sequence. Numerical checks test the formulas, including the zero mode and sharp lower-bound witnesses. An unrelated solver should check the Bessel transfer against the ODE.

## Geometric repair

Keep the stored Gowdy potential with positive quadratures. Derive the metric and physical stress tensor in that convention. Require the homogeneous vacuum and stiff-scalar constraints, explicit test of an independent late-window metric anchor, and comparisons on the original smooth grids and depths. Distinguish the analytic-input Bessel arm from endpoint-assisted measurements, and check an actual Bessel holdout separately. Report full scale only in a marked torus and limit numerical claims to the scored quantities.

## Nonlinear check

Release two P samples, endpoint Q and endpoint Q velocity to each decoder. These data are equivalent to the existing four-function finite-depth record. Run a trailing-slope baseline and a finite-window shooting iteration using the same inputs. Truth at the endpoint and earlier evaluation slices is reserved for scoring. A privileged Cauchy restart measures the solver floor. Use the legacy numerical tape and independently evaluate the inverse with fourth-order spatial differences at matched finite resolution. Report all of the fixed full-circle score and the original detector mask. Compare at depths 12 and 18 on 1024, 2048, and 4096 grids. An exact periodically transformed Bessel snapshot pair checks the method without using the forward discretization to manufacture data. A no-convergence or failure remains in the record. These tests concern the specified finite windows and IVP. Recovering a generic asymptotic spike state requires a separate result.

## Additional nonlinear falsification

Use a distinct off-geodesic periodic vacuum IVP with $P=0$, $P_\tau=2+1.4\cos\theta$, $Q=0.9\cos\theta+0.1\sin2\theta$, and $Q_\tau=0.025\cos3\theta$. Orthogonality makes the integrated initial momentum source zero. Evolve with the fourth-order stencil and reconstruct with the historical second-order stencil at depth 12 and resolutions 512 and 1024. At 1024 add an absolute perturbation $10^{-7}\cos(k\theta)$ to the earlier field snapshot only, for $k=1$ and $k=128$. Report both convergence residual and the additional target and inferred velocity errors. This probes precision dependence and independent field generation for a second IVP. It supplies no universal robustness estimate.

The initial auxiliary metric in the first distinct-IVP trial was zero despite a nonzero pointwise momentum source. Its maximum residual exceeded the historical $10^{-3}$ threshold. The result is preserved in the preconstraint record. Initialize the corrected metric by a periodic primitive of $2Q_\tau Q_\theta$ and anchor it at zero. The wave-map evolution and score are independent of the metric potential. The rerun checks that the independence holds in the results and that the reference passes the historical residual threshold.
