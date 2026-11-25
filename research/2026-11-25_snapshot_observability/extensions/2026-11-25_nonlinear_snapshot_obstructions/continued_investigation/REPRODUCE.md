# Reproducing the continuation

Scientific computation used the Linux aarch64 guest with NumPy 2.2.6, SciPy 1.15.3 and mpmath 1.3.0. Run in this directory within the checkout. The enclosing checkout can be named or placed anywhere. The archive certifications use the exact binary values of the depth-18 data in the bundled CBE3b trace. The smooth benchmark loads the bundled historical CBE2i Python source. See the [frozen inputs and environment](../../../../../reproduction/README.md).

| Evidence | Driver |
|---|---|
| Full nonlinear delayed observations and stencil-sensitive reconstruction | `check_nonlinear_sampling.py` |
| Direct null transport, high modes, and fixed-slice modal singular values | `check_null_transport.py` |
| Open two-slice examples, polarized metric gradient, and analytic background path | `check_open_two_snapshot.py` |
| Fixed-slice nonlinear reconstruction and measured-field noise | `reconstruct_fixed_snapshots.py` |
| Numerical dispersion correction in mild and strongly coupled cases | `check_dispersion_correction.py` |
| Both missing velocity components in the circular-loop family | `check_both_velocity_holonomy.py` |
| Both missing velocities with the full metric at an irrational gap | `check_full_metric_holonomy.py` |
| Longer circular-loop delay with weak modes even for two full metrics | `check_long_delay_metric.py` |
| Continuation from the recovered class to the low-velocity asymptotic state | `check_low_velocity_tail.py` |
| Modal sidebands and second-order holonomy in the exact Ernst pair | `check_ernst_both_velocities.py` |
| Early-anchor production metric modulus, Bessel tail and continuum reference | `check_production_modulus.py` |
| Circular-loop detuning and limiting subcritical speed | `check_circular_loop.py` |
| Interval brackets for the fundamental and sixth-mode phase | `certify_asymmetry.py` |
| Phase and metric asymptotics, and mean-metric boundary identity | `check_asymptotic_geometry.py` |
| Sharp modulus multiplier probes, amplitude scaling and the twisted-bundle example | `check_extended_predictions.py` |
| Archive velocity-ball enclosure and perturbation of stored measurements | `certify_archive_tangent.py`, `certify_perturbed_archive_record.py` |
| Secant-based archive iteration with and without noise | `replay_archive_contraction.py` |

Run `OPENBLAS_NUM_THREADS=1 python -B run_replay.py` to execute the new checks and the selected parent checks. The script writes Linux JSON evidence and a summary under `build/replay`. [Completed replay](REPLAY_SUMMARY.json)

The longer-delay check is a separate supplemental run. Its floating-point timing illustrates the closed-transport mechanism. Irrationality and the loss of an energy bound for the specified class use the proof. [Supplemental evidence](long_delay_metric_checks.json)

The drivers share helpers within the research tree. The tests produce JSON on stdout. To regenerate scientific tables, run `build_figure_data.py` after the associated evidence is generated. The LaTeX root is `finite_time_observability.tex`, compiled from this directory with pdfLaTeX and PGFPlots. It inputs earlier sections from the parent directory. All main artifacts and figure tables are indexed in the integrity manifest.

Observable comparisons and numerical conditioning trials are illustrations. The archive enclosures, annulus root brackets and explicitly interval-bounded comparisons have separate certificates. The full-vacuum continuum recovery estimates follow from the proof text.

[Working manuscript](output/pdf/Gowdy_Finite_Time_Observability.pdf) · [Readout](READOUT.md)
