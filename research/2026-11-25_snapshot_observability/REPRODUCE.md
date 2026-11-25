# Reproduce the new checks

The source files in this directory replay the new measurements in a Linux Python environment with NumPy and SciPy. The investigation used Linux aarch64, NumPy 2.2.6, and SciPy 1.15.3. The drivers print JSON to standard output. The original solvers are imported read-only. The drivers load the frozen sources from the repository’s `reproduction/inputs/` directory.

| Driver | Consumes | Output |
|---|---|---|
| `check_snapshots.py` | analytic Gowdy and Euler–Poisson–Darboux modal formulas and an independent DOP853 integration | `snapshot_checks.json` |
| `check_geometry.py` | historical smooth Gowdy solver and exact Bessel construction, with revised geometric scoring | `geometry_checks.json` |
| `check_nonlinear.py` | historical vacuum solver, a fourth-order variant, analytic transformed periodic data, and local Lim formulas | `nonlinear_checks.json` |
| `check_nonlinear_robustness.py` | the fourth-order variant for forward data and the historical solver for the inverse, with a distinct constrained IVP | `nonlinear_robustness.json` |

Run the following in the research directory on Linux. The three shorter studies take seconds to minutes. The multi-resolution nonlinear comparison runs substantially longer on a single core.

```sh
python -B check_snapshots.py > snapshot_checks.json
python -B check_geometry.py > geometry_checks.json
python -B check_nonlinear.py > nonlinear_checks.json
python -B check_nonlinear_robustness.py > nonlinear_robustness.json
```

Generate the plot data directly from the recorded modal checks.

```sh
jq -r '["k","fixed_aligned","fixed_offset","band_offset"], (.gain_curve[]|[.k,.fixed_aligned,.fixed_offset,.band_offset]) | join(",")' snapshot_checks.json > mode_gains.csv
```

Compile `snapshot_observability.tex` with pdfLaTeX and PGFPlots. The main proof uses continuum equations. The snapshot script illustrates its predictions at finite bandwidth. The geometry script uses the historical discretization for the reduced equations and corrected Einstein–scalar observables. The nonlinear measurements use a finite-resolution wave-map score. In the main primary experiment the different restart stencil still uses observations generated with the historical solver. In the separate distinct-IVP test the forward reference uses fourth-order spatial differences and the inverse uses second-order differences. These routes share a time integrator. No claim of a certified continuum nonlinear inverse is made.

The auxiliary metric initialization error in the first distinct-IVP test is retained in `nonlinear_robustness_preconstraint.json`. The corrected run initializes the periodic primitive of the analytic momentum source. It reproduces the same field measurements and restores the fine-grid momentum check.

[Protocol](PROTOCOL.md) · [Research summary](RESEARCH_BRIEF.md)

The frozen-input path change updates the verifier checksum in `source_hashes.txt`. The original verifier is preserved in commit `60e34cd` with SHA-256 `240a490fbb2b3124d5b7eb7506de9e03744d7ed25da380f7a5b07d9d4cb89aad`.
