# Paper E1 reproducibility manifest

The [November continuation](research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/READOUT.md) contains the current working manuscript, evidence, and self-contained research replay. This manifest also records the historical E1 workflows.

## Environment

- Python 3.14.4
- NumPy 2.5.1
- Matplotlib 3.11.1
- pypdf 6.14.2
- XeTeX 3.141592653-2.6-0.999998 (TeX Live 2026)

Pinned Python packages are listed in `requirements-lock.txt`.

## Primary sources

- `PAPER_E1_MANUSCRIPT_DRAFT_01_2026-08-09.md`: manuscript source.
- `build_review_pdf.py`: submission-candidate PDF builder.
- `CQG_COVER_LETTER_2026-08-09.tex`: cover-letter source.
- `figures/build_figures.py`: deterministic figure pipeline.
- `verification/`: exact benchmark, constraint, metric-route, modulus,
  depth/norm-sensitivity, unpolarized-extension admission, and result-chain
  records and scripts.
- `reconstruction-adjudication/`: predictive holdout, ablation, global-mode,
  and Jacobian specifications, scripts, and outputs.
- `reproduction/inputs/bb-cbe2i/`: frozen production evolution used read-only
  by the Paper E1 verification scripts.
- Under `reproduction/inputs/`, the `bb-cbe3/`, `bb-cbe3t/`, and `bb-cbe3b/`
  directories contain nonlinear unpolarized production, exact spike controls,
  and independent footprint verification.
- Under `reproduction/inputs/`,
  `bb-cbe3c/CBE3C_R1_RESULTS.json` and `bb-cbe-deep40-cp1/CP1_GATE_B_RESULT.json`: preserved architectural and
  numerical stop records; neither is treated as favorable physical evidence.

## Reproduction commands

From the repository root, using the [documented reproduction environment](reproduction/README.md):

```sh
python -B verification/run_e1_depth_norm_sensitivity.py
python -B verification/verify_e1_unpolarized_extension.py
python -B verification/verify_e1_results.py
python -B figures/build_figures.py
python -B build_review_pdf.py
```

The first command reproduces the shallow-depth and norm-sensitivity extension.
The two verifiers read the bundled frozen records. The figure and manuscript
builders consume machine-readable records. To regenerate the original
nonlinear results, run these from a scratch copy of `reproduction/inputs/`:

```sh
python -B bb-cbe3/run_cbe3.py
python -B bb-cbe3t/run_cbe3t_r1.py
python -B bb-cbe3b/verify_cbe3b.py
```

## Gate edge case

- Original record: `verification/E1_FG1_RESULT.json`.
- Original predicate implementation: `verification/run_e1_fg1.py`.
- Separate Boolean correction: `verification/adjudicate_e1_fg1.py`.
- Corrected record: `verification/E1_FG1_ADJUDICATION.json`.

The original record is retained unchanged. Appendix B of the manuscript gives
the original and corrected predicates and their effect on the aggregate gate
state.

## Archive status

The manuscript source, scripts, specifications, machine-readable results,
figure sources, dependency lock, and frozen production records are public at
`https://github.com/RonBibb/line-e-boundary-state-paper`. The linked nonlinear
repositories `bb-cbe3`, `bb-cbe3t`, `bb-cbe3b`, `bb-cbe3c`, and
`bb-cbe-deep40-cp1` are also public under the `RonBibb` GitHub account. The
Zenodo concept DOI `10.5281/zenodo.21873511` resolves to the latest published
version; the initial `v0.9` record has version DOI
`10.5281/zenodo.21873512`.
