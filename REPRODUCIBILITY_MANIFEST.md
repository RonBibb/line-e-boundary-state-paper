# Paper E1 reproducibility manifest

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
- sibling directory `bb-cbe2i/`: frozen production evolution used read-only by
  the Paper E1 verification scripts.
- sibling directories `bb-cbe3/`, `bb-cbe3t/`, and `bb-cbe3b/`: nonlinear
  unpolarized production, exact spike controls, and independent footprint
  verification.
- `bb-cbe3c/CBE3C_R1_RESULTS.json` and
  `bb-cbe-deep40-cp1/CP1_GATE_B_RESULT.json`: preserved architectural and
  numerical stop records; neither is treated as favorable physical evidence.

## Reproduction commands

Run from the common directory containing `line-e-boundary-state-paper/` and
`bb-cbe2i/`:

```text
python3 line-e-boundary-state-paper/verification/run_e1_depth_norm_sensitivity.py
python3 bb-cbe3/run_cbe3.py
python3 bb-cbe3t/run_cbe3t_r1.py
python3 bb-cbe3b/verify_cbe3b.py
python3 line-e-boundary-state-paper/verification/verify_e1_unpolarized_extension.py
python3 line-e-boundary-state-paper/verification/verify_e1_results.py
python3 line-e-boundary-state-paper/figures/build_figures.py
python3 line-e-boundary-state-paper/build_review_pdf.py
```

The first command is the shallow-depth and norm-sensitivity extension. The
next three reproduce the nonlinear unpolarized result, unseen exact-family
qualification, and independent footprint checks. The two Paper-E1 verifiers
admit the bounded extension and then verify the complete result chain. The
figure and manuscript builders consume preserved machine-readable records;
they do not rerun the physical evolution unless an explicit verification
script is invoked.

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
`bb-cbe-deep40-cp1` are also public under the `RonBibb` GitHub account. A
version-specific archival DOI may be added after Zenodo deposition.
