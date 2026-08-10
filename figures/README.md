# Paper E1 Figure Pipeline

Run from the Academic workspace root:

```bash
MPLCONFIGDIR=/tmp/e1-matplotlib-cache \
python3 line-e-boundary-state-paper/figures/build_figures.py
```

The pipeline reads preserved JSON evidence only; it does not rerun the physical
models. It writes PDF, SVG, and 300-dpi PNG versions to `figures/output/` and
records source provenance in `figure_manifest.json`.

Figures 1--4 are intended for the main manuscript. Figure 5 is optional or
supplementary if the journal page budget is restrictive.
