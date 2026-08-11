# CQG / IOP Overleaf source

The Overleaf main document is `CQG_OVERLEAF_SUBMISSION_2026-08-10.tex`.

The project uses IOP Publishing's official `iopjournal.cls` from the July 2025
template package. The generic display placeholders (`Journal Name`, `Journal
vv (yyyy) aaaaaa`, and `Author et al`) are specialized to *Classical and
Quantum Gravity* and `Ronald Bibb`; no layout or functional class code is changed.
Compile with pdfLaTeX. The repository-relative figure paths are already
configured.

Required files are:

- `CQG_OVERLEAF_SUBMISSION_2026-08-10.tex`
- `iopjournal.cls`
- `orcid.pdf`
- the six PDF figures under `figures/output/`

To use the prepared archive, upload
`output/archive/CQG_OVERLEAF_PACKAGE_2026-08-10.zip` to Overleaf as a new
project and select `CQG_OVERLEAF_SUBMISSION_2026-08-10.tex` as the main
document if Overleaf does not select it automatically.

Regenerate the journal source from the manuscript Markdown with:

```text
python3 build_cqg_overleaf_tex.py
```
