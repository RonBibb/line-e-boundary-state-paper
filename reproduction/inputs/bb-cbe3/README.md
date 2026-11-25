# BB-CBE-3: Unpolarized Gowdy Spike-Sector Test

This project tests the four-function velocity-dominated reconstruction in unpolarized Gowdy T³ data containing localized spike candidates. It is the first adverse stage in the Line-E collapse–bang equivalence program.

## Result

**Aggregate: `intermediate/CBE-D` with a converged adverse localized spike result.**

- Localized spike candidates show a resolution-consistent mismatch of approximately `D ≈ 0.82` under the four-function VTD reconstruction.
- The ordinary region remains formally unclassified because the original mask is contaminated by the broader spike footprint and the backward reconstruction is ill-conditioned.
- Transient spikes were absent in the production data.
- Later exact-family controls showed that the mismatch diagnoses information missing from the four-function data set, not necessarily failure of equivalence under an enriched spike state.

## Contents

- `run_cbe3.py` — unpolarized evolution, detector, reconstruction, and scoring.
- `CBE3.json` — complete machine-readable classification.
- `cbe3_run_log*.txt` — preserved execution logs.

## Reproduction

Requires Python 3 with NumPy and SciPy:

```bash
python3 run_cbe3.py
```

## Scope

The detector classified localized candidates; this repository alone does not certify exact Lim-family membership or almost-everywhere equivalence.
