# BB-CBE-3b: Almost-Everywhere and Centered-Spike Architecture

This project evaluates the first footprint-aware, centered-spike successor to the unpolarized Gowdy test. It preserves three production resolutions, fitted diagnostics, figures, and an independent verifier.

## Result

**`3b-D` — architectural nonclosure; no physical equivalence verdict.**

- The original 99.9%-mass footprint saturated almost the entire circle and therefore failed as a feature localizer.
- A P-only window velocity retained a convergent drift, while the preregistered hyperbolic-velocity variant was substantially more stable.
- Forward prediction remained intermediate rather than conclusively passing or failing.
- The Lim fit was not a valid membership test because its training mask spanned the whole inhomogeneous field.
- Independent verification passes all 32 checks.

## Contents

- `run_cbe3b.py` — production test.
- `verify_cbe3b.py` — independent recomputation of gates and outcome.
- `CBE3B_RESULTS.json`, `CBE3B_VERIFICATION.json` — result records.
- `cbe3b_slices_N*.npz` — preserved resolution data.
- `fig_*.png` — diagnostic figures.

## Reproduction

Clone `bb-cbe3` and `bb-cbe3t` beside this repository. Requires Python 3, NumPy, and SciPy.

```bash
python3 run_cbe3b.py
python3 verify_cbe3b.py
```

## Scope

The outcome is a test-architecture failure, not evidence for or against physical collapse–bang equivalence.
