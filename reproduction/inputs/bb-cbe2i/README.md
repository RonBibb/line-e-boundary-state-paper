# BB-CBE-2I: Polarized Inhomogeneous Equivalence

This project extends the collapse–bang equivalence test to polarized Gowdy T³ gravity with a massless scalar. Independently constructed contracting and expanding branches are compared through spatially varying invariant asymptotic data.

## Result

**`CBE-2I-EQ` — equivalence survives gradients in the polarized-inhomogeneous sector.**

- Manufactured-solution, homogeneous-limit, and constraint-convergence controls pass.
- The comparison residual converges at second order toward the numerical floor.
- The full spatial functions describing the velocity-dominated data reproduce the invariant collapse profiles.
- Unpolarized dynamics and true spikes are outside this result.

## Contents

- `run_cbe2i.py` — frozen Gowdy evolution and comparison campaign.
- `CBE2I.json` — controls, convergence data, and formal classification.

## Reproduction

Requires Python 3 with NumPy and SciPy:

```bash
python3 run_cbe2i.py
```

## Scope

This supports a pointwise invariant equivalence in the polarized sector only. It does not establish a continuation through a singular boundary or equivalence in generic unpolarized/spiky data.
