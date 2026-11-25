# BB-CBE-3T: Exact True-Spike Known-Answer Controls

This project validates the unpolarized spike detector and reconstruction against exact true- and false-spike Gowdy families. It also preserves the unseen-parameter confirmation that qualified the first favorable control result.

## Result

**`CBE-3T-R1-A` — the amended analytic formulas, implementation-equivariance checks, and half-cell unseen reconstructions pass.**

- Re-derived true- and false-spike formulas solve the displayed field equations to numerical precision; the literal printed alternatives fail.
- Exact-family tests show that four-function VTD data loses the spike while enriched spike parameters restore it in centered controls.
- Integer cyclic shifts give identical centered and shifted errors, as required
  by discrete shift equivariance. This is an implementation check rather than
  a positional-robustness result.
- A separate Fourier half-cell shift samples the spike off the grid nodes.
  Both unseen families converge at approximately second order and fall below
  the \(10^{-2}\) reconstruction threshold on the finest grid.

## Contents

- `run_cbe3t.py` — exact-family control battery.
- `run_cbe3t_r1.py` — unseen-parameter confirmation.
- `CBE3T.json`, `CBE3T_R1.json` — machine-readable results.
- `cbe3t_log*.txt` — preserved logs.

## Reproduction

Place this repository beside `bb-cbe3`, then install Python 3, NumPy, SciPy, and SymPy:

```text
parent/
├── bb-cbe3/
└── bb-cbe3t/
```

```bash
python3 run_cbe3t.py
python3 run_cbe3t_r1.py
```

## Scope

This validates a class of exact spike controls. It does not establish that the production feature in `bb-cbe3` belongs to that family.
