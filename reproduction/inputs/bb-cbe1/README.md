# BB-CBE-1: Homogeneous Collapse–Bang Equivalence Controls

This project tests whether contracting and expanding homogeneous Einstein–massless-scalar solutions approach the same invariant asymptotic state under time reversal. It contains the Stage 0 analytic controls, Stage 1 numerical control, and the homogeneous portion of Stage 2 for Research Line E (singularity ↔ Big-Bang equivalence).

## Result

**`CBE-B_homogeneous` — exact-control equivalence with a conserved-shear-stratified basin.**

- Flat stiff-fluid collapse and independently constructed expansion agree below the frozen tolerance on two invariant clocks.
- Homogeneous perturbations do not return to one universal isotropic state. Each conserved Hubble-normalized shear value defines its own equivalence stratum.
- Radiation controls become shear dominated, confirming that quiescent behavior is special to stiff/scalar matter.
- The work does not define a continuation law, child universe, topology change, or inhomogeneous equivalence.

## Contents

- `run_cbe1_stage0_stage1.py` — analytic and numerical controls.
- `CBE1_STAGE0_STAGE1.json` — machine-readable results and defect ledger.

## Reproduction

Requires Python 3 with NumPy and SciPy:

```bash
python3 run_cbe1_stage0_stage1.py
```

The saved result was checked for JSON validity and the source was byte-compiled before archival.

## Scope

This is a symmetry-reduced mathematical control. It does not test inhomogeneous Gowdy dynamics, black-hole interiors, reheating, or any physical parent-to-child map.
