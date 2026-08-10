"""Verify the bounded unpolarized/spike evidence admitted to Paper E1."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).with_name("E1_UNPOLARIZED_EXTENSION_RESULT.json")


def load(path: str):
    return json.loads((ROOT / path).read_text())


def main():
    production = load("bb-cbe3/CBE3.json")
    exact = load("bb-cbe3t/CBE3T.json")
    unseen = load("bb-cbe3t/CBE3T_R1.json")
    footprint = load("bb-cbe3c/CBE3C_R1_RESULTS.json")
    deep = load("bb-cbe-deep40-cp1/CP1_GATE_B_RESULT.json")
    independent = load("bb-cbe3b/CBE3B_VERIFICATION.json")

    resolutions = (1024, 2048, 4096)
    spike_d = []
    physical_widths = []
    for n in resolutions:
        row = production["campaign"][f"depth18_N{n}"]
        spike_d.append(row["surfaces"]["2.0"]["spike"])
        widths = [c["width_cells"] for c in row["spike_components"]
                  if c["kind"] == "TRUE"]
        physical_widths.append(float(np.mean(widths) * 2 * np.pi / n))

    analytic_max = unseen["analytic"]["amended_max"]
    printed_max = unseen["analytic"]["printed_max"]
    ordinary_vtd_loss = max(
        row["D_center_max"] for row in exact["T3_information_loss"].values())
    enriched_max = max(exact["T4_spike_aware"]["sup_D_scored"].values())

    checks = {
        "nonlinear_controls_pass": bool(production["controls_pass"]),
        "localized_mismatch_above_1e-2": bool(min(spike_d) > 1e-2),
        "localized_mismatch_resolution_stable_1pct": bool(
            (max(spike_d) - min(spike_d)) / np.mean(spike_d) < 0.01),
        "spike_resolved_at_all_resolutions": bool(
            min(c["width_cells"]
                for n in resolutions
                for c in production["campaign"][f"depth18_N{n}"]["spike_components"]
                if c["kind"] == "TRUE") >= 4),
        "physical_width_resolution_stable_10pct": bool(
            (max(physical_widths) - min(physical_widths))
            / np.mean(physical_widths) < 0.10),
        "amended_exact_formulas_pass": bool(analytic_max < 1e-12),
        "printed_formula_negative_control_fails": bool(printed_max > 1e-4),
        "four_function_exact_spike_loss": bool(ordinary_vtd_loss > 1e-2),
        "enriched_centered_spike_recovery": bool(enriched_max < 1e-2),
        "translated_exact_controls_pass": bool(
            unseen["outcome"] == "CBE-3T-R1-A"
            and all(c["passes"] for c in unseen["unseen_cases"])
            and all(c["finest_translation_gap"] == 0.0
                    for c in unseen["unseen_cases"])),
        "footprint_nonclosure_retained": bool(
            footprint["outcome"] == "CBE-D-architecture"),
        "deep_numerical_stop_retained": bool(
            deep["outcome"] == "CBE-D-numerical"),
        "independent_footprint_verifier_32_of_32": bool(
            len(independent["checks"]) == 32
            and not independent["failures"]
            and all(c["ok"] for c in independent["checks"])),
    }

    result = {
        "classification": (
            "E1-U1-ADVERSE-LOCALIZED: four-function VTD is insufficient "
            "near the resolved unpolarized spike candidate; enriched exact "
            "centered-spike data recover the tested exact family; generic "
            "almost-everywhere sufficiency remains unresolved"
        ),
        "all_pass": bool(all(checks.values())),
        "checks": checks,
        "production": {
            "resolutions": list(resolutions),
            "spike_D_tau2": spike_d,
            "physical_width_radians": physical_widths,
            "fine_spike_D": production["classification"]["D_spike_fine"],
            "fine_ordinary_D": production["classification"]["D_ordinary_fine"],
            "per_class": production["classification"]["per_class"],
        },
        "exact_family": {
            "analytic_amended_max": analytic_max,
            "printed_negative_control_max": printed_max,
            "four_function_center_loss_max": ordinary_vtd_loss,
            "enriched_centered_recovery_max": enriched_max,
            "unseen_outcome": unseen["outcome"],
            "translated_finest_errors": {
                c["case"]: c["shifted"][-1]["backward_sup_D"]
                for c in unseen["unseen_cases"]
            },
            "translated_finest_gaps": {
                c["case"]: c["finest_translation_gap"]
                for c in unseen["unseen_cases"]
            },
        },
        "standing_limits": {
            "footprint": footprint["outcome"],
            "deep_campaign": deep["outcome"],
            "generic_family_membership": "not established",
            "almost_everywhere_sufficiency": "not established",
        },
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
