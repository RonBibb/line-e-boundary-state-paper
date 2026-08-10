"""Transparent adjudication of the malformed strict-decrease FG1 gate."""

from __future__ import annotations

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "E1_FG1_RESULT.json"
OUTPUT = HERE / "E1_FG1_ADJUDICATION.json"


def main():
    original = json.loads(SOURCE.read_text())
    fuchsian_errors = [
        max(original["cases"][f"FUCHSIAN_depth{depth}_N512"]["data_errors"].values())
        for depth in (8, 12, 16)
    ]
    exact_zero = all(v == 0.0 for v in fuchsian_errors)
    corrected_gates = dict(original["gates"])
    corrected_gates["FUCHSIAN_data_depth_decreases_or_exact"] = exact_zero or original["gates"][
        "FUCHSIAN_data_depth_decreases"
    ]
    corrected_gates.pop("FUCHSIAN_data_depth_decreases")
    result = {
        "test_id": "E1-FG1",
        "original_classification": original["classification"],
        "original_all_pass": original["all_pass"],
        "malformed_gate": "FUCHSIAN_data_depth_decreases",
        "reason": (
            "The exact Fuchsian arm is initialized from analytic asymptotic data, "
            "so its data error is identically zero at every depth and cannot "
            "strictly decrease. Exact zero is stronger than the intended gate."
        ),
        "fuchsian_data_errors_by_depth": fuchsian_errors,
        "corrected_gates": corrected_gates,
        "corrected_all_pass": all(corrected_gates.values()),
        "governing_scientific_classification": (
            "E1-FG1-FUCHSIAN-GEOMETRIC-A-R" if all(corrected_gates.values())
            else "E1-FG1-NONCLOSURE"
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
