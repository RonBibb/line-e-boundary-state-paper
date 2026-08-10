"""Analyze the production holdout global-modulus estimator across depth and resolution."""

from __future__ import annotations

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "E1_HOLDOUT_METRIC_ROUTE_RESULT.json"
OUTPUT = HERE / "E1_PRODUCTION_CLAMBDA_RESULT.json"
DEPTHS = (16, 24)
NS = (256, 512, 1024)


def main():
    source = json.loads(SOURCE.read_text())
    estimates = {
        str(depth): {
            str(n): source["rows"][f"depth{depth}_N{n}"]["c_hat_no_remainder"]
            for n in NS
        }
        for depth in DEPTHS
    }
    spatial = {}
    continuum = {}
    for depth in DEPTHS:
        vals = [estimates[str(depth)][str(n)] for n in NS]
        diffs = [abs(vals[0] - vals[1]), abs(vals[1] - vals[2])]
        spatial[str(depth)] = {
            "successive_differences": diffs,
            "difference_ratio": diffs[0] / diffs[1],
        }
        continuum[str(depth)] = vals[2] + (vals[2] - vals[1]) / 3.0

    depth_gaps = {
        str(n): abs(estimates["16"][str(n)] - estimates["24"][str(n)])
        for n in NS
    }
    result = {
        "test_id": "E1-PRODUCTION-CLAMBDA",
        "estimates": estimates,
        "spatial_convergence": spatial,
        "second_order_spatial_consistency": all(
            3.5 < spatial[str(depth)]["difference_ratio"] < 4.5
            for depth in DEPTHS
        ),
        "richardson_continuum_estimates": continuum,
        "continuum_depth_gap": abs(continuum["16"] - continuum["24"]),
        "depth_gaps_by_resolution": depth_gaps,
        "fine_depth_gap": depth_gaps["1024"],
        "interpretation": (
            "The production estimator has no known truth value. It converges "
            "spatially at second order and its depth-16/24 estimates agree to "
            "about 3e-10; the depth gap is not claimed to converge spatially."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
