"""Independent verifier for the admitted Paper-E1 result chain."""

from __future__ import annotations

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
ACADEMIC = PAPER.parent
OUT = HERE / "E1_RESULT_VERIFICATION.json"


def read(path):
    return json.loads(path.read_text())


def check(name, condition, detail):
    return {"name": name, "passes": bool(condition), "detail": detail}


def main():
    cbe1 = read(ACADEMIC / "bb-cbe1" / "CBE1_STAGE0_STAGE1.json")
    cbe2 = read(ACADEMIC / "bb-cbe2i" / "CBE2I.json")
    ra1 = read(PAPER / "reconstruction-adjudication" / "E1_RA1_RESULT.json")
    ra2 = read(PAPER / "reconstruction-adjudication" / "E1_RA2_RESULT.json")
    ra3 = read(PAPER / "reconstruction-adjudication" / "E1_RA3_RESULT.json")
    cj1 = read(HERE / "E1_CJ1_RESULT.json")
    fg1 = read(HERE / "E1_FG1_ADJUDICATION.json")
    fg1x = read(HERE / "E1_FG1_EXTENSION_RESULT.json")
    holdout_metric = read(HERE / "E1_HOLDOUT_METRIC_ROUTE_RESULT.json")
    clambda_est = read(HERE / "E1_CLAMBDA_ESTIMATOR_RESULT.json")
    production_clambda = read(HERE / "E1_PRODUCTION_CLAMBDA_RESULT.json")
    depth_norm = read(HERE / "E1_DEPTH_NORM_SENSITIVITY_RESULT.json")
    unpolarized = read(HERE / "E1_UNPOLARIZED_EXTENSION_RESULT.json")

    checks = []
    stage1 = cbe1["stage1"]
    checks.append(check(
        "CBE1 stage-0 controls",
        cbe1["stage0"]["flrw_pair"]["passes"]
        and cbe1["stage0"]["schwarzschild_control"]["passes"]
        and cbe1["stage0"]["passes"],
        cbe1["stage0"],
    ))
    checks.append(check(
        "CBE1 two-clock fine gate",
        stage1["fine"]["D_lam_rho"] < 1.0e-6
        and stage1["fine"]["D_lam_R"] < 1.0e-6
        and stage1["passes"],
        stage1["fine"],
    ))
    checks.append(check(
        "CBE2I controls",
        cbe2["controls_pass"],
        cbe2["controls"],
    ))
    for depth in (16, 24):
        values = []
        for n in (256, 512, 1024):
            rows = cbe2["campaign"][f"depth{depth}_N{n}"]["surfaces"].values()
            values.append(max(v["D"] for v in rows))
        ratios = (values[0] / values[1], values[1] / values[2])
        checks.append(check(
            f"CBE2I depth-{depth} second-order trend",
            all(3.5 < r < 4.5 for r in ratios) and values[-1] < 1.0e-3,
            {"D": values, "ratios": ratios},
        ))
    checks.append(check(
        "RA1 baseline and holdout",
        ra1["summary"]["baseline_reproduction"]["passes"]
        and ra1["summary"]["holdout_V1_A0"]["below_gate"]
        and ra1["summary"]["holdout_V1_A0"]["decreases_with_resolution"],
        ra1["summary"],
    ))
    checks.append(check(
        "RA1 anchor ablation fails as classified",
        not ra1["summary"]["anchor_ablation_V0_A1_below_gate"]
        and not ra1["summary"]["strict_V1_A1"]["below_gate"],
        ra1["summary"]["D_max_by_arm"],
    ))
    checks.append(check(
        "RA2 modulus scaling",
        ra2["scaling_laws_pass"]
        and not ra2["period_preserving_transform_exists_for_nonzero_c"],
        {"classification": ra2["classification"], "rows": ra2["rows"]},
    ))
    checks.append(check(
        "RA3 rank and null directions",
        ra3["expected_ranks_pass"] and ra3["direction_gates_pass"]
        and ra3["resolution_stability_pass"] and ra3["linearity_pass"],
        {k: ra3[k] for k in (
            "classification", "expected_ranks_pass", "direction_gates_pass",
            "resolution_stability_pass", "linearity_pass"
        )},
    ))
    checks.append(check(
        "CJ1 constrained reduction and a6 correction",
        cj1["a6_scaling_supported"] and not cj1["a4_explanation_supported"]
        and not cj1["comparison_statistic_used_by_reduction"],
        {k: cj1[k] for k in (
            "classification", "median_fitted_slope",
            "a6_scaling_supported", "a4_explanation_supported"
        )},
    ))
    checks.append(check(
        "FG1 exact Fuchsian and geometric benchmark",
        fg1["corrected_all_pass"]
        and fg1["governing_scientific_classification"]
        == "E1-FG1-FUCHSIAN-GEOMETRIC-A-R"
        and fg1["original_classification"] == "E1-FG1-NONCLOSURE",
        {
            "governing_classification": fg1["governing_scientific_classification"],
            "original_classification": fg1["original_classification"],
            "malformed_gate": fg1["malformed_gate"],
            "corrected_gates": fg1["corrected_gates"],
        },
    ))
    fg_orders = fg1x["convergence"]["exact_fuchsian"]
    checks.append(check(
        "FG1 third resolution and metric-route decomposition",
        all(
            1.8 < order < 2.2
            for key in ("full_D", "lambda", "momentum_residual")
            for order in fg_orders[key]["orders"]
        )
        and abs(
            fg1x["rows"]["512"]["window_asymptotic_metric"]["full_D"]
            - fg1x["rows"]["512"]["exact_fuchsian"]["full_D"]
        ) < 1.0e-10
        and fg1x["multimode_total_field_constraint_with_cross_terms"],
        {
            "exact_orders": {
                key: fg_orders[key]["orders"]
                for key in ("full_D", "lambda", "momentum_residual")
            },
            "N512_full_distances": {
                arm: fg1x["rows"]["512"][arm]["full_D"]
                for arm in (
                    "exact_fuchsian",
                    "window_asymptotic_metric",
                    "window_finite_metric",
                )
            },
        },
    ))
    d16_metric = holdout_metric["rows"]["depth16_N1024"]["D_max"]
    d24_metric = holdout_metric["rows"]["depth24_N1024"]["D_max"]
    checks.append(check(
        "Production holdout depth-24 metric-route closure",
        d16_metric < 1.0e-5
        and d24_metric < 1.0e-5
        and abs(d24_metric - d16_metric) / d16_metric < 0.01,
        {
            "depth16_N1024_D_max": d16_metric,
            "depth24_N1024_D_max": d24_metric,
            "relative_depth_gap": abs(d24_metric - d16_metric) / d16_metric,
            "interpretation": holdout_metric["interpretation"],
        },
    ))
    checks.append(check(
        "Known-answer global-modulus estimator",
        clambda_est["zero_remainder_error_decreases_with_depth"]
        and clambda_est["exact_remainder_recovers_to_roundoff"]
        and clambda_est["fine_zero_remainder_errors"][-1] < 1.0e-12,
        {
            "known_c_lambda": clambda_est["known_c_lambda"],
            "depths": clambda_est["depths"],
            "zero_remainder_errors": clambda_est["fine_zero_remainder_errors"],
            "exact_remainder_recovers_to_roundoff":
                clambda_est["exact_remainder_recovers_to_roundoff"],
        },
    ))
    checks.append(check(
        "Production global-modulus depth consistency",
        production_clambda["second_order_spatial_consistency"]
        and production_clambda["fine_depth_gap"] < 1.0e-8
        and production_clambda["continuum_depth_gap"] < 1.0e-8,
        {
            "spatial_convergence": production_clambda["spatial_convergence"],
            "richardson_estimates":
                production_clambda["richardson_continuum_estimates"],
            "fine_depth_gap": production_clambda["fine_depth_gap"],
            "continuum_depth_gap": production_clambda["continuum_depth_gap"],
            "interpretation": production_clambda["interpretation"],
        },
    ))
    deep_norms = [
        depth_norm["summary"][str(depth)][norm]["values"][-1]
        for depth in (16, 24)
        for norm in (
            "identity_l2", "shape_emphasis_l2",
            "clock_emphasis_l2", "component_linf",
        )
    ]
    deep_ratios = [
        depth_norm["summary"][str(depth)][norm]["difference_ratio"]
        for depth in (16, 24)
        for norm in (
            "identity_l2", "shape_emphasis_l2",
            "clock_emphasis_l2", "component_linf",
        )
    ]
    shallow_identity = depth_norm["summary"]["8"]["identity_l2"]["values"][-1]
    checks.append(check(
        "Production shallow-depth and norm sensitivity",
        shallow_identity > 1.0e-3
        and max(deep_norms) < 2.0e-6
        and all(3.9 < ratio < 4.1 for ratio in deep_ratios),
        {
            "depth8_fine_identity": shallow_identity,
            "deep_fine_norm_max": max(deep_norms),
            "deep_difference_ratios": deep_ratios,
            "interpretation": depth_norm["interpretation"],
        },
    ))
    checks.append(check(
        "Unpolarized localized insufficiency and exact-spike controls",
        unpolarized["all_pass"]
        and unpolarized["production"]["fine_spike_D"] > 1.0e-2
        and unpolarized["exact_family"]["four_function_center_loss_max"] > 1.0e-2
        and unpolarized["exact_family"]["enriched_centered_recovery_max"] < 1.0e-2
        and unpolarized["standing_limits"]["almost_everywhere_sufficiency"]
        == "not established",
        unpolarized,
    ))

    result = {
        "verifier": "E1 result-chain verifier",
        "checks": checks,
        "passed": sum(c["passes"] for c in checks),
        "total": len(checks),
        "all_pass": all(c["passes"] for c in checks),
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
