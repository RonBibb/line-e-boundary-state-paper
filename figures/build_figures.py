"""Build the complete Paper-E1 publication figure set from preserved results.

Outputs each figure as PDF, SVG, and 300-dpi PNG. No numerical evolution is
performed here; the pipeline consumes archived machine-readable evidence.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/e1-matplotlib-cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import NullFormatter


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
VERIFY = PAPER / "verification"
RA = PAPER / "reconstruction-adjudication"
OUTPUT = HERE / "output"

COLORS = {
    "blue": "#3266A8",
    "orange": "#D97721",
    "green": "#2D8A63",
    "red": "#B84343",
    "purple": "#7453A6",
    "gray": "#666666",
    "light": "#E8EDF3",
}


def read(path: Path):
    return json.loads(path.read_text())


def setup():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.5,
        "lines.markersize": 5,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def save(fig, stem, caption, sources):
    paths = []
    for suffix in ("pdf", "svg", "png"):
        path = OUTPUT / f"{stem}.{suffix}"
        fig.savefig(path, metadata={"Title": caption} if suffix == "pdf" else None)
        paths.append(str(path.relative_to(PAPER)))
    plt.close(fig)
    return {"figure": stem, "caption": caption, "files": paths, "sources": sources}


def figure_1():
    fig, ax = plt.subplots(figsize=(7.2, 2.7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = [
        (0.02, 0.58, 0.17, 0.25, "Contracting\nCauchy evolution"),
        (0.23, 0.58, 0.17, 0.25, "Predictive holdout\nextraction"),
        (0.44, 0.58, 0.17, 0.25, "Boundary data\n+ global modulus"),
        (0.65, 0.58, 0.14, 0.25, "Independent\nreconstruction"),
        (0.83, 0.58, 0.15, 0.25, "Geometric\nscoring"),
    ]
    for i, (x, y, w, h, label) in enumerate(boxes):
        patch = FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.015",
            facecolor=COLORS["light"] if i != 2 else "#DDEBDD",
            edgecolor=COLORS["gray"], linewidth=0.9,
        )
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=7.6)
        if i < len(boxes) - 1:
            nx = boxes[i + 1][0]
            ax.add_patch(FancyArrowPatch(
                (x + w, y + h / 2), (nx, y + h / 2), arrowstyle="-|>",
                mutation_scale=10, linewidth=1.0, color=COLORS["blue"],
            ))
    ax.annotate("bulk history discarded", xy=(0.315, 0.56), xytext=(0.315, 0.39),
                ha="center", arrowprops=dict(arrowstyle="-[", color=COLORS["red"]))
    ax.text(0.525, 0.30,
            "smooth or spike-aware asymptotic data"
            "\nconstraint compatibility enforced", ha="center", va="center", fontsize=8)
    ax.text(0.905, 0.28,
            "time-parity map\nareal-clock comparison\ndensity-scale consistency",
            ha="center", va="center", fontsize=8)
    ax.text(0.12, 0.08, "compatible control", color=COLORS["green"], ha="center", fontsize=8)
    ax.text(0.50, 0.08, "anchor / route ablation", color=COLORS["orange"], ha="center", fontsize=8)
    ax.text(0.86, 0.08, "incompatible control", color=COLORS["red"], ha="center", fontsize=8)
    return save(fig, "figure_1_method_architecture",
                "Extraction, holdout, reconstruction, and geometric-scoring protocol.",
                ["PAPER_E1_MANUSCRIPT_DRAFT_01_2026-08-09.md"])


def figure_2():
    data = read(VERIFY / "E1_DEPTH_NORM_SENSITIVITY_RESULT.json")
    ns = np.array([256, 512, 1024])
    depths = np.array([8, 12, 16, 24])
    styles = ((8, COLORS["red"], "D"), (12, COLORS["green"], "^"),
              (16, COLORS["blue"], "o"), (24, COLORS["orange"], "s"))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1))
    ax = axes[0]
    for depth, color, marker in styles:
        vals = data["summary"][str(depth)]["identity_l2"]["values"]
        ax.loglog(ns, vals, marker=marker, color=color, label=fr"$\tau_f={depth}$")
    ax.axhline(1e-3, color=COLORS["red"], linestyle=":",
               label=r"acceptance threshold $10^{-3}$")
    ax.set_xticks(ns, [str(n) for n in ns]); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Periodic grid points, $N$")
    ax.set_ylabel(r"Geometric reconstruction distance $D_{\rm geo}$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, ncol=2, fontsize=8.2)
    ax.text(0.02, 0.95, "(a) spatial refinement", transform=ax.transAxes,
            va="top", fontweight="bold")

    ax = axes[1]
    norm_styles = (("identity_l2", "identity $L_2$", COLORS["blue"], "o"),
                   ("shape_emphasis_l2", "shape emphasized", COLORS["green"], "^"),
                   ("clock_emphasis_l2", "clock emphasized", COLORS["orange"], "s"),
                   ("component_linf", r"component $L_\infty$", COLORS["purple"], "D"))
    for key, label, color, marker in norm_styles:
        vals = [data["summary"][str(d)][key]["values"][-1] for d in depths]
        ax.semilogy(depths, vals, marker=marker, color=color, label=label)
    ax.axhline(1e-3, color=COLORS["red"], linestyle=":")
    ax.set_xticks(depths)
    ax.set_xlabel(r"Extraction depth $\tau_f$ at $N=1024$")
    ax.set_ylabel(r"Weighted reconstruction distance $D_W$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=7.8)
    ax.text(0.02, 0.95, "(b) norm sensitivity", transform=ax.transAxes,
            va="top", fontweight="bold")
    fig.tight_layout(w_pad=1.1)
    return save(fig, "figure_2_production_convergence",
                "Common-surface depth and norm-sensitivity audit: shallow truncation fails while deep cases converge under all declared weights.",
                ["verification/E1_DEPTH_NORM_SENSITIVITY_RESULT.json"])


def figure_3():
    data = read(VERIFY / "E1_FG1_EXTENSION_RESULT.json")
    ns = np.array([256, 512, 1024])
    arms = [
        ("exact_fuchsian", "Exact Fuchsian", COLORS["blue"], "o"),
        ("window_asymptotic_metric", "Window + asymptotic metric", COLORS["green"], "s"),
        ("window_finite_metric", "Window + finite-slice metric", COLORS["orange"], "^"),
    ]
    metrics = [("full_D", "Full benchmark distance"),
               ("lambda", r"$\|\delta\lambda\|_\infty$"),
               ("momentum_residual", "Momentum residual")]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), sharex=True)
    for ax, (metric, ylabel) in zip(axes, metrics):
        for arm, label, color, marker in arms:
            vals = [data["rows"][str(n)][arm][metric] for n in ns]
            ax.loglog(ns, vals, marker=marker, color=color, label=label)
        ax.set_xticks(ns, [str(n) for n in ns])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("$N$")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.22)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="upper center", ncol=3,
               fontsize=8.5, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.88), w_pad=1.0)
    return save(fig, "figure_3_metric_error_budget",
                "Exact and window/asymptotic routes coincide; finite-slice metric error is larger but second order.",
                ["verification/E1_FG1_EXTENSION_RESULT.json"])


def _shape_component_sup(case, arm):
    rows = case["arms"][arm]["surfaces"].values()
    return max(max(row["component_sup"][k] for k in ("Omega", "Sigma2", "s_P"))
               for row in rows)


def figure_4():
    data = read(RA / "E1_RA1_RESULT.json")
    ns = np.array([256, 512, 1024])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25))
    for ax, depth in zip(axes, (16, 24)):
        shape_with, shape_without, full_with, full_without = [], [], [], []
        for n in ns:
            case = data["cases"][f"depth{depth}_N{n}"]
            shape_with.append(_shape_component_sup(case, "V1-A0"))
            shape_without.append(_shape_component_sup(case, "V1-A1"))
            full_with.append(case["arms"]["V1-A0"]["D_max"])
            full_without.append(case["arms"]["V1-A1"]["D_max"])
        ax.loglog(ns, shape_with, "o-", color=COLORS["blue"], label="shape, anchor retained")
        ax.loglog(ns, shape_without, "s--", color=COLORS["green"], label="shape, anchor ablated")
        ax.loglog(ns, full_with, "^-", color=COLORS["purple"], label="production, anchor retained")
        ax.loglog(ns, full_without, "D--", color=COLORS["red"], label="production, anchor ablated")
        ax.axhline(1e-3, color=COLORS["gray"], linestyle=":")
        ax.set_xticks(ns, [str(n) for n in ns])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("$N$")
        ax.set_title(fr"$\tau_f={depth}$")
        ax.grid(True, which="both", alpha=0.22)
    axes[0].text(0.03, 0.08, "(a)", transform=axes[0].transAxes, fontweight="bold")
    axes[1].text(0.03, 0.08, "(b)", transform=axes[1].transAxes, fontweight="bold")
    axes[0].set_ylabel("Maximum residual")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="upper center", ncol=2,
               fontsize=8.5, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    return save(fig, "figure_4_global_modulus_ablation",
                "Ablating the global metric datum leaves normalized shape unchanged but destroys marked-geometric closure.",
                ["reconstruction-adjudication/E1_RA1_RESULT.json"])


def figure_5():
    exact = read(VERIFY / "E1_CLAMBDA_ESTIMATOR_RESULT.json")
    prod = read(VERIFY / "E1_PRODUCTION_CLAMBDA_RESULT.json")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.05))
    depths = np.array(exact["depths"], dtype=float)
    errs = np.array(exact["fine_zero_remainder_errors"])
    model = errs[0] * (depths * np.exp(-2 * depths)) / (depths[0] * np.exp(-2 * depths[0]))
    axes[0].semilogy(depths, errs, "o-", color=COLORS["blue"], label="analytic estimator error")
    axes[0].semilogy(depths, model, "--", color=COLORS["gray"], label=r"$\tau_f e^{-2\tau_f}$")
    axes[0].set_xlabel(r"Extraction depth $\tau_f$")
    axes[0].set_ylabel(r"$|\widehat c_\lambda^{(0)}-c_\lambda|$")
    axes[0].grid(True, which="both", alpha=0.22)
    axes[0].legend(frameon=False)
    axes[0].text(0.03, 0.08, "(a) known-answer remainder",
                 transform=axes[0].transAxes, fontweight="bold")

    ns = np.array([256, 512, 1024])
    x_scaled = 1.0e6 / ns**2
    for depth, color, marker in ((16, COLORS["blue"], "o"), (24, COLORS["orange"], "s")):
        vals = [prod["estimates"][str(depth)][str(n)] for n in ns]
        axes[1].plot(x_scaled, vals, marker=marker, color=color,
                     label=fr"$\tau_f={depth}$")
        axes[1].plot(0, prod["richardson_continuum_estimates"][str(depth)],
                     marker=marker, color=color, fillstyle="none")
    axes[1].set_xlabel(r"$10^6 N^{-2}$")
    axes[1].set_ylabel(r"Production $\widehat c_\lambda^{(0)}$")
    axes[1].ticklabel_format(axis="both", style="plain", useOffset=False)
    axes[1].set_xticks([0.0, 4.0, 8.0, 12.0, 16.0])
    axes[1].grid(True, alpha=0.22)
    axes[1].legend(frameon=False)
    axes[1].text(0.03, 0.92, "(b) production refinement",
                 transform=axes[1].transAxes, va="top", fontweight="bold")
    for ax in axes:
        ax.tick_params(labelsize=9)
    fig.tight_layout(w_pad=1.2)
    return save(fig, "figure_5_modulus_recovery",
                "Known-answer remainder decay and production modulus recovery under spatial refinement.",
                ["verification/E1_CLAMBDA_ESTIMATOR_RESULT.json",
                 "verification/E1_PRODUCTION_CLAMBDA_RESULT.json"])


def figure_6():
    data = read(VERIFY / "E1_UNPOLARIZED_EXTENSION_RESULT.json")
    ns = np.array(data["production"]["resolutions"])
    distances = np.array(data["production"]["spike_D_tau2"])
    widths = np.array(data["production"]["physical_width_radians"])
    exact = data["exact_family"]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15))
    ax = axes[0]
    ax.plot(ns, distances, "o-", color=COLORS["red"], label="localized distance")
    ax.axhline(1e-2, color=COLORS["gray"], linestyle=":", label=r"adverse criterion $10^{-2}$")
    ax.set_xticks(ns, [str(n) for n in ns])
    ax.set_ylim(0, 0.9)
    ax.set_xlabel("Periodic grid points, $N$")
    ax.set_ylabel(r"Localized reconstruction distance at $\tau=2$")
    ax.grid(True, alpha=0.22)
    ax.legend(frameon=False, loc="lower right")
    twin = ax.twinx()
    twin.plot(ns, widths, "s--", color=COLORS["blue"], label="physical width")
    twin.set_ylabel("Detected width (radians)", color=COLORS["blue"])
    twin.tick_params(axis="y", colors=COLORS["blue"])
    twin.set_ylim(0, 0.04)
    ax.set_title("(a) Generic localized failure", loc="left", fontsize=9.5)

    ax = axes[1]
    labels = ["Four-function\nVTD", "Enriched centered\nspike state"]
    values = [exact["four_function_center_loss_max"],
              exact["enriched_centered_recovery_max"]]
    ax.bar(labels, values, color=[COLORS["red"], COLORS["green"]], width=0.58)
    ax.set_yscale("log")
    ax.set_ylim(1.0e-4, 1.2)
    ax.axhline(1e-2, color=COLORS["gray"], linestyle=":")
    ax.set_ylabel("Exact-family maximum distance")
    ax.grid(True, axis="y", which="both", alpha=0.22)
    for i, value in enumerate(values):
        ax.text(i, value * 1.25, f"{value:.2e}", ha="center", fontsize=8.5)
    ax.set_title("(b) Known-answer state comparison", loc="left", fontsize=9.5)

    fig.tight_layout(w_pad=1.2)
    return save(fig, "figure_6_unpolarized_spike_extension",
                "The generic unpolarized mismatch is resolution stable and localized; enriched spike data repair the exact family and translated controls close covariantly.",
                ["verification/E1_UNPOLARIZED_EXTENSION_RESULT.json"])


def main():
    setup()
    manifest = {
        "pipeline": "Paper E1 publication figures",
        "formats": ["pdf", "svg", "png"],
        "figures": [figure_1(), figure_2(), figure_3(), figure_4(), figure_5(), figure_6()],
    }
    (OUTPUT / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
