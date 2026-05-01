"""Aggregate Phase 1 variance decomposition outputs into figures and tables.

Reads:  tables/phase1_variance_decomp_targeted_only.csv
        tables/phase1_variance_decomp_with_controls.csv
        tables/phase1_mixedlm_validation.csv

Produces:
  figures/phase1_variance_components_boxplot.png
  figures/phase1_target_vs_batch_scatter.png
  figures/phase1_anova_vs_mixedlm.png
  tables/phase1_summary.csv
  tables/phase1_summary.md
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)


def main() -> None:
    targeted = pd.read_csv(TABLES / "phase1_variance_decomp_targeted_only.csv")
    with_ctrl = pd.read_csv(TABLES / "phase1_variance_decomp_with_controls.csv")

    rename = {
        "frac_target": "Target gene",
        "frac_guide_within_target": "Guide-within-target",
        "frac_batch": "Batch (gemgroup)",
        "frac_residual": "Residual (cell-level)",
    }

    # Boxplot
    long = targeted.melt(
        id_vars=["gene"],
        value_vars=list(rename.keys()),
        var_name="component",
        value_name="fraction",
    )
    long["component"] = long["component"].map(rename)
    fig, ax = plt.subplots(figsize=(8, 5))
    order = list(rename.values())
    sns.boxplot(
        data=long,
        x="component",
        y="fraction",
        order=order,
        ax=ax,
        showfliers=False,
        color="#5B8FF9",
    )
    ax.set_title("Replogle K562 essential perturb-seq\nHierarchical variance decomposition (per-gene)")
    ax.set_ylabel("Fraction of variance")
    ax.set_xlabel("")
    ax.set_ylim(0, 1)
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    fig.tight_layout()
    fig.savefig(FIGS / "phase1_variance_components_boxplot.png", dpi=150)
    fig.savefig(FIGS / "phase1_variance_components_boxplot.svg")
    plt.close(fig)
    print("wrote figures/phase1_variance_components_boxplot.png")

    # Scatter: target vs batch, per gene
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(
        targeted["frac_batch"],
        targeted["frac_target"],
        s=4,
        alpha=0.35,
        edgecolors="none",
    )
    lim = max(targeted["frac_batch"].max(), targeted["frac_target"].max())
    ax.plot([0, lim], [0, lim], color="black", linestyle=":", linewidth=1)
    ax.set_xlabel("Variance fraction: batch (technical)")
    ax.set_ylabel("Variance fraction: target gene (biological)")
    ax.set_title("Replogle K562: per-gene biological vs. batch variance")
    fig.tight_layout()
    fig.savefig(FIGS / "phase1_target_vs_batch_scatter.png", dpi=150)
    fig.savefig(FIGS / "phase1_target_vs_batch_scatter.svg")
    plt.close(fig)
    print("wrote figures/phase1_target_vs_batch_scatter.png")

    # MixedLM cross-check
    mixed_path = TABLES / "phase1_mixedlm_validation.csv"
    if mixed_path.exists():
        mixed = pd.read_csv(mixed_path)
        ok = mixed.dropna(subset=["reml_frac_target"]) if "reml_frac_target" in mixed.columns else pd.DataFrame()
        if len(ok) > 0:
            merged = ok.merge(targeted, on="gene", how="left", suffixes=("", "_anova"))
            fig, axes = plt.subplots(1, 3, figsize=(13, 4))
            for ax, (mcol, acol, label) in zip(
                axes,
                [("reml_frac_target", "frac_target", "Target gene"),
                 ("reml_frac_guide", "frac_guide_within_target", "Guide-within-target"),
                 ("reml_frac_batch", "frac_batch", "Batch")],
            ):
                ax.scatter(merged[acol], merged[mcol], s=20, alpha=0.7)
                m = max(merged[acol].max(), merged[mcol].max(), 0.01)
                ax.plot([0, m], [0, m], "k:", linewidth=1)
                ax.set_xlabel("ANOVA Type-I fraction")
                ax.set_ylabel("MixedLM REML fraction")
                ax.set_title(label)
            fig.suptitle("ANOVA vs. REML MixedLM (20 random HVGs)")
            fig.tight_layout()
            fig.savefig(FIGS / "phase1_anova_vs_mixedlm.png", dpi=150)
            plt.close(fig)
            print("wrote figures/phase1_anova_vs_mixedlm.png")

    # Summary table
    summary_rows = []
    for label, df in [("targeted-only", targeted), ("with controls", with_ctrl)]:
        for col, name in rename.items():
            summary_rows.append({
                "subset": label,
                "component": name,
                "median_fraction": df[col].median(),
                "mean_fraction": df[col].mean(),
                "p25": df[col].quantile(0.25),
                "p75": df[col].quantile(0.75),
                "p95": df[col].quantile(0.95),
            })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TABLES / "phase1_summary.csv", index=False)

    md_lines = ["# Phase 1: Replogle K562 essential variance decomposition\n",
                f"- Source: scverse mirror of Replogle 2022 (`replogle_2022_k562_essential.h5ad`)",
                f"- Single cells used: see pseudobulk script (>= 5 cells per (guide x batch))",
                f"- Pseudobulk samples: {len(targeted) and len(targeted)}",
                f"- Genes: {len(targeted)}",
                "",
                "## Variance fractions across genes (median across {} genes)".format(len(targeted)),
                "",
                "| Subset | Component | Median | Mean | P25 | P75 | P95 |",
                "|---|---|---:|---:|---:|---:|---:|"]
    for r in summary_rows:
        md_lines.append(
            f"| {r['subset']} | {r['component']} | {r['median_fraction']:.3f} "
            f"| {r['mean_fraction']:.3f} | {r['p25']:.3f} | {r['p75']:.3f} | {r['p95']:.3f} |"
        )
    (TABLES / "phase1_summary.md").write_text("\n".join(md_lines))
    print("wrote tables/phase1_summary.csv and phase1_summary.md")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
