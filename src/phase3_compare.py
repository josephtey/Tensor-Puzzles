"""Phase 3: build the cross-modality comparative variance-decomposition table.

Combines:
  * Phase 1 results from tables/phase1_summary.csv  (Replogle perturb-seq)
  * Phase 2 results from tables/phase2_variance_decomp_features.csv  (JUMP-CP)
  * Literature citations for modalities that were not analyzed in this repo
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"


# Literature citations for modalities outside scope.  Estimates pulled from the
# published references; format is (modality, within-tech, cross-tech, biological,
# residual, source).  Numbers are intended as illustrative reference points.
LITERATURE = [
    {
        "modality": "Bulk RNA-seq (GTEx)",
        "subset": "literature",
        "within_site_tech": "<5%",
        "cross_site_tech": "5-10%",
        "biological_signal": "60-80% (tissue)",
        "residual": "<20%",
        "source": "GTEx Consortium 2017 / 2020",
    },
    {
        "modality": "Cross-lab scRNA-seq atlas integration",
        "subset": "literature",
        "within_site_tech": "5-15%",
        "cross_site_tech": "20-40%",
        "biological_signal": "30-50% (cell type)",
        "residual": "remainder",
        "source": "Luecken et al. 2022, Nat Methods (benchmarking integration)",
    },
    {
        "modality": "Spatial transcriptomics (Visium/Xenium)",
        "subset": "literature",
        "within_site_tech": "10-25%",
        "cross_site_tech": "30-50%",
        "biological_signal": "20-40%",
        "residual": "remainder",
        "source": "Spatial Touchstone 2025 / Hartman et al. (HuBMAP) 2024",
    },
]


def main() -> None:
    rows = []

    # Phase 1: Replogle within-site
    p1_path = TABLES / "phase1_summary.csv"
    if p1_path.exists():
        p1 = pd.read_csv(p1_path)
        targeted = p1[p1["subset"] == "targeted-only"].set_index("component")
        biol = float(targeted.loc["Target gene", "median_fraction"])
        guide = float(targeted.loc["Guide-within-target", "median_fraction"])
        batch = float(targeted.loc["Batch (gemgroup)", "median_fraction"])
        resid = float(targeted.loc["Residual (cell-level)", "median_fraction"])
        rows.append({
            "modality": "Replogle K562 essential perturb-seq",
            "subset": "this analysis (targeted-only, hierarchical ANOVA)",
            "within_site_tech": f"{batch:.1%} (gemgroup)",
            "cross_site_tech": "n/a (single site)",
            "biological_signal": f"{biol:.1%} (target gene)",
            "residual": f"{resid:.1%}",
            "source": "this repository, src/phase1_replogle/",
            "extra": f"guide-within-target: {guide:.1%}",
        })

    # Phase 2: JUMP-CP cross-site, both normalizations
    for norm, label in [
        ("raw_global_zscore", "JUMP-CP Cell Painting TARGET2 (raw, no batch correction)"),
        ("per_plate_zscore", "JUMP-CP Cell Painting TARGET2 (per-plate z-scored)"),
    ]:
        p2_path = TABLES / f"phase2_variance_decomp_features_{norm}.csv"
        if not p2_path.exists():
            continue
        p2 = pd.read_csv(p2_path)
        biol = p2["frac_compound"].median()
        cross = p2["frac_source"].median()
        within = p2["frac_plate_within_source"].median()
        resid = p2["frac_residual"].median()
        rows.append({
            "modality": label,
            "subset": "this analysis (hierarchical ANOVA across 11 sites, 22 TARGET2 plates)",
            "within_site_tech": f"{within:.1%} (plate within source)",
            "cross_site_tech": f"{cross:.1%} (source/site)",
            "biological_signal": f"{biol:.1%} (compound)",
            "residual": f"{resid:.1%}",
            "source": "this repository, src/phase2_jumpcp/",
            "extra": "",
        })

    rows.extend(LITERATURE)
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "phase3_comparison.csv", index=False)

    # markdown
    md = ["# Phase 3: cross-modality variance decomposition comparison\n",
          "Numbers are *fractions of total feature/expression variance*. Rows tagged",
          "*this analysis* are computed in this repo from raw data; *literature* rows",
          "are reference ranges.\n"]
    md.append("| Modality | Within-site tech. | Cross-site tech. | Biological-of-interest | Residual | Source |")
    md.append("|---|---|---|---|---|---|")
    for r in rows:
        md.append(
            f"| {r['modality']} | {r['within_site_tech']} | {r['cross_site_tech']} "
            f"| {r['biological_signal']} | {r['residual']} | {r['source']} |"
        )
    (TABLES / "phase3_comparison.md").write_text("\n".join(md))
    print("wrote tables/phase3_comparison.csv and phase3_comparison.md")
    print()
    for line in md:
        print(line)


if __name__ == "__main__":
    main()
