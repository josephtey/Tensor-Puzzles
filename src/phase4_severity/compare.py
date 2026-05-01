"""Phase 4 deliverable: combine retrieval and effect-size results into a single
severity-focused comparative table.

Reads:
  tables/phase4_retrieval_metrics.json
  tables/phase4_retrieval_summary.csv
  tables/phase4_effect_sizes.json

Writes:
  tables/phase4_severity.csv
  tables/phase4_severity.md
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "tables"


# Reference numbers from the literature.  Ranges are illustrative; tighten with
# primary analysis if needed.
LITERATURE = [
    {
        "modality": "Bulk RNA-seq, cross-lab (GTEx-style)",
        "mAP_or_proxy": "~0.95+",
        "biological_effect_size": "tissue effects: 5-20x fold change",
        "technical_noise_size": "<2x cross-lab",
        "snr_qualitative": "biology dominates",
        "source": "GTEx 2017/2020; Lappalainen et al.",
    },
    {
        "modality": "Cross-lab scRNA-seq atlas",
        "mAP_or_proxy": "0.5-0.8 (cell-type retrieval, post-integration)",
        "biological_effect_size": "cell-type markers: 5-50x",
        "technical_noise_size": "lab + protocol: 2-5x",
        "snr_qualitative": "biology > tech but messy",
        "source": "Luecken et al. 2022, Nat Methods",
    },
]


def main() -> None:
    retrieval_csv = TABLES / "phase4_retrieval_summary.csv"
    eff_json = TABLES / "phase4_effect_sizes.json"

    rows = []

    if retrieval_csv.exists() and eff_json.exists():
        retrieval = pd.read_csv(retrieval_csv)
        with open(eff_json) as f:
            eff = json.load(f)

        # Replogle within-site
        rep_within = retrieval[retrieval["key"] == "within_site"].iloc[0]
        rep_cross = retrieval[retrieval["key"] == "cross_batch"].iloc[0]
        rep_eff = eff["replogle"]
        rows.append({
            "modality": "Replogle K562 perturb-seq (within site)",
            "mAP_or_proxy": f"{rep_within['mAP']:.3f}",
            "knn_purity_at_5": f"{rep_within['knn_purity_k5']:.3f}",
            "same_class_cosine": f"{rep_within['mean_same_class_cosine']:.3f}",
            "cross_class_cosine": f"{rep_within['mean_cross_class_cosine']:.3f}",
            "biological_effect_size": f"|on-target log-FC| median {rep_eff['biological_effect_log2_fc_median']:.2f} ({rep_eff['biological_effect_fold_change_median']:.2f}x)",
            "technical_noise_size": f"control batch-drift std median {rep_eff['technical_noise_log_std_median']:.2f} log",
            "snr_qualitative": f"median SNR {rep_eff['snr_median']:.1f}x  (P10={rep_eff['snr_p10']:.1f}, P90={rep_eff['snr_p90']:.1f})",
            "source": "this repo, src/phase4_severity/",
        })
        rows.append({
            "modality": "Replogle K562 perturb-seq (cross-batch within site)",
            "mAP_or_proxy": f"{rep_cross['mAP']:.3f}",
            "knn_purity_at_5": f"{rep_cross['knn_purity_k5']:.3f}",
            "same_class_cosine": f"{rep_cross['mean_same_class_cosine']:.3f}",
            "cross_class_cosine": f"{rep_cross['mean_cross_class_cosine']:.3f}",
            "biological_effect_size": "(same as above)",
            "technical_noise_size": "(same as above)",
            "snr_qualitative": "harder query: replicates from different gemgroups",
            "source": "this repo, src/phase4_severity/",
        })

        for key in ["raw_within_anywhere", "per_plate_within_anywhere", "raw_cross_site", "per_plate_cross_site"]:
            r = retrieval[retrieval["key"] == key].iloc[0]
            norm = "raw" if "raw" in key else "per-plate z-scored"
            scope = "any well" if "anywhere" in key else "different site"
            row = {
                "modality": f"JUMP-CP TARGET2 ({norm}, retrieve same compound at {scope})",
                "mAP_or_proxy": f"{r['mAP']:.3f}",
                "knn_purity_at_5": f"{r['knn_purity_k5']:.3f}",
                "same_class_cosine": f"{r['mean_same_class_cosine']:.3f}",
                "cross_class_cosine": f"{r['mean_cross_class_cosine']:.3f}",
                "source": "this repo, src/phase4_severity/",
            }
            # absolute effect size from eff json
            label_key = "raw_global_zscore" if "raw" in key else "per_plate_zscore"
            ej = eff["jumpcp"].get(label_key, {})
            if ej:
                row["biological_effect_size"] = f"compound-vs-DMSO median |delta| = {ej['median_compound_vs_dmso_z_delta']:.2f} z"
                if "cross_site" in key:
                    row["technical_noise_size"] = f"DMSO cross-site std = {ej['median_cross_site_dmso_drift_std']:.2f} z"
                    row["snr_qualitative"] = f"SNR vs cross-site noise: {ej['snr_vs_cross_site']:.2f}x"
                else:
                    row["technical_noise_size"] = f"DMSO within-plate std = {ej['median_within_plate_dmso_std']:.2f} z"
                    row["snr_qualitative"] = f"SNR vs within-plate noise: {ej['snr_vs_within_plate']:.2f}x"
            rows.append(row)

    rows.extend(LITERATURE)

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "phase4_severity.csv", index=False)

    md = ["# Phase 4: severity of variance — can biology recognize itself under noise?\n",
          "*mAP = mean Average Precision for retrieving same-class neighbours;",
          "1.0 = perfect, 1/n_classes = chance.* Effect sizes are absolute units",
          "(log expression for Replogle; z-units for Cell Painting). SNR =",
          "biological-effect-size / matched-noise-size.\n",
          "| Modality | mAP | k-NN@5 purity | same-class cosine | cross-class cosine | Biological effect | Technical noise | SNR |",
          "|---|---:|---:|---:|---:|---|---|---|"]
    for r in rows:
        md.append(
            f"| {r['modality']} | {r.get('mAP_or_proxy','')} | {r.get('knn_purity_at_5','')} "
            f"| {r.get('same_class_cosine','')} | {r.get('cross_class_cosine','')} "
            f"| {r.get('biological_effect_size','')} | {r.get('technical_noise_size','')} "
            f"| {r.get('snr_qualitative','')} |"
        )
    md.append("")
    md.append(f"_Source: {rows[0]['source']} (rows tagged 'this repo'); literature otherwise._")
    (TABLES / "phase4_severity.md").write_text("\n".join(md))
    print("wrote tables/phase4_severity.csv and phase4_severity.md")
    print()
    for line in md:
        print(line)


if __name__ == "__main__":
    main()
