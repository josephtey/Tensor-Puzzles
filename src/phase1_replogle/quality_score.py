"""Compute the perturbation-quality score per target gene.

Per the project plan:

    perturbation_quality(target) =
        var_attributable_to(target_effect)
        / (var(target_effect) + var(batch) + var(guide_within_gene))

Operationalization: for each (target_gene, batch, guide) decomposition we take the
absolute SS contributions per gene, sum them across genes (i.e. integrate over the
expression-space response), then form the ratio per target gene.

Concretely we compute, for each *target gene* t:

    SS_target(t)  = sum over genes g of  n_t * (mean_t(g) - mean(g))^2
                                          (between-target SS contributed by t)
    SS_batch_t    = sum over genes g of  Σ_b n_{t,b} (mean_{t,b}(g) - mean_t(g))^2
                                          (batch SS within target t)
    SS_guide_t    = sum over genes g of  Σ_h n_{t,h} (mean_{t,h}(g) - mean_t(g))^2
                                          (guide SS within target t)

Then quality(t) = SS_target(t) / (SS_target(t) + SS_batch_t + SS_guide_t).

Outputs: tables/phase1_perturbation_quality.csv
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
IN = ROOT / "data" / "replogle_pseudobulk.h5ad"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"


def main() -> None:
    adata = ad.read_h5ad(IN)
    obs = adata.obs.copy()
    obs["target_gene"] = obs["target_gene"].astype(str)
    obs["guide_id"] = obs["guide_id"].astype(str)
    obs["batch"] = obs["batch"].astype(int)

    # restrict to targeted perturbations only
    if "nperts" in obs.columns:
        mask = obs["nperts"].astype(int).values == 1
    else:
        mask = np.ones(len(obs), bool)
    obs = obs.loc[mask].reset_index(drop=True)
    Y = adata.X[mask, :].astype(np.float64)
    n, p = Y.shape
    print(f"Pseudobulk samples: {n} | genes: {p}")

    grand_mean = Y.mean(axis=0, keepdims=True)
    centered = Y - grand_mean

    target_codes, target_levels = pd.factorize(obs["target_gene"])
    guide_codes, _ = pd.factorize(obs["guide_id"])
    batch_codes, _ = pd.factorize(obs["batch"])

    # SS_target per target: sum over genes
    n_targets = len(target_levels)
    sums_t = np.zeros((n_targets, p))
    counts_t = np.bincount(target_codes, minlength=n_targets)
    np.add.at(sums_t, target_codes, Y)
    means_t = sums_t / np.maximum(counts_t[:, None], 1)
    ss_target_per_t = (counts_t[:, None] * (means_t - grand_mean) ** 2).sum(axis=1)

    # Within-target residual: y - mean_t
    resid_t = Y - means_t[target_codes]

    # SS_guide per target: how much of within-target residual is explained by guide
    # (guide is nested in target). Decompose by aggregating per-(target, guide).
    pair_codes, pair_levels = pd.factorize(obs["target_gene"] + "::" + obs["guide_id"])
    n_pairs = len(pair_levels)
    sums_tg = np.zeros((n_pairs, p))
    counts_tg = np.bincount(pair_codes, minlength=n_pairs)
    np.add.at(sums_tg, pair_codes, resid_t)
    means_tg = sums_tg / np.maximum(counts_tg[:, None], 1)
    ss_guide_per_pair = (counts_tg[:, None] * (means_tg ** 2)).sum(axis=1)
    pair_target = np.array([s.split("::")[0] for s in pair_levels])
    pair_target_codes = pd.Series(pair_target).map(
        {t: i for i, t in enumerate(target_levels)}
    ).values
    ss_guide_per_t = np.bincount(pair_target_codes, weights=ss_guide_per_pair, minlength=n_targets)

    # SS_batch per target: same pattern but with batch grouping inside target
    pair_tb_codes, pair_tb_levels = pd.factorize(obs["target_gene"] + "::" + obs["batch"].astype(str))
    n_tb = len(pair_tb_levels)
    sums_tb = np.zeros((n_tb, p))
    counts_tb = np.bincount(pair_tb_codes, minlength=n_tb)
    np.add.at(sums_tb, pair_tb_codes, resid_t)
    means_tb = sums_tb / np.maximum(counts_tb[:, None], 1)
    ss_batch_per_pair = (counts_tb[:, None] * (means_tb ** 2)).sum(axis=1)
    tb_target = np.array([s.split("::")[0] for s in pair_tb_levels])
    tb_target_codes = pd.Series(tb_target).map(
        {t: i for i, t in enumerate(target_levels)}
    ).values
    ss_batch_per_t = np.bincount(tb_target_codes, weights=ss_batch_per_pair, minlength=n_targets)

    quality = ss_target_per_t / np.maximum(
        ss_target_per_t + ss_guide_per_t + ss_batch_per_t, 1e-12
    )

    df = pd.DataFrame({
        "target_gene": target_levels,
        "n_pseudobulk_samples": counts_t,
        "ss_target": ss_target_per_t,
        "ss_guide_within_target": ss_guide_per_t,
        "ss_batch_within_target": ss_batch_per_t,
        "perturbation_quality": quality,
    })
    # Drop singleton targets (n=1) where within-target SS is trivially 0
    df = df[df["n_pseudobulk_samples"] >= 4].copy()
    df = df.sort_values("perturbation_quality", ascending=False).reset_index(drop=True)
    print(f"\nKept {len(df)} target genes with >= 4 pseudobulk samples")

    out_path = TABLES / "phase1_perturbation_quality.csv"
    df.to_csv(out_path, index=False)
    print(f"wrote {out_path}")

    # Plot distribution + bottom decile callout
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["perturbation_quality"], bins=50, color="#5B8FF9", edgecolor="white")
    p10 = df["perturbation_quality"].quantile(0.10)
    ax.axvline(p10, color="red", linestyle="--", label=f"10th pct = {p10:.2f}")
    ax.set_xlabel("Perturbation quality score (target SS / (target+guide+batch SS))")
    ax.set_ylabel("Number of target genes")
    ax.set_title("Replogle K562: perturbation-quality across {} target genes".format(len(df)))
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGS / "phase1_perturbation_quality_hist.png", dpi=150)
    plt.close(fig)
    print("wrote figures/phase1_perturbation_quality_hist.png")

    bottom = df.tail(20)
    top = df.head(20)
    print("\n=== Top 20 perturbations by quality ===")
    print(top.to_string(index=False))
    print("\n=== Bottom 20 perturbations by quality ===")
    print(bottom.to_string(index=False))


if __name__ == "__main__":
    main()
