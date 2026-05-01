"""Absolute effect-size analysis: how big is biology vs technical noise in raw units?

For Replogle (gene expression, log-normalized):
  * biological_effect_per_gene  = |median(log_expr | targeted) - median(log_expr | non-targeting controls)|
                                  where the perturbation directly knocks down the target gene's
                                  own expression -- the ON-TARGET effect.
  * technical_noise_per_gene    = std across batches of mean(log_expr | non-targeting only)
                                  -- how much a control well drifts day-to-day.
  * snr_per_gene                = biological_effect / technical_noise

For JUMP-CP (Cell Painting features):
  * biological_effect_per_feature = mean of |median(z | compound c) - median(z | DMSO)|
                                    averaged over compounds -- typical compound-vs-DMSO swing.
  * cross_site_noise_per_feature  = std across sources of median(z | DMSO)
                                    -- how much a "no-treatment" well changes between labs.
  * snr_per_feature               = biological_effect / cross_site_noise

  Computed both on raw (global zscore) and per-plate-zscored features.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "tables"


def replogle_effect_sizes() -> dict:
    print("=== Replogle effect sizes ===")
    adata = ad.read_h5ad(ROOT / "data" / "replogle_pseudobulk.h5ad")
    obs = adata.obs.copy()
    obs["target_gene"] = obs["target_gene"].astype(str)
    obs["batch"] = obs["batch"].astype(int)
    nperts = obs["nperts"].astype(int).values if "nperts" in obs.columns else np.ones(len(obs), int)
    var_names = adata.var.index.astype(str)

    # On-target effect: for each gene g that is also a knockdown target,
    # compute median(log_expr | target == g) - median(log_expr | non-targeting).
    X = adata.X.astype(np.float32)
    is_ctrl = nperts == 0
    print(f"  control pseudobulk samples: {is_ctrl.sum()};  targeted: {(~is_ctrl).sum()}")
    if is_ctrl.sum() == 0:
        # fall back: define controls as median across all
        ctrl_median = np.median(X, axis=0)
    else:
        ctrl_median = np.median(X[is_ctrl], axis=0)

    # Map var_names -> column index
    var_to_idx = {g: i for i, g in enumerate(var_names)}

    # Build per-(target_gene) median expression vector
    target_to_idx_in_X = {}
    targets = obs.loc[~is_ctrl, "target_gene"].unique()

    on_target_log_fc = []
    for t in tqdm(targets, desc="on-target FCs"):
        if t not in var_to_idx:
            continue
        col_for_target = var_to_idx[t]
        rows_for_target = (obs["target_gene"].values == t) & (~is_ctrl)
        if rows_for_target.sum() < 1:
            continue
        med = np.median(X[rows_for_target, col_for_target])
        on_target_log_fc.append({
            "target_gene": t,
            "n_pseudobulk": int(rows_for_target.sum()),
            "median_log_expr_target_when_perturbed": float(med),
            "median_log_expr_target_in_controls": float(ctrl_median[col_for_target]),
            "log_fold_change_on_target": float(med - ctrl_median[col_for_target]),
        })
    on_target_df = pd.DataFrame(on_target_log_fc)
    on_target_df["abs_log_fc"] = on_target_df["log_fold_change_on_target"].abs()
    on_target_df.to_csv(TABLES / "phase4_replogle_on_target_fc.csv", index=False)
    print(f"  wrote tables/phase4_replogle_on_target_fc.csv  ({len(on_target_df)} target genes)")

    # Technical noise: std across batches of mean(log_expr | non-targeting controls)
    if is_ctrl.sum() > 0:
        ctrl_obs = obs.loc[is_ctrl].reset_index(drop=True)
        ctrl_X = X[is_ctrl]
        # mean per batch per gene
        batch_codes, batch_levels = pd.factorize(ctrl_obs["batch"])
        n_batches = len(batch_levels)
        batch_means = np.zeros((n_batches, X.shape[1]))
        batch_counts = np.bincount(batch_codes, minlength=n_batches)
        np.add.at(batch_means, batch_codes, ctrl_X)
        batch_means = batch_means / np.maximum(batch_counts[:, None], 1)
        tech_std = batch_means.std(axis=0, ddof=1)
        tech_std_summary = pd.Series(tech_std).describe()
        print(f"  per-gene control batch-drift std (log units):\n    median {np.median(tech_std):.3f}, p90 {np.percentile(tech_std, 90):.3f}")

    # SNR: |on-target log-FC| / matched-gene technical noise (when we have it)
    on_target_df["batch_drift_log_std"] = np.nan
    if is_ctrl.sum() > 0:
        col_idx = on_target_df["target_gene"].map(var_to_idx)
        on_target_df["batch_drift_log_std"] = col_idx.map(
            dict(enumerate(tech_std))
        ).where(col_idx.notna()).values
        on_target_df["snr"] = on_target_df["abs_log_fc"] / on_target_df["batch_drift_log_std"].replace(0, np.nan)
    on_target_df.to_csv(TABLES / "phase4_replogle_on_target_fc.csv", index=False)

    # Headline numbers
    biol = on_target_df["abs_log_fc"].median()
    tech = float(np.median(tech_std)) if is_ctrl.sum() > 0 else float("nan")
    median_snr = on_target_df["snr"].median() if "snr" in on_target_df.columns else float("nan")

    print(f"\n  median |on-target log-FC|:  {biol:.3f} log units  ({2**biol:.2f}x fold change)")
    print(f"  median per-gene control batch-drift:  {tech:.3f} log units")
    print(f"  median signal/noise ratio:  {median_snr:.2f}")

    return {
        "biological_effect_log2_fc_median": biol,
        "biological_effect_fold_change_median": 2 ** biol,
        "technical_noise_log_std_median": tech,
        "snr_median": float(median_snr),
        "snr_p10": float(on_target_df["snr"].quantile(0.10)) if "snr" in on_target_df.columns else float("nan"),
        "snr_p90": float(on_target_df["snr"].quantile(0.90)) if "snr" in on_target_df.columns else float("nan"),
        "n_target_genes": int(len(on_target_df)),
    }


def jumpcp_effect_sizes() -> dict:
    print("\n=== JUMP-CP effect sizes ===")
    files = sorted((ROOT / "data" / "jump_cp" / "profiles").glob("*.parquet"))
    frames = [pd.read_parquet(f) for f in files]
    common = set(frames[0].columns)
    for fr in frames[1:]:
        common &= set(fr.columns)
    full = pd.concat([fr[list(common)] for fr in frames], axis=0, ignore_index=True)
    well_meta = pd.read_csv(ROOT / "data" / "jump_cp" / "well.csv")
    full = full.merge(
        well_meta[["Metadata_Source", "Metadata_Plate", "Metadata_Well", "Metadata_JCP2022"]],
        on=["Metadata_Source", "Metadata_Plate", "Metadata_Well"], how="left",
    )
    full = full.dropna(subset=["Metadata_JCP2022"]).reset_index(drop=True)
    feature_cols = [c for c in full.columns
                    if c.startswith(("Cells_", "Nuclei_", "Cytoplasm_", "Image_"))
                    and not c.startswith(("Image_FileName", "Image_PathName"))]
    Y_raw = full[feature_cols].astype(np.float64).values
    print(f"  matrix: {Y_raw.shape}")

    # JUMP-CP DMSO id
    dmso_id = "JCP2022_033924"  # negcon DMSO
    is_dmso = (full["Metadata_JCP2022"] == dmso_id).values
    print(f"  DMSO wells: {is_dmso.sum()}")

    # Global per-feature standardization
    med_g = np.nanmedian(Y_raw, axis=0)
    mad_g = np.nanmedian(np.abs(Y_raw - med_g), axis=0) * 1.4826 + 1e-9
    Y_global = np.clip(np.nan_to_num((Y_raw - med_g) / mad_g, nan=0.0, posinf=0.0, neginf=0.0), -20, 20)

    # Per-plate
    plates = pd.factorize(full["Metadata_Plate"])[0]
    Y_plate = np.empty_like(Y_raw)
    for p in tqdm(np.unique(plates), desc="per-plate z"):
        sel = plates == p
        Y_p = Y_raw[sel]
        med_p = np.nanmedian(Y_p, axis=0)
        mad_p = np.nanmedian(np.abs(Y_p - med_p), axis=0) * 1.4826 + 1e-9
        Y_plate[sel] = (Y_p - med_p) / mad_p
    Y_plate = np.clip(np.nan_to_num(Y_plate, nan=0.0, posinf=0.0, neginf=0.0), -20, 20)

    out = {}
    for label, Y in [("raw_global_zscore", Y_global), ("per_plate_zscore", Y_plate)]:
        # biological_effect: mean over compounds of |median(z) - median(z|DMSO)|
        compounds = full["Metadata_JCP2022"].astype(str).values
        cmp_codes, cmp_levels = pd.factorize(compounds)
        n_cmp = len(cmp_levels)
        # per-compound median
        idx_dmso = list(cmp_levels).index(dmso_id) if dmso_id in cmp_levels else None
        comp_meds = np.zeros((n_cmp, Y.shape[1]))
        for c in range(n_cmp):
            sel = cmp_codes == c
            if sel.sum():
                comp_meds[c] = np.median(Y[sel], axis=0)
        if idx_dmso is None:
            dmso_med = np.median(Y, axis=0)
        else:
            dmso_med = comp_meds[idx_dmso]
        # |compound_median - dmso_median|
        biol_per_cmp_per_feat = np.abs(comp_meds - dmso_med)  # (n_cmp, p)
        # average over compounds (excluding dmso itself), then take median across features
        keep = np.array([i for i in range(n_cmp) if i != idx_dmso])
        biol_per_feat = biol_per_cmp_per_feat[keep].mean(axis=0)
        biol_median = np.median(biol_per_feat)

        # cross-site noise: std across sources of median(z | DMSO)
        if is_dmso.sum() > 0:
            sources = pd.factorize(full["Metadata_Source"])[0]
            src_levels = np.unique(sources)
            dmso_per_source = []
            for s in src_levels:
                sel = (sources == s) & is_dmso
                if sel.sum():
                    dmso_per_source.append(np.median(Y[sel], axis=0))
            dmso_per_source = np.array(dmso_per_source)
            cross_site_noise_per_feat = dmso_per_source.std(axis=0, ddof=1)
            cross_noise_median = np.median(cross_site_noise_per_feat)
        else:
            cross_noise_median = float("nan")

        # within-plate noise: std across replicate DMSO wells WITHIN a plate, then median across plates
        within_plate_noise_per_feat = []
        for p in np.unique(plates):
            sel = (plates == p) & is_dmso
            if sel.sum() >= 2:
                within_plate_noise_per_feat.append(Y[sel].std(axis=0, ddof=1))
        if within_plate_noise_per_feat:
            within_plate_noise_per_feat = np.array(within_plate_noise_per_feat)
            within_noise_median = np.median(within_plate_noise_per_feat)
        else:
            within_noise_median = float("nan")

        snr_cross = biol_median / cross_noise_median if cross_noise_median > 0 else float("nan")
        snr_within = biol_median / within_noise_median if within_noise_median > 0 else float("nan")

        print(f"\n  [{label}]")
        print(f"     median compound-vs-DMSO |delta|:        {biol_median:.3f}  (z units)")
        print(f"     median within-plate DMSO replicate std: {within_noise_median:.3f}")
        print(f"     median cross-site DMSO median std:      {cross_noise_median:.3f}")
        print(f"     SNR vs. within-plate noise:  {snr_within:.2f}")
        print(f"     SNR vs. cross-site noise:    {snr_cross:.2f}")

        out[label] = {
            "median_compound_vs_dmso_z_delta": float(biol_median),
            "median_within_plate_dmso_std": float(within_noise_median),
            "median_cross_site_dmso_drift_std": float(cross_noise_median),
            "snr_vs_within_plate": float(snr_within),
            "snr_vs_cross_site": float(snr_cross),
        }

    return out


def main() -> None:
    results = {
        "replogle": replogle_effect_sizes(),
        "jumpcp": jumpcp_effect_sizes(),
    }
    with open(TABLES / "phase4_effect_sizes.json", "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nwrote tables/phase4_effect_sizes.json")


if __name__ == "__main__":
    main()
