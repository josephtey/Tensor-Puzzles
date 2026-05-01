"""JUMP-CP TARGET2 cross-site variance decomposition.

Loads the downloaded TARGET2 plate parquets, joins per-well JCP2022 compound IDs,
standardizes features per-plate (so feature scaling differences across sites don't
dominate), then runs a hierarchical Type-I ANOVA in the order:

    compound  ->  source (lab/site)  ->  plate within source  ->  residual

Aggregates the per-feature fractions by (compartment x channel x feature_type).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[2]
PROFILES = ROOT / "data" / "jump_cp" / "profiles"
WELL_META = ROOT / "data" / "jump_cp" / "well.csv"
PLATE_META = ROOT / "data" / "jump_cp" / "plate.csv"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
TABLES.mkdir(exist_ok=True)
FIGS.mkdir(exist_ok=True)


def parse_feature_name(col: str) -> dict:
    """Decompose a CellProfiler feature name into compartment / type / channel."""
    parts = col.split("_")
    compartment = parts[0]
    feature_type = parts[1] if len(parts) > 1 else None
    channel = None
    for c in ("AGP", "DNA", "ER", "Mito", "RNA", "BFAGP"):
        if c in parts[2:]:
            channel = c
            break
    return {"compartment": compartment, "feature_type": feature_type, "channel": channel}


def hierarchical_anova_2level(
    Y: np.ndarray,
    factors: list[tuple[str, np.ndarray, int]],
) -> pd.DataFrame:
    """Sequential Type-I decomposition. factors: list of (name, codes, n_levels).

    SS allocated to each factor in order; residual is the leftover.
    """
    Y = Y.astype(np.float64, copy=False)
    grand = Y.mean(axis=0, keepdims=True)
    centered = Y - grand
    ss_total = (centered ** 2).sum(axis=0)
    resid = centered.copy()

    breakdown = {}
    for name, codes, k in factors:
        sums = np.zeros((k, Y.shape[1]))
        counts = np.bincount(codes, minlength=k).astype(np.float64)
        np.add.at(sums, codes, resid)
        means = np.where(counts[:, None] > 0, sums / np.maximum(counts[:, None], 1), 0.0)
        proj = means[codes]
        ss_factor = (proj ** 2).sum(axis=0)
        breakdown[f"frac_{name}"] = ss_factor
        resid = resid - proj

    breakdown["frac_residual"] = (resid ** 2).sum(axis=0)
    df = pd.DataFrame(breakdown)
    totals = df.sum(axis=1)
    for c in df.columns:
        df[c] = df[c] / np.maximum(totals, 1e-12)
    df["ss_total"] = ss_total
    return df


def main() -> None:
    well_meta = pd.read_csv(WELL_META)
    print(f"well metadata: {well_meta.shape}")

    files = sorted(PROFILES.glob("*.parquet"))
    print(f"loading {len(files)} parquet files ...")

    frames = []
    for f in tqdm(files):
        df = pd.read_parquet(f)
        # ensure source/plate cols exist
        if "Metadata_Source" not in df.columns:
            df["Metadata_Source"] = f.stem.split("__")[0]
        frames.append(df)

    # Use the smallest set of feature columns shared across plates
    common_cols = set(frames[0].columns)
    for fr in frames[1:]:
        common_cols &= set(fr.columns)
    print(f"common columns: {len(common_cols)}")

    full = pd.concat([fr[list(common_cols)] for fr in frames], axis=0, ignore_index=True)
    print(f"concatenated: {full.shape}")

    # Merge in JCP2022 compound id
    full = full.merge(
        well_meta[["Metadata_Source", "Metadata_Plate", "Metadata_Well", "Metadata_JCP2022"]],
        on=["Metadata_Source", "Metadata_Plate", "Metadata_Well"],
        how="left",
    )
    print(f"merged: {full.shape}")
    print(f"compounds: {full['Metadata_JCP2022'].nunique()} unique, {full['Metadata_JCP2022'].isna().sum()} missing")

    # Drop wells without compound assignment
    full = full.dropna(subset=["Metadata_JCP2022"]).reset_index(drop=True)

    feature_cols = [c for c in full.columns if not c.startswith(("Metadata_", "Image_FileName", "Image_PathName"))]
    feature_cols = [c for c in feature_cols if c.startswith(("Cells_", "Nuclei_", "Cytoplasm_", "Image_"))]
    print(f"feature columns: {len(feature_cols)}")

    # We compute the variance decomposition under two normalization strategies
    # to give the *honest* range across Cell Painting analysis pipelines:
    #
    #   (a) GLOBAL z-score per feature: keeps plate/site differences as
    #       variance to be partitioned.  This matches what the *raw* feature
    #       matrix looks like and reveals how dominant cross-site effects are
    #       before correction.
    #
    #   (b) PER-PLATE z-score (median/MAD): standard pipeline for compound
    #       retrieval.  Removes plate-level location/scale so what's left to
    #       partition is the *post-correction* signal.
    #
    # We report both to make the story comparable to literature.
    print("per-feature global standardization ...")
    Y = full[feature_cols].astype(np.float64).values
    med = np.nanmedian(Y, axis=0)
    mad = np.nanmedian(np.abs(Y - med), axis=0) * 1.4826 + 1e-9
    Y_global = np.clip(np.nan_to_num((Y - med) / mad, nan=0.0, posinf=0.0, neginf=0.0), -20, 20)

    print("per-plate (within-source) standardization ...")
    plate_codes_for_z, plate_levels_for_z = pd.factorize(full["Metadata_Plate"])
    Y_plate = np.empty_like(Y)
    for p in tqdm(range(len(plate_levels_for_z)), desc="plate-z"):
        mask = plate_codes_for_z == p
        if mask.sum() < 2:
            Y_plate[mask] = 0
            continue
        Y_p = Y[mask]
        med_p = np.nanmedian(Y_p, axis=0)
        mad_p = np.nanmedian(np.abs(Y_p - med_p), axis=0) * 1.4826 + 1e-9
        Y_plate[mask] = (Y_p - med_p) / mad_p
    Y_plate = np.clip(np.nan_to_num(Y_plate, nan=0.0, posinf=0.0, neginf=0.0), -20, 20)

    compound_codes, compound_levels = pd.factorize(full["Metadata_JCP2022"])
    source_codes, source_levels = pd.factorize(full["Metadata_Source"])
    full["plate_within_source"] = full["Metadata_Source"].astype(str) + "::" + full["Metadata_Plate"].astype(str)
    pls_codes, pls_levels = pd.factorize(full["plate_within_source"])
    print(f"compounds={len(compound_levels)}  sources={len(source_levels)}  plates={len(pls_levels)}")

    rename = {
        "frac_compound": "Compound (biological)",
        "frac_source": "Source/site (cross-site technical)",
        "frac_plate_within_source": "Plate (within-site technical)",
        "frac_residual": "Residual (well-level)",
    }

    summaries = {}
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for ax, (norm_name, Y_z, color) in zip(axes, [
        ("raw_global_zscore", Y_global, "#FF8C53"),
        ("per_plate_zscore", Y_plate, "#5B8FF9"),
    ]):
        print(f"\n=== Normalization: {norm_name} ===")
        decomp = hierarchical_anova_2level(
            Y_z,
            [
                ("compound", compound_codes, len(compound_levels)),
                ("source", source_codes, len(source_levels)),
                ("plate_within_source", pls_codes, len(pls_levels)),
            ],
        )
        decomp.insert(0, "feature", feature_cols)
        parsed = decomp["feature"].apply(parse_feature_name).apply(pd.Series)
        decomp = pd.concat([decomp, parsed], axis=1)
        out_path = TABLES / f"phase2_variance_decomp_features_{norm_name}.csv"
        decomp.to_csv(out_path, index=False)
        summaries[norm_name] = decomp
        print(f"wrote {out_path}")
        print("Median fractions across all features:")
        print(decomp[["frac_compound", "frac_source", "frac_plate_within_source", "frac_residual"]].median().to_string())

        cat_summary = (
            decomp.dropna(subset=["channel"])
            .groupby(["compartment", "feature_type", "channel"], observed=True)
            [["frac_compound", "frac_source", "frac_plate_within_source", "frac_residual"]]
            .median()
            .reset_index()
        )
        cat_summary.to_csv(TABLES / f"phase2_variance_by_category_{norm_name}.csv", index=False)

        long = decomp.melt(
            id_vars=["feature"],
            value_vars=["frac_compound", "frac_source", "frac_plate_within_source", "frac_residual"],
            var_name="component",
            value_name="fraction",
        )
        long["component"] = long["component"].map(rename)
        sns.boxplot(data=long, x="component", y="fraction", order=list(rename.values()), ax=ax,
                    showfliers=False, color=color)
        ax.set_title(f"JUMP-CP TARGET2 ({norm_name})")
        ax.set_ylabel("Fraction of variance")
        ax.set_xlabel("")
        ax.set_ylim(0, 1)
        plt.setp(ax.get_xticklabels(), rotation=25, ha="right")

    fig.suptitle("JUMP-CP TARGET2 plates: variance decomposition before vs. after per-plate normalization")
    fig.tight_layout()
    fig.savefig(FIGS / "phase2_variance_components_boxplot.png", dpi=150)
    fig.savefig(FIGS / "phase2_variance_components_boxplot.svg")
    plt.close(fig)
    print("wrote figures/phase2_variance_components_boxplot.png")

    # Combined summary table
    rows = []
    for norm_name, df in summaries.items():
        for col in ["frac_compound", "frac_source", "frac_plate_within_source", "frac_residual"]:
            rows.append({
                "normalization": norm_name,
                "component": rename[col],
                "median": df[col].median(),
                "mean": df[col].mean(),
                "p25": df[col].quantile(0.25),
                "p75": df[col].quantile(0.75),
            })
    pd.DataFrame(rows).to_csv(TABLES / "phase2_summary.csv", index=False)
    print("wrote tables/phase2_summary.csv")


if __name__ == "__main__":
    main()
