"""Severity metrics: can same-class profiles recognize themselves under noise?

Computes for each modality:
  1. mAP  (mean average precision for class retrieval -- "given query X, how
     well do its same-class siblings rank above other classes?")
  2. k-NN purity at k=5,10  (fraction of nearest neighbours sharing class)
  3. Replicate correlation  (mean pairwise cosine sim, same-class vs cross-class)

Two datasets:
  Replogle: class = target_gene, samples = (guide x batch) pseudobulk profiles
  JUMP-CP:  class = compound (Metadata_JCP2022), samples = wells; computed
            twice -- raw global zscore vs per-plate zscore -- to show the
            severity gap before/after a standard "batch correction".
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


def cosine_normalize(X: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms = np.where(norms < 1e-12, 1.0, norms)
    return X / norms


def retrieval_metrics(
    X: np.ndarray,
    labels: np.ndarray,
    k_values: tuple[int, ...] = (1, 5, 10),
    chunk: int = 512,
    exclude_self_for_class: np.ndarray | None = None,
    desc: str = "retrieval",
) -> dict:
    """Compute mAP, k-NN purity, and same-vs-cross-class similarity stats.

    X: (n, d) feature matrix
    labels: (n,) integer class labels
    exclude_self_for_class: optional (n,) integer "block" id; when given, two
        rows with the same label AND same block are not considered a
        same-class pair (e.g. exclude same-plate replicates for cross-site mAP).
    """
    X = cosine_normalize(X.astype(np.float32))
    n = X.shape[0]
    labels = np.asarray(labels)

    # Class size lookup
    label_counts = pd.Series(labels).value_counts().to_dict()

    # Aggregate streaming statistics
    ap_total = 0.0
    ap_count = 0
    knn_hits = {k: 0 for k in k_values}
    knn_total = {k: 0 for k in k_values}
    same_sim_sum = 0.0
    same_sim_n = 0
    cross_sim_sum = 0.0
    cross_sim_n = 0

    for start in tqdm(range(0, n, chunk), desc=desc):
        end = min(start + chunk, n)
        Q = X[start:end]
        sims = Q @ X.T  # (b, n)
        b = sims.shape[0]
        # mask self
        for i in range(b):
            sims[i, start + i] = -np.inf
        # optionally also mask same-block pairs
        if exclude_self_for_class is not None:
            block_q = exclude_self_for_class[start:end]
            for i in range(b):
                same_block = exclude_self_for_class == block_q[i]
                same_label = labels == labels[start + i]
                drop = same_block & same_label
                drop[start + i] = False  # already dropped
                sims[i, drop] = -np.inf

        for i in range(b):
            qi = start + i
            row = sims[i]
            valid = np.isfinite(row)
            if valid.sum() < 2:
                continue
            row_v = row[valid]
            ranks = np.argsort(-row_v)
            row_labels = labels[valid][ranks]
            same = row_labels == labels[qi]
            n_same = same.sum()
            if n_same == 0:
                continue
            # Average precision
            cum_tp = np.cumsum(same)
            precision_at_k = cum_tp / np.arange(1, len(same) + 1)
            ap = (precision_at_k[same].sum()) / n_same
            ap_total += ap
            ap_count += 1
            # k-NN purity
            for k in k_values:
                if len(same) < k:
                    continue
                knn_hits[k] += int(same[:k].sum())
                knn_total[k] += k
            # similarity stats: same vs cross (top of distribution counts)
            same_sims = row_v[ranks[:200]][same[:200]] if len(same) >= 1 else np.array([])
            cross_sims = row_v[ranks[:200]][~same[:200]] if len(same) >= 1 else np.array([])
            same_sim_sum += same_sims.sum()
            same_sim_n += same_sims.size
            cross_sim_sum += cross_sims.sum()
            cross_sim_n += cross_sims.size

    map_score = ap_total / max(ap_count, 1)
    out = {
        "mAP": map_score,
        "n_queries": ap_count,
        "n_classes": int(pd.Series(labels).nunique()),
        "n_samples": n,
        "median_class_size": float(pd.Series(label_counts).median()),
        "knn_purity": {k: knn_hits[k] / max(knn_total[k], 1) for k in k_values},
        "mean_same_class_cosine": same_sim_sum / max(same_sim_n, 1),
        "mean_cross_class_cosine": cross_sim_sum / max(cross_sim_n, 1),
    }
    return out


def run_replogle() -> dict:
    print("=== Replogle K562 essential ===")
    adata = ad.read_h5ad(ROOT / "data" / "replogle_pseudobulk.h5ad")
    obs = adata.obs.copy()
    obs["target_gene"] = obs["target_gene"].astype(str)
    obs["batch"] = obs["batch"].astype(int)
    nperts = obs["nperts"].astype(int).values if "nperts" in obs.columns else np.ones(len(obs), int)

    # Restrict to targeted perturbations with >= 3 pseudobulk replicates
    counts = pd.Series(obs.loc[nperts == 1, "target_gene"]).value_counts()
    keep = counts[counts >= 3].index
    mask = (nperts == 1) & obs["target_gene"].isin(keep).values
    print(f"  retained {mask.sum()} samples for {len(keep)} target genes")

    X = adata.X[mask, :]
    labels = pd.factorize(obs.loc[mask, "target_gene"])[0]
    blocks = obs.loc[mask, "batch"].values

    # Within-site retrieval (do not exclude same-batch pairs -- the question
    # is "can a profile recognize its same-target sibling at all").
    res_overall = retrieval_metrics(X, labels, desc="replogle within-site")
    res_overall["scenario"] = "Replogle K562 perturb-seq (within site, no correction)"

    # Cross-batch retrieval: only count same-target pairs that come from
    # DIFFERENT batches, simulating "can it generalize across days".
    res_cross = retrieval_metrics(X, labels, exclude_self_for_class=blocks,
                                   desc="replogle cross-batch")
    res_cross["scenario"] = "Replogle K562 perturb-seq (cross-batch, no correction)"

    return {"within_site": res_overall, "cross_batch": res_cross}


def run_jumpcp() -> dict:
    print("\n=== JUMP-CP TARGET2 ===")
    # Reload the full feature matrix the same way decompose.py did, but here
    # we don't run ANOVA -- we run retrieval.
    from glob import glob
    files = sorted((ROOT / "data" / "jump_cp" / "profiles").glob("*.parquet"))
    print(f"  loading {len(files)} parquet files")
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
    print(f"  matrix: {Y_raw.shape}; compounds: {full['Metadata_JCP2022'].nunique()}; sources: {full['Metadata_Source'].nunique()}")

    # Filter to compounds with replicates at >= 2 sites for cross-site retrieval
    rep_ct = full.groupby("Metadata_JCP2022")["Metadata_Source"].nunique()
    multi_site = rep_ct[rep_ct >= 2].index
    keep_mask = full["Metadata_JCP2022"].isin(multi_site).values
    print(f"  compounds present at >= 2 sites: {len(multi_site)}; wells kept: {keep_mask.sum()}")

    Y_raw = Y_raw[keep_mask]
    full_kept = full.loc[keep_mask].reset_index(drop=True)
    labels = pd.factorize(full_kept["Metadata_JCP2022"])[0]
    sources = pd.factorize(full_kept["Metadata_Source"])[0]
    plates = pd.factorize(full_kept["Metadata_Plate"])[0]

    # Two normalisations
    med_g = np.nanmedian(Y_raw, axis=0)
    mad_g = np.nanmedian(np.abs(Y_raw - med_g), axis=0) * 1.4826 + 1e-9
    Y_global = np.clip(np.nan_to_num((Y_raw - med_g) / mad_g, nan=0.0, posinf=0.0, neginf=0.0), -20, 20)

    Y_plate = np.empty_like(Y_raw)
    for p in tqdm(np.unique(plates), desc="per-plate z"):
        sel = plates == p
        if sel.sum() < 2:
            Y_plate[sel] = 0
            continue
        Y_p = Y_raw[sel]
        med_p = np.nanmedian(Y_p, axis=0)
        mad_p = np.nanmedian(np.abs(Y_p - med_p), axis=0) * 1.4826 + 1e-9
        Y_plate[sel] = (Y_p - med_p) / mad_p
    Y_plate = np.clip(np.nan_to_num(Y_plate, nan=0.0, posinf=0.0, neginf=0.0), -20, 20)

    out = {}

    # Within-plate retrieval (don't exclude any pairs)
    out["raw_within_anywhere"] = retrieval_metrics(
        Y_global, labels, desc="jumpcp raw, any-pair"
    )
    out["raw_within_anywhere"]["scenario"] = "JUMP-CP raw, retrieve same compound (any well)"

    out["per_plate_within_anywhere"] = retrieval_metrics(
        Y_plate, labels, desc="jumpcp per-plate, any-pair"
    )
    out["per_plate_within_anywhere"]["scenario"] = "JUMP-CP per-plate z-scored, retrieve same compound (any well)"

    # Cross-site retrieval: exclude same-source replicates so the only positive
    # matches are at OTHER sites
    out["raw_cross_site"] = retrieval_metrics(
        Y_global, labels, exclude_self_for_class=sources,
        desc="jumpcp raw, cross-site"
    )
    out["raw_cross_site"]["scenario"] = "JUMP-CP raw, retrieve same compound at a different site"

    out["per_plate_cross_site"] = retrieval_metrics(
        Y_plate, labels, exclude_self_for_class=sources,
        desc="jumpcp per-plate, cross-site"
    )
    out["per_plate_cross_site"]["scenario"] = "JUMP-CP per-plate z-scored, retrieve same compound at a different site"

    return out


def main() -> None:
    results = {}
    results["replogle"] = run_replogle()
    results["jumpcp"] = run_jumpcp()

    # Persist
    out_json = TABLES / "phase4_retrieval_metrics.json"
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nwrote {out_json}")

    # Flatten into a CSV/markdown summary
    rows = []
    for modality, sub in results.items():
        for key, r in sub.items():
            rows.append({
                "modality": modality,
                "key": key,
                "scenario": r.get("scenario", key),
                "mAP": round(r["mAP"], 4),
                "knn_purity_k1": round(r["knn_purity"][1], 4),
                "knn_purity_k5": round(r["knn_purity"][5], 4),
                "knn_purity_k10": round(r["knn_purity"][10], 4),
                "mean_same_class_cosine": round(r["mean_same_class_cosine"], 4),
                "mean_cross_class_cosine": round(r["mean_cross_class_cosine"], 4),
                "n_queries": r["n_queries"],
                "n_classes": r["n_classes"],
            })
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "phase4_retrieval_summary.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
