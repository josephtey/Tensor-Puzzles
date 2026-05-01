"""Variance decomposition on Replogle pseudobulk.

Hierarchical Type-I ANOVA partition with order:
    target_gene  ->  guide_within_target  ->  batch  ->  residual

This is the same ordering the user's plan specifies for the random-effect nesting.
The decomposition is computed vectorized across all genes via group-mean
projections; total time on 21k x 8.5k is a couple of minutes.

For validation, also fit a small number of statsmodels MixedLM models with
REML on a random sample of genes and report variance components alongside.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[2]
IN_PATH = ROOT / "data" / "replogle_pseudobulk.h5ad"
OUT_DIR = ROOT / "tables"
OUT_DIR.mkdir(exist_ok=True)


def _project_groupmean(y: np.ndarray, group_idx: np.ndarray, n_groups: int) -> np.ndarray:
    """For each row of y (samples x genes), subtract the group mean.

    group_idx: (n_samples,) integer group labels in [0, n_groups)
    Returns the *residual* (y - group_mean) and the SS removed (per gene).
    """
    n_samples, n_genes = y.shape
    sums = np.zeros((n_groups, n_genes), dtype=np.float64)
    counts = np.bincount(group_idx, minlength=n_groups).astype(np.float64)
    np.add.at(sums, group_idx, y)
    means = np.where(counts[:, None] > 0, sums / np.maximum(counts[:, None], 1), 0.0)
    proj = means[group_idx]  # (n_samples, n_genes)
    ss_removed = ((proj - y.mean(axis=0, keepdims=True)) ** 2 * 0).sum(axis=0)  # placeholder
    # Recompute SS_removed correctly: SS of the projection minus SS of grand mean
    grand = y.mean(axis=0, keepdims=True)
    ss_removed = ((proj - grand) ** 2).sum(axis=0)
    resid = y - proj
    return resid, ss_removed


def hierarchical_anova(
    Y: np.ndarray,
    target_idx: np.ndarray,
    guide_idx: np.ndarray,
    batch_idx: np.ndarray,
    n_targets: int,
    n_guides: int,
    n_batches: int,
) -> pd.DataFrame:
    """Returns per-gene variance-fraction breakdown.

    Order:
      1. Fit grand mean -> residual_0  (sum of (y - mean)^2 = SS_total)
      2. Subtract target-gene mean (within-grand)  -> ss_target
      3. Subtract guide mean within whatever residual remains -> ss_guide_within_target
      4. Subtract batch mean from remaining residual -> ss_batch
      5. Whatever's left -> ss_residual
    """
    Y = Y.astype(np.float64, copy=False)
    grand = Y.mean(axis=0, keepdims=True)
    centered = Y - grand
    ss_total = (centered ** 2).sum(axis=0)

    # Step 1: target_gene.  Subtract group mean per target.
    resid1, ss_target = _project_groupmean(centered, target_idx, n_targets)

    # Step 2: guide.  Now residualize on guide (which is finer than target).
    resid2, ss_guide = _project_groupmean(resid1, guide_idx, n_guides)

    # Step 3: batch.
    resid3, ss_batch = _project_groupmean(resid2, batch_idx, n_batches)

    ss_residual = (resid3 ** 2).sum(axis=0)

    # Renormalize so sums match ss_total exactly (clean up tiny float drift)
    breakdown = np.stack([ss_target, ss_guide, ss_batch, ss_residual], axis=1)
    totals = breakdown.sum(axis=1)
    fractions = breakdown / np.maximum(totals[:, None], 1e-12)
    df = pd.DataFrame(
        fractions,
        columns=["frac_target", "frac_guide_within_target", "frac_batch", "frac_residual"],
    )
    df["ss_total"] = ss_total
    df["ss_total_check"] = totals
    return df


def _fit_mixedlm_single_re(y: np.ndarray, group_labels: np.ndarray, name: str) -> dict | None:
    """Fit y ~ 1 + (1 | group); return REML variance for the random effect and residual."""
    import statsmodels.formula.api as smf

    df = pd.DataFrame({"y": y, "g": group_labels})
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            md = smf.mixedlm("y ~ 1", data=df, groups="g")
            res = md.fit(reml=True, method="lbfgs", maxiter=200)
        if not res.converged:
            return None
        var_re = float(res.cov_re.iloc[0, 0])
        var_resid = float(res.scale)
        return {f"reml_var_{name}": var_re, f"reml_var_resid_{name}_model": var_resid}
    except Exception as exc:  # pragma: no cover
        return {f"error_{name}": str(exc)}


def validate_with_mixedlm(
    Y: np.ndarray,
    sample_genes: Iterable[int],
    var_names: pd.Index,
    target_labels: np.ndarray,
    guide_labels: np.ndarray,
    batch_labels: np.ndarray,
) -> pd.DataFrame:
    """Fit one-way REML mixed models for each random effect separately.

    With ~4400 random-effect levels in a single VC model the REML solver hangs.
    Single-RE fits converge in seconds and the variance estimates can be
    compared in ratio to the hierarchical-ANOVA components for a sanity check.
    """
    rows = []
    for gi in tqdm(list(sample_genes), desc="MixedLM validation"):
        y = Y[:, gi]
        row = {"gene_idx": gi, "gene": var_names[gi]}
        for nm, labs in [("target", target_labels), ("guide", guide_labels), ("batch", batch_labels)]:
            out = _fit_mixedlm_single_re(y, labs, nm)
            if out is None:
                row[f"converged_{nm}"] = False
            else:
                row.update(out)
                row[f"converged_{nm}"] = True
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    print(f"Loading {IN_PATH}")
    adata = ad.read_h5ad(IN_PATH)
    print(f"  pseudobulk shape: {adata.shape}")

    # Use only targeted-perturbation rows (skip non-targeting controls so the
    # target-gene factor reflects real biological signal).  But also report a
    # version including controls for comparison.
    obs = adata.obs.copy()
    obs["target_gene"] = obs["target_gene"].astype(str)
    obs["guide_id"] = obs["guide_id"].astype(str)
    obs["batch"] = obs["batch"].astype(int)

    nperts = obs["nperts"].astype(int).values if "nperts" in obs.columns else np.ones(len(obs), int)

    Y_full = adata.X.astype(np.float64)
    var_names = adata.var.index

    targeted_anova_obj = None
    targeted_data = None

    for tag, mask in [("targeted_only", nperts == 1), ("with_controls", np.ones(len(obs), bool))]:
        print(f"\n--- decomposition: {tag} ({mask.sum()} samples) ---")
        sub_obs = obs.loc[mask].reset_index(drop=True)
        Y = Y_full[mask, :]

        target_codes, target_levels = pd.factorize(sub_obs["target_gene"])
        guide_codes, guide_levels = pd.factorize(sub_obs["guide_id"])
        batch_codes, batch_levels = pd.factorize(sub_obs["batch"])

        df = hierarchical_anova(
            Y,
            target_idx=target_codes,
            guide_idx=guide_codes,
            batch_idx=batch_codes,
            n_targets=len(target_levels),
            n_guides=len(guide_levels),
            n_batches=len(batch_levels),
        )
        df.insert(0, "gene", var_names)
        out_path = OUT_DIR / f"phase1_variance_decomp_{tag}.csv"
        df.to_csv(out_path, index=False)
        print(f"  wrote {out_path}")
        print("  median fractions:")
        print(df[["frac_target", "frac_guide_within_target", "frac_batch", "frac_residual"]].median().to_string())

        if tag == "targeted_only":
            targeted_anova_obj = df
            targeted_data = (Y, sub_obs)

    # Validate hierarchical ANOVA against single-RE REML on 5 random HVGs
    Y, sub_obs = targeted_data
    rng = np.random.default_rng(42)
    var_per_gene = np.var(Y, axis=0)
    top_idx = np.argsort(-var_per_gene)[:200]
    sample = rng.choice(top_idx, size=5, replace=False)
    print("\n--- Single-RE REML validation (5 HVGs) ---")
    mixed_df = validate_with_mixedlm(
        Y, sample, var_names,
        target_labels=sub_obs["target_gene"].values,
        guide_labels=sub_obs["guide_id"].values,
        batch_labels=sub_obs["batch"].astype(str).values,
    )
    mixed_path = OUT_DIR / "phase1_mixedlm_validation.csv"
    mixed_df.to_csv(mixed_path, index=False)
    print(f"  wrote {mixed_path}")
    print(mixed_df.to_string(index=False))


if __name__ == "__main__":
    main()
