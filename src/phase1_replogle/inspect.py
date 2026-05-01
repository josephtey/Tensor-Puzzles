"""Inspect Replogle K562 essential h5ad metadata schema.

Loads the file via backed mode (no full materialization into memory) and reports
the obs/var schema, obs column unique counts, and a guess at which columns are
plausible random effects for the variance decomposition.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

DATA_PATH = Path(__file__).resolve().parents[2] / "data" / "replogle_2022_k562_essential.h5ad"
OUT_DIR = Path(__file__).resolve().parents[2] / "tables"
OUT_DIR.mkdir(exist_ok=True)


def _summarize_obs(obs: pd.DataFrame, max_unique: int = 20) -> pd.DataFrame:
    rows = []
    for col in obs.columns:
        s = obs[col]
        n_unique = s.nunique(dropna=True)
        sample_vals = s.dropna().unique()[:max_unique].tolist()
        sample_repr = ", ".join(repr(v) for v in sample_vals)
        if len(sample_repr) > 200:
            sample_repr = sample_repr[:200] + " ..."
        rows.append({
            "column": col,
            "dtype": str(s.dtype),
            "n_unique": int(n_unique),
            "n_missing": int(s.isna().sum()),
            "sample_values": sample_repr,
        })
    return pd.DataFrame(rows)


def main() -> None:
    print(f"Loading {DATA_PATH} (backed='r')")
    adata = ad.read_h5ad(DATA_PATH, backed="r")

    print(f"\nshape: {adata.shape}  (n_cells x n_genes)")
    print(f"X dtype: {adata.X.dtype if hasattr(adata.X, 'dtype') else type(adata.X)}")
    print(f"layers: {list(adata.layers.keys())}")
    print(f"obsm keys: {list(adata.obsm.keys())}")
    print(f"uns keys: {list(adata.uns.keys())}")

    print("\n=== obs schema ===")
    obs_summary = _summarize_obs(adata.obs)
    print(obs_summary.to_string(index=False))
    obs_summary.to_csv(OUT_DIR / "phase1_obs_schema.csv", index=False)

    print("\n=== var schema ===")
    var_summary = _summarize_obs(adata.var)
    print(var_summary.to_string(index=False))
    var_summary.to_csv(OUT_DIR / "phase1_var_schema.csv", index=False)

    candidates = {
        "perturbation_target": ["gene", "target", "target_gene", "perturbation", "gene_target"],
        "guide": ["guide", "guide_identity", "sgRNA", "grna", "guide_id"],
        "batch": ["gemgroup", "gem_group", "lane", "batch", "channel", "sample"],
        "cell_cycle": ["cell_cycle", "phase", "cc_phase", "cycle"],
    }
    print("\n=== candidate metadata columns ===")
    found = {}
    for role, names in candidates.items():
        hits = [c for c in adata.obs.columns if c.lower() in {n.lower() for n in names}]
        partial = [c for c in adata.obs.columns
                   if any(n.lower() in c.lower() for n in names) and c not in hits]
        found[role] = {"exact": hits, "partial": partial}
        print(f"  {role:>20s}  exact={hits}  partial={partial}")

    with open(OUT_DIR / "phase1_metadata_candidates.json", "w") as f:
        json.dump(found, f, indent=2)

    print(f"\nWrote tables to {OUT_DIR}")


if __name__ == "__main__":
    main()
