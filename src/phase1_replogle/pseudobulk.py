"""Pseudobulk Replogle K562 essential by (guide_id x batch).

Output: an AnnData saved to data/replogle_pseudobulk.h5ad whose:
  - X            : mean log-normalized expression per (guide x batch)
  - obs          : guide_id, target_gene, batch, n_cells
  - var          : copied from input (per-gene chr/start/end/cv/fano/etc)
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[2]
IN_PATH = ROOT / "data" / "replogle_2022_k562_essential.h5ad"
OUT_PATH = ROOT / "data" / "replogle_pseudobulk.h5ad"

# minimum cells per (guide x batch) cell to retain
MIN_CELLS = 5


def main() -> None:
    print(f"Loading {IN_PATH}")
    adata = ad.read_h5ad(IN_PATH)
    print(f"  shape: {adata.shape}")

    obs = adata.obs.copy()
    obs["target_gene"] = obs["gene"].astype(str)
    obs["guide_id"] = obs["guide_id"].astype(str)
    obs["batch"] = obs["batch"].astype(int)

    obs["group_key"] = obs["guide_id"].astype(str) + "||" + obs["batch"].astype(str)
    groups = obs.groupby("group_key", observed=True, sort=False)
    n_groups = groups.ngroups
    print(f"  number of (guide x batch) groups: {n_groups}")

    sizes = groups.size()
    keep_keys = sizes[sizes >= MIN_CELLS].index
    print(f"  groups with >= {MIN_CELLS} cells: {len(keep_keys)}")

    # Build mapping group_key -> integer row index
    keep_keys = pd.Index(keep_keys)
    key_to_idx = {k: i for i, k in enumerate(keep_keys)}

    n_out = len(keep_keys)
    n_genes = adata.n_vars

    obs["row_idx"] = obs["group_key"].map(key_to_idx)
    mask = obs["row_idx"].notna()
    print(f"  retaining {int(mask.sum())} of {len(obs)} cells")

    obs_keep = obs.loc[mask].copy()
    row_idx = obs_keep["row_idx"].astype(int).to_numpy()
    cell_idx = np.where(mask.to_numpy())[0]

    X = adata.X
    is_sparse = sparse.issparse(X)
    print(f"  X is sparse: {is_sparse}")

    print("  computing per-group sums and counts ...")
    sum_mat = np.zeros((n_out, n_genes), dtype=np.float64)
    counts = np.zeros(n_out, dtype=np.int64)

    chunk = 20000
    for start in tqdm(range(0, len(cell_idx), chunk)):
        sl = slice(start, start + chunk)
        rows = cell_idx[sl]
        rgroup = row_idx[sl]
        block = X[rows, :]
        if is_sparse:
            block = block.toarray()
        np.add.at(sum_mat, rgroup, block)
        np.add.at(counts, rgroup, 1)

    means = sum_mat / counts[:, None]

    # Build obs frame for pseudobulk
    out_obs_rows = []
    for k in keep_keys:
        guide_id, batch_str = k.split("||", 1)
        out_obs_rows.append({
            "group_key": k,
            "guide_id": guide_id,
            "batch": int(batch_str),
            "n_cells": 0,  # filled below
        })
    out_obs = pd.DataFrame(out_obs_rows)
    out_obs["n_cells"] = counts
    # attach target gene from any cell in the group (deterministic)
    guide_to_target = (
        adata.obs.assign(_g=adata.obs["gene"].astype(str))
        .groupby(adata.obs["guide_id"].astype(str), observed=True, sort=False)["_g"]
        .first()
        .to_dict()
    )
    out_obs["target_gene"] = out_obs["guide_id"].map(guide_to_target)
    # attach perturbation status (1=targeted, 0=non-targeting) when it varies
    if "nperts" in adata.obs.columns:
        guide_to_nperts = (
            adata.obs.groupby(adata.obs["guide_id"].astype(str), observed=True, sort=False)["nperts"]
            .first()
            .to_dict()
        )
        out_obs["nperts"] = out_obs["guide_id"].map(guide_to_nperts).astype(int)

    out_obs.set_index("group_key", inplace=True)

    out = ad.AnnData(X=means.astype(np.float32), obs=out_obs, var=adata.var.copy())
    print(f"  pseudobulk shape: {out.shape}")
    out.write_h5ad(OUT_PATH, compression="gzip")
    print(f"  wrote {OUT_PATH}  ({OUT_PATH.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
