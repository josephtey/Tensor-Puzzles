# Phase 3: cross-modality variance decomposition comparison

Numbers are *fractions of total feature/expression variance*. Rows tagged
*this analysis* are computed in this repo from raw data; *literature* rows
are reference ranges.

| Modality | Within-site tech. | Cross-site tech. | Biological-of-interest | Residual | Source |
|---|---|---|---|---|---|
| Replogle K562 essential perturb-seq | 4.6% (gemgroup) | n/a (single site) | 14.3% (target gene) | 80.2% | this repository, src/phase1_replogle/ |
| JUMP-CP Cell Painting TARGET2 (raw, no batch correction) | 3.1% (plate within source) | 77.4% (source/site) | 2.5% (compound) | 12.4% | this repository, src/phase2_jumpcp/ |
| JUMP-CP Cell Painting TARGET2 (per-plate z-scored) | 0.2% (plate within source) | 4.3% (source/site) | 24.7% (compound) | 71.3% | this repository, src/phase2_jumpcp/ |
| Bulk RNA-seq (GTEx) | <5% | 5-10% | 60-80% (tissue) | <20% | GTEx Consortium 2017 / 2020 |
| Cross-lab scRNA-seq atlas integration | 5-15% | 20-40% | 30-50% (cell type) | remainder | Luecken et al. 2022, Nat Methods (benchmarking integration) |
| Spatial transcriptomics (Visium/Xenium) | 10-25% | 30-50% | 20-40% | remainder | Spatial Touchstone 2025 / Hartman et al. (HuBMAP) 2024 |