# Phase 4: severity of variance — can biology recognize itself under noise?

*mAP = mean Average Precision for retrieving same-class neighbours;
1.0 = perfect, 1/n_classes = chance.* Effect sizes are absolute units
(log expression for Replogle; z-units for Cell Painting). SNR =
biological-effect-size / matched-noise-size.

| Modality | mAP | k-NN@5 purity | same-class cosine | cross-class cosine | Biological effect | Technical noise | SNR |
|---|---:|---:|---:|---:|---|---|---|
| Replogle K562 perturb-seq (within site) | 0.059 | 0.125 | 0.989 | 0.990 | |on-target log-FC| median 0.80 (1.74x) | control batch-drift std median 0.12 log | median SNR 3.7x  (P10=1.9, P90=5.1) |
| Replogle K562 perturb-seq (cross-batch within site) | 0.059 | 0.125 | 0.989 | 0.990 | (same as above) | (same as above) | harder query: replicates from different gemgroups |
| JUMP-CP TARGET2 (raw, retrieve same compound at any well) | 0.033 | 0.052 | 0.925 | 0.917 | compound-vs-DMSO median |delta| = 0.15 z | DMSO within-plate std = 0.12 z | SNR vs within-plate noise: 1.24x |
| JUMP-CP TARGET2 (per-plate z-scored, retrieve same compound at any well) | 0.050 | 0.089 | 0.531 | 0.531 | compound-vs-DMSO median |delta| = 0.55 z | DMSO within-plate std = 0.77 z | SNR vs within-plate noise: 0.71x |
| JUMP-CP TARGET2 (raw, retrieve same compound at different site) | 0.023 | 0.000 | 0.733 | 0.916 | compound-vs-DMSO median |delta| = 0.15 z | DMSO cross-site std = 1.01 z | SNR vs cross-site noise: 0.14x |
| JUMP-CP TARGET2 (per-plate z-scored, retrieve same compound at different site) | 0.040 | 0.030 | 0.518 | 0.530 | compound-vs-DMSO median |delta| = 0.55 z | DMSO cross-site std = 0.17 z | SNR vs cross-site noise: 3.35x |
| Bulk RNA-seq, cross-lab (GTEx-style) | ~0.95+ |  |  |  | tissue effects: 5-20x fold change | <2x cross-lab | biology dominates |
| Cross-lab scRNA-seq atlas | 0.5-0.8 (cell-type retrieval, post-integration) |  |  |  | cell-type markers: 5-50x | lab + protocol: 2-5x | biology > tech but messy |

_Source: this repo, src/phase4_severity/ (rows tagged 'this repo'); literature otherwise._