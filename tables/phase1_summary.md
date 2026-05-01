# Phase 1: Replogle K562 essential variance decomposition

- Source: scverse mirror of Replogle 2022 (`replogle_2022_k562_essential.h5ad`)
- Single cells used: see pseudobulk script (>= 5 cells per (guide x batch))
- Pseudobulk samples: 8563
- Genes: 8563

## Variance fractions across genes (median across 8563 genes)

| Subset | Component | Median | Mean | P25 | P75 | P95 |
|---|---|---:|---:|---:|---:|---:|
| targeted-only | Target gene | 0.143 | 0.164 | 0.118 | 0.192 | 0.294 |
| targeted-only | Guide-within-target | 0.004 | 0.004 | 0.003 | 0.005 | 0.005 |
| targeted-only | Batch (gemgroup) | 0.046 | 0.064 | 0.026 | 0.084 | 0.178 |
| targeted-only | Residual (cell-level) | 0.802 | 0.768 | 0.720 | 0.846 | 0.875 |
| with controls | Target gene | 0.140 | 0.160 | 0.115 | 0.188 | 0.289 |
| with controls | Guide-within-target | 0.007 | 0.007 | 0.007 | 0.008 | 0.009 |
| with controls | Batch (gemgroup) | 0.046 | 0.064 | 0.026 | 0.083 | 0.179 |
| with controls | Residual (cell-level) | 0.803 | 0.769 | 0.721 | 0.846 | 0.874 |