# Variance Decomposition Across Bio-Data Modalities

Quantifies how much of the observed variance in high-dimensional biological assays
is biological-of-interest vs. within-site technical vs. cross-site technical.

## Headline result

### Variance partition (Phase 1-3): how much variance comes from each source

| Modality | Within-site tech. | Cross-site tech. | Biological signal | Residual |
|---|---:|---:|---:|---:|
| Replogle K562 perturb-seq (single site) | 4.6% (gemgroup) | n/a | **14.3%** (target gene) | 80.2% |
| JUMP-CP Cell Painting TARGET2 — raw | 3.1% (plate) | **77.4%** (site) | 2.5% (compound) | 12.4% |
| JUMP-CP Cell Painting TARGET2 — per-plate z-scored | 0.2% (plate) | 4.3% (site) | **24.7%** (compound) | 71.3% |

### Severity (Phase 4): can biology recognize itself under noise?

Variance partitions tell you the *fraction* of variance that's technical, but not whether
biology survives in practice. Phase 4 answers that with three concrete metrics: nearest-
neighbor retrieval (mAP), absolute biological effect size, and signal-to-noise ratio.

| Modality | mAP | 1-NN purity | Biological effect | Technical noise | SNR |
|---|---:|---:|---|---|---:|
| Replogle K562 (within site) | **0.06** (84× chance) | 16% | on-target log2-FC = 0.80 (1.74×) | control batch drift = 0.12 log | **3.7×** |
| JUMP-CP raw, **same compound at different site** | 0.02 | **0.0%** | 0.15 z | 1.01 z (cross-site DMSO drift) | **0.14×** |
| JUMP-CP per-plate z-scored, same compound at different site | 0.04 | 2.5% | 0.55 z | 0.17 z | 3.35× |

For raw Cell Painting, **a same-compound replicate at another lab is never the nearest
neighbor** — the cross-site DMSO drift is 7× the size of the typical compound effect.
Per-plate normalization recovers most of this — the compound-vs-DMSO swing becomes 3×
the residual cross-site noise — but only by inflating residual variance (Phase 3 table).

The story:
- **Within a single site**, perturb-seq target-gene effects are about 3x the size of
  gemgroup batch effects. Biology is workable, mAP=6% (84× better than chance) is low
  but mostly because cosine over the full transcriptome is housekeeping-dominated.
- **Across 11 imaging sites, raw**, biology is **7× smaller than the cross-site drift**
  of a no-treatment control. You cannot recover compound identity from a different lab's
  data without a correction step.
- **A standard per-plate normalization** flips the SNR from 0.14× to 3.35× — a 24×
  improvement — but this is a property of the design (TARGET2 plates have many DMSO
  controls per plate). Cell-painting datasets without such controls don't enjoy the
  same fix.

Numbers are *median fraction of variance per feature/gene*. Full per-modality and
per-feature breakdowns in `tables/`. Reference comparators (bulk RNA-seq,
cross-lab scRNA-seq atlases, spatial transcriptomics) in
`tables/phase3_comparison.md`.

The story:
- **Within a single site**, perturb-seq target-gene effects are about 3x the
  size of gemgroup batch effects, but residual cell-level noise still dominates
  (80% of variance).
- **Across 11 imaging sites**, raw Cell Painting features have almost
  all (77%) of their variance explained by *which site* generated the data —
  literal site effects swamp biology by ~30x.
- A standard per-plate normalization (median/MAD) collapses the site effect to
  4% but at the cost of inflating residual variance: it rescues compound signal
  (2.5% → 24.7%) but does not actually recover the lost variance — it just
  hides it in the residual.

This is the core "observability" argument: between-site variance is the largest
single source of variance in image-based assays, and conventional per-plate
normalization papers over rather than addresses it.

## Layout

```
data/        Raw and processed datasets (gitignored)
notebooks/   (currently empty)
src/
  phase1_replogle/   Replogle 2022 K562 perturb-seq
  phase2_jumpcp/     JUMP-CP Cell Painting TARGET2 across 11 sites
  phase3_compare.py  Cross-modality comparative table
figures/     Plots (PNG + SVG)
tables/      Summary tables (CSV + Markdown)
```

## Reproducing

```bash
pip install -r requirements.txt

# Phase 1: Replogle K562 essential perturb-seq
python -m src.phase1_replogle.inspect          # metadata schema (no download)
python -m src.phase1_replogle.pseudobulk       # 310k cells -> 21k pseudobulk samples
python -m src.phase1_replogle.fit_variance     # hierarchical ANOVA + REML spot-check
python -m src.phase1_replogle.aggregate        # tables + figures
python -m src.phase1_replogle.quality_score    # per-target perturbation-quality score

# Phase 2: JUMP-CP TARGET2 across sites
python -m src.phase2_jumpcp.download           # 22 plates from 11 sources
python -m src.phase2_jumpcp.decompose          # variance decomp under 2 normalizations

# Phase 3: cross-modality variance-partition table
python -m src.phase3_compare

# Phase 4: severity metrics (retrieval, absolute effect sizes, SNR)
python -m src.phase4_severity.retrieval        # mAP / k-NN / cosine for both modalities
python -m src.phase4_severity.effect_sizes     # log-FC, z-units, signal/noise ratios
python -m src.phase4_severity.compare          # combined severity table
```

The Replogle h5ad (~1.5 GB) is fetched from the scverse S3 mirror; the
JUMP-CP plate parquets (~344 MB total for 22 plates) are pulled from the
cellpainting-gallery AWS Open Data bucket. Total disk: ~2 GB. Total wall-clock:
~10 min on a single 8-core machine.

## Methodology

### Phase 1 (Replogle)

- **Source**: `replogle_2022_k562_essential.h5ad` from the scverse mirror
  (originally Replogle et al. 2022, *Cell*).
- **Pseudobulk**: mean log-normalized expression aggregated by
  (`guide_id` × `batch`). Min 5 cells per group → 21,544 pseudobulk samples
  retained from 152,061 of 310,385 single cells.
- **Decomposition**: hierarchical Type-I ANOVA in the order
  `target_gene → guide_within_target → batch → residual`. Components sum to 1
  by construction; vectorized over all 8,563 genes (a few seconds total).
- **Validation**: single-random-effect REML mixed models (`statsmodels.MixedLM`)
  on 5 random high-variance genes for each random effect separately. The full
  joint 3-RE REML model with 4,400+ vc levels does not converge in
  reasonable time; this is documented as a limitation. The single-RE REML
  fits all converged.

### Phase 2 (JUMP-CP)

- **Source**: 22 TARGET2 control plates (2 per source) across 11 of 13 JUMP-CP
  sources, downloaded from the `cellpainting-gallery` AWS Open Data bucket.
  Per-well JCP2022 compound IDs joined from
  `jump-cellpainting/datasets/metadata/well.csv.gz`.
- **Features**: 4,762 CellProfiler features after restricting to columns
  shared across all plates. Parsed into compartment (Cells / Cytoplasm /
  Nuclei / Image) × feature type (AreaShape / Intensity / Texture /
  RadialDistribution / Granularity / Correlation) × channel (AGP / DNA /
  ER / Mito / RNA).
- **Two normalizations** reported side-by-side:
  - **`raw_global_zscore`**: median/MAD standardization per feature globally.
    Plate and site location effects are preserved as variance to partition.
  - **`per_plate_zscore`**: median/MAD standardization per plate. The standard
    Cell Painting pipeline before compound retrieval. Site-level variance is
    largely removed by construction.
- **Decomposition**: hierarchical ANOVA, `compound → source → plate_within_source → residual`.

### Phase 3

- Reads phase 1 and phase 2 summary CSVs, fills in literature reference ranges
  for modalities not analyzed here (bulk RNA-seq, cross-lab scRNA-seq atlas
  integration, spatial transcriptomics), and writes a single comparative table.

### Phase 4 (severity)

Three complementary measures of "how big is the noise relative to biology":

- **`retrieval.py`** — mAP (mean average precision for class retrieval), k-NN
  purity, and same-class vs cross-class mean cosine. Cosine-normalizes the
  feature matrix, then for each query computes the rank of its same-class
  siblings. For JUMP-CP, the cross-site variant masks same-source replicates so
  the only positives are at *other* sites.
- **`effect_sizes.py`** — absolute biological and technical magnitudes in
  native units. For Replogle: |median(log2_expr | targeted) - median(log2_expr |
  non-targeting controls)| for the on-target gene; cross-batch drift std of
  control wells. For JUMP-CP: median compound-vs-DMSO z-delta; DMSO replicate
  std within plate; DMSO-median std across sources.
- **`compare.py`** — joins the two into a single severity table.

## Limitations and caveats

- **Hierarchical ANOVA depends on factor order.** I use
  `target → guide → batch` for Replogle and `compound → source → plate` for
  JUMP-CP, which absorbs biological signal first. Reversing the order would
  shift fractions; the ratios are interpretation-dependent.
- **REML is a sanity check, not the headline.** For the Replogle data, full
  3-random-effect REML did not converge in tractable time. Single-RE REML on
  high-variance genes returned near-zero variance components, which is
  consistent with the ANOVA decomposition that residual dominates for the
  highest-variance genes (mostly ribosomal). Genes with strong target-gene
  ANOVA signal (e.g. Mediator subunits) would be a better test set; this is
  documented as a follow-up.
- **JUMP-CP scope is limited.** Only TARGET2 plates from each source's first
  two batches were used (22 plates total). The full JUMP-CP catalog is much
  larger and would tighten estimates.
- **The "per-plate z-score" interpretation is subtle.** It does not "fix" site
  effects — it merely projects them out. The lost variance shows up in
  residual rather than being recovered as biology.
- **Literature numbers in the comparison table are reference ranges**, not
  point estimates. They're plausible bounds drawn from the cited papers and
  should be tightened with primary analysis if used in a final deliverable.
- **Phase 4 retrieval mAP for Replogle is conservative**: cosine similarity on
  full log-normalized profiles is dominated by housekeeping genes that move
  together across all perturbations. Restricting to differentially expressed
  genes per perturbation, or pre-residualizing on controls, would push the
  mAP much higher. The number reported is the floor, not the ceiling.
- **Phase 4 cross-site mAP for JUMP-CP raw being below chance for k=1** is a
  real finding, not a bug: the cosine of two same-compound wells from
  different sites is on average *lower* than the cosine of two different-
  compound wells from the same site, because site identity dominates the
  feature-space geometry.
