# Adipogenic Drift of the Osteoblast Lineage Peaks at the Osteopenic Stage: A Single-Cell Study of Human Bone-Mass Decline

**Single-cell evidence that adipogenic drift of the osteoblast lineage forms a continuum of hybrid intermediate states and is most active at the osteopenic stage — benchmarked against a healthy adult-mouse reference.**

This repository contains the reproducible analysis pipeline (R), the per-panel figure scripts (R), and the multi-panel figure-assembly scripts (Python) for the study:

> *Adipogenic Drift of the Osteoblast Lineage Peaks at the Osteopenic Stage: A Single-Cell Study of Human Bone-Mass Decline*

## Summary of the study

Using single-cell RNA-seq of human femoral heads (GSE169396) together with a healthy adult mouse metaphysis reference (GSE317069), we show that:

1. **45.8%** of human bone-lineage cells (3,850 cells) occupy a hybrid (osteogenic⁺ / adipogenic⁺) intermediate state — the largest single category — indicating that adipogenic drift proceeds through a **continuum** of intermediate states rather than a discrete binary switch;
2. Along the osteoblast differentiation trajectory the cell-level **DriftIndex** (= UCell adipogenic-program score − UCell osteogenic-program score) is **biphasic**: it rises to an intermediate maximum and then declines monotonically (segmented-regression breakpoint at normalized pseudotime **0.448**, 95% CI 0.390–0.506; overall Spearman ρ = −0.321, p = 3.3 × 10⁻⁹³);
3. Drift is **most active at the osteopenic stage** (DriftIndex 0.022; 56.0% hybrid) — exceeding the osteoporotic-range donor (DriftIndex 0.013; 37.1%) — and this ordering is robust across classification thresholds;
4. Healthy mouse bone-lineage cells show only **9.6%** hybrid cells and a negative DriftIndex (−0.059), providing a cross-species healthy baseline.

## Datasets (all public, from GEO)

| Dataset | Content | Use |
|---|---|---|
| **GSE169396** | Human femoral-head scRNA-seq, 4 donors (10x v3). Lumbar spine T-scores from the original study: S1 = 61 y F, T = −3.0 (osteoporotic range); S2 = 45 y F, T = −1.3 (osteopenic range); S3 = 66 y M, T = NA; S4 = 31 y M, T = +0.6. **S1 and S2 are the only donors densitometrically stratifiable under current criteria; S3 and S4 are retained as descriptive references only (see Sample classification).** | Main human analysis (steps 1–3) |
| **GSE317069** | Healthy adult C57BL/6 mouse metaphysis scRNA-seq, 5 replicates (GSM6552950–54, `Metaphysis_1–5`) | Cross-species healthy reference (step 6) |

Download (example for the mouse metaphysis samples):

```bash
# Human data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169396/suppl/GSE169396_RAW.tar
# Mouse metaphysis data (5 samples × 3 files)
for i in 50 51 52 53 54; do
  base=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6552nnn/GSM65529${i}/suppl
  rep=$((i - 49))
  wget $base/GSM65529${i}_Metaphysis_${rep}_barcodes.tsv.gz
  wget $base/GSM65529${i}_Metaphysis_${rep}_genes.tsv.gz
  wget $base/GSM65529${i}_Metaphysis_${rep}_matrix.mtx.gz
done
```

## Directory layout required to run

This repository ships the **code, the pre-computed result tables, and the rendered figures**. The raw sequencing data are *not* bundled (they are public GEO downloads — see above); create `data/raw/` yourself and the scripts will populate `data/processed/`.

```
<repo-root>/
├── data/                                       # NOT bundled — create it, then download the GEO data here
│   ├── raw/
│   │   ├── GSM5201883_S1_barcodes.tsv.gz       # GSE169396 RAW.tar extracted,
│   │   ├── ...                                  # GSM<id>_<S1..S4>_*.gz files
│   │   └── GSE317069_metaphysis/
│   │       ├── GSM6552950_Metaphysis_1_*.gz     # mouse 10x files
│   │       └── ...
│   └── processed/                               # created by the scripts
├── scripts/                                     # analysis + per-panel figure scripts (R)
├── figure_assembly/                             # multi-panel figure assembly (Python)
└── results/
    ├── figures/
    │   ├── composite/                           # final multi-panel figures (Fig 1–4, Fig S1–S6)
    │   └── PNG_hires/                           # individual rendered panels (300 dpi)
    └── tables/                                  # pre-computed result tables (CSV)
```

The scripts handle all preprocessing automatically (per-sample standard-10x directories, QC, integration). **No path editing is required** — see *Paths* below.

## Dependencies

- R ≥ 4.4 (tested on R 4.6.1)
- [Seurat](https://satijalab.org/seurat/) v5 (`Seurat`, `sctransform`, `glmGamPoi`)
- [harmony](https://github.com/immunogenomics/harmony)
- [UCell](https://github.com/carmonalab/UCell) (Bioconductor)
- [slingshot](https://github.com/kstreet13/slingshot) (Bioconductor)
- `mgcv` (GAM smoother, base R) and `segmented` (breakpoint regression) — used by `08_pseudotime_shape.R` and `09_robustness_lodo.R`
- `dplyr`, `tidyr`, `ggplot2`, `patchwork`
- Python ≥ 3.9 with `Pillow` (figure assembly only)

## Pipeline (run in order)

Scripts are numbered by **execution order**; comments in each script header map it back to the corresponding step of the parent study protocol (protocol steps 1, 2, 2b/sensitivity, 3, 6). Numbering follows the parent protocol, so **steps `06`–`07` are reserved** (DriftIndex formula freeze; permutation-based threshold) and are not yet part of this release — which is why the shape analysis is numbered `08`.

| Script | Analysis step | Input | Key output |
|---|---|---|---|
| `01_qc_cluster.R` | QC + SCTransform + Harmony integration + clustering + **bone-lineage extraction** | GSE169396 raw 10x (4 samples) | `seurat_gse169396_all.rds`, `seurat_gse169396_bone.rds`; QC/UMAP/DotPlot figures |
| `02_drift_analysis.R` | **Core drift analysis**: UCell scoring, cell-level DriftIndex, hybrid-state classification, per-sample ordering + age-dimension trends | `seurat_gse169396_bone.rds` | `seurat_gse169396_drift.rds`; `02_drift_summary.csv`; drift figures |
| `03_sensitivity.R` | Threshold-sweep sensitivity of the hybrid-state classification (informative thresholds ≥ 0) | `seurat_gse169396_drift.rds` | `03_sensitivity_summary.csv`, `03_sensitivity_trend.pdf` |
| `04_trajectory.R` | slingshot pseudotime; DriftIndex–pseudotime Spearman correlation; decile trend; lineage-gene dynamics | `seurat_gse169396_drift.rds` | `seurat_gse169396_traj.rds`; `04_trajectory_summary.csv`; trajectory figures |
| `05_cross_species.R` | Healthy mouse metaphysis processing + same scoring; cross-species DriftIndex & hybrid-proportion comparison | GSE317069 mouse 10x (5 samples) + `seurat_gse169396_drift.rds` | `seurat_mouse_metaphysis_bone.rds`; `05_cross_species_summary.csv`; cross-species figures |
| `08_pseudotime_shape.R` | **Shape of DriftIndex along pseudotime**: GAM smoother (`mgcv`) + single- and double-breakpoint segmented regression (`segmented`) + Davies test — establishes the biphasic (rise-then-decline) trajectory | `seurat_gse169396_traj.rds` | `08_pseudotime_shape_summary.csv`, `08_pseudotime_deciles.csv`, `08_shape_smoother.png`, `08_shape_deciles.png` |
| `09_robustness_lodo.R` | **Donor-aware robustness** of the DriftIndex–pseudotime association: leave-one-donor-out re-fitting, donor-cluster bootstrap, donor-block permutation test, and a within-donor-centred correlation — replaces the cell-level (pseudoreplicated) p-values with donor-level inference | `seurat_gse169396_traj.rds` | `09_robustness_summary.csv`, `09_lodo_summary.csv`, `09_lodo_rho_forest.png` |

```bash
Rscript scripts/01_qc_cluster.R
Rscript scripts/02_drift_analysis.R
Rscript scripts/03_sensitivity.R
Rscript scripts/04_trajectory.R
Rscript scripts/05_cross_species.R
Rscript scripts/08_pseudotime_shape.R
Rscript scripts/09_robustness_lodo.R
```

### Figure generation

Per-panel figures are rendered by the `fig*.R` scripts in `scripts/`, which write panels into `results/figures/PNG_hires/`. **Each panel has exactly one producing script** (no overlaps), and script names follow `fig<Figure><Panel>_<purpose>.R`:

| Figure script | Panels produced |
|---|---|
| `fig1_render_uniform.R` | Fig 1 UMAPs → `fig1_panels/{A1,A3,B}.png` |
| `fig_strata_rebuild.R` | Fig 2A/2B (by densitometric stratum) and Fig 4A/4B (by donor S1–S4) |
| `fig_scale_uniform.R` | Fig 2C age dimension and the Fig 3 panels; documents the shared axis-limit / colour / label conventions in its header |
| `figS1_scale_uniform.R`, `figS2_stratum_panel.R`, `figS3_scale_uniform.R`, `figS5_axis_uniform.R` | the corresponding Fig S1 / S2 / S3 / S5 panels |
| `fig_scale_diag.R` | helper only — prints axis ranges, renders nothing |

Multi-panel figures (Fig 1–4, Fig S1–S6) are then assembled by the two Python scripts in `figure_assembly/`:

```bash
python figure_assembly/compose_fig1_uniform.py     # Fig 1 (from fig1_panels/)
python figure_assembly/compose_rest_figures.py     # Fig 2–4 and Fig S1–S6 (from PNG_hires/)
```

Both scripts are idempotent and read only the panels listed above.

> The uniform axis limits, colour scheme, and label conventions applied across figures are documented in the header of `scripts/fig_scale_uniform.R`.

### Paths (nothing to edit — detected from the script location)

Every script is **portable out of the box**. Instead of a hard-coded absolute path, each one locates the project root *relative to its own file* (a script in `scripts/` or `figure_assembly/` looks one level up). They therefore work whether you invoke them from the repository root, from inside their own folder, or by absolute path:

```bash
Rscript scripts/04_trajectory.R                  # from the repo root
cd scripts && Rscript 04_trajectory.R            # from inside scripts/
python figure_assembly/compose_rest_figures.py   # from the repo root
```

Override any of them with an environment variable if your layout differs:

| Variable | Used by | Default |
|---|---|---|
| `DRIFT_PROJECT_DIR` | all R and Python scripts | repository root (inferred from the script's own location) |
| `DRIFT_PANELS_DIR` | Python assembly scripts | `<root>/results/figures/PNG_hires` |
| `DRIFT_COMPOSITE_DIR` | Python assembly scripts | `<root>/results/figures/composite` |

For example, to keep the outputs outside the checkout:

```bash
DRIFT_PROJECT_DIR=/data/drift Rscript scripts/04_trajectory.R
```

**Fonts** — the assembly scripts draw the panel letters with a bold sans-serif TTF. They probe a cross-platform candidate list in order (Arial Bold → Liberation Sans Bold → DejaVu Sans Bold) and fall back to Pillow's built-in bitmap font if none is present, so **no edit is needed on any OS**.

## Sample classification (human, GSE169396)

Per the original study's Supplementary Table 1 (Aging 2021, PMID 34111027), using lumbar spine T-scores and WHO criteria:

| Sample | Age/Sex | Lumbar T | Hip T | Densitometric stratum |
|---|---|---|---|---|
| S1 | 61 y / F | −3.0 | −1.9 | **Osteoporotic range** |
| S2 | 45 y / F | −1.3 | −1.2 | **Osteopenic range** |
| S3 | 66 y / M | NA | NA | Not stratifiable (no usable T-score) — descriptive only |
| S4 | 31 y / M | +0.6 | −1.1 | Not stratifiable (T-score not applied for diagnosis in men < 50 y) — descriptive only |

Only **S1 vs S2** supports a densitometric contrast; S3 and S4 are reported descriptively and no cross-group inference is drawn from them. All donors still contribute to the cell-level pseudotime and hybrid-state analyses.

## Key results tables (pre-computed, in `results/tables/`)

- `02_drift_summary.csv` — per-sample DriftIndex (mean ± SD), hybrid/osteo/adipo proportions, MD & SASP scores
- `03_sensitivity_summary.csv` — hybrid % across 7 thresholds
- `04_trajectory_summary.csv` — mean DriftIndex/AD/OS/MD per pseudotime decile
- `05_cross_species_summary.csv` — human vs healthy-mouse comparison
- `08_pseudotime_shape_summary.csv` — GAM fit, breakpoint estimates, pre/post-breakpoint slopes, Davies test
- `08_pseudotime_deciles.csv` — DriftIndex mean ± SD and cell counts per pseudotime decile
- `08_pseudotime_shape_report.txt` — full console report of the shape analysis
- `09_robustness_summary.csv` — donor-aware robustness metrics (leave-one-donor-out ρ range, donor-cluster bootstrap CI, donor-block permutation p, within-donor-centred ρ)
- `09_lodo_summary.csv` — per-fold leave-one-donor-out statistics (ρ, GAM edf, breakpoint, pre/post slopes, Davies p)
- `09_robustness_report.txt` — full console report of the donor-aware robustness analysis

## Figures (pre-rendered, in `results/figures/`)

**Final multi-panel figures** — `results/figures/composite/` (300 dpi):

| File | Content |
|---|---|
| `Fig1.png` | Study overview — DriftIndex UMAP, per-sample UMAP and cell-state composition (panels A–C) |
| `Fig2.png` | DriftIndex and cell-state proportions by densitometric stratum, plus the age dimension |
| `Fig3.png` | Osteoblast pseudotime — trajectory UMAP, the biphasic DriftIndex–pseudotime relation, lineage-gene dynamics |
| `Fig4.png` | Cross-species comparison (human femoral head vs healthy mouse metaphysis) |
| `FigS1.png` | QC violin plots (pre- and post-filter) |
| `FigS2.png` | Data integration — UMAP by sample / densitometric stratum / cluster, plus a lineage-marker DotPlot |
| `FigS3.png` | Threshold-sweep sensitivity of the hybrid-state classification |
| `FigS4.png` | Mouse metaphysis cell-state composition |
| `FigS5.png` | Per-donor DriftIndex and hybrid proportion for all four donors (descriptive) |
| `FigS6.png` | Shape of the DriftIndex along pseudotime — GAM smoother and segmented-regression fit |
| `Fig1_nolabel_preview.png` | Label-free preview of Fig 1 (layout check only) |

**Individual panels** — `results/figures/PNG_hires/` holds the 24 single-panel renders (300 dpi); the assembly scripts stitch most of them into the multi-panel figures above. Five panels (`01_UMAP_bone_lineage`, `02_UMAP_DriftIndex`, `02_UMAP_cell_state`, `02_hybrid_signature_heatmap`, `09_lodo_rho_forest`) are exploratory or analysis-only and are not stitched. Fig 1 is assembled from its own three panels in `results/figures/fig1_panels/` (`A1` DriftIndex, `A3` Sample, `B` Cell state).

## Methodological notes

- **DriftIndex** (cell level) = `UCell(AD_drift)` − `UCell(OS_identity)`; hybrid = both UCell scores > 0 (rank-based baseline).
- **Biphasic trajectory**: `08_pseudotime_shape.R` shows that the DriftIndex–pseudotime relationship is non-linear (GAM edf = 8.2; linear-vs-GAM likelihood-ratio p = 5.8 × 10⁻⁴⁵). A single-breakpoint segmented regression places the breakpoint at normalized pseudotime 0.448 (95% CI 0.390–0.506); the pre-breakpoint slope is positive (+0.187, p = 7.9 × 10⁻⁶) and the post-breakpoint slope is negative (−0.268, p = 2.3 × 10⁻⁴²) (Davies test p = 3.2 × 10⁻³¹). AIC further favours a two-breakpoint model, consistent with rise-then-decline.
- UCell scores are rank-normalized; **cross-species** comparisons therefore rely on hybrid proportions and relative trends, not absolute scores.
- **Donor-aware inference**: the 3,850 bone-lineage cells derive from only four donors, so cell-level p-values are pseudoreplicated. `09_robustness_lodo.R` therefore re-evaluates the DriftIndex–pseudotime association with the donor as the resampling unit — leave-one-donor-out (ρ range −0.33 to −0.28, all folds negative), donor-cluster bootstrap (95% CI −0.34 to −0.26), a donor-block permutation test that shuffles pseudotime within donors (p = 5 × 10⁻⁴, the resolution floor of 2,000 permutations), and a within-donor-centred correlation that removes between-donor confounding (ρ = −0.31, p = 5 × 10⁻⁴). The trajectory shape is equally stable across all four leave-one-donor-out folds (breakpoint 0.43–0.47).
- Cells within a donor are non-independent; per-sample summaries (mean ± SD) are reported alongside cell-level statistics. For the two-donor densitometric comparison the unit of inference (the donor) gives n = 1 per group, so that comparison is reported descriptively and without a p-value.

## License & citation

Data: GSE169396 (Aging 2021, PMID 34111027); GSE317069 (Nat Genet 2026, PMID 42432248). Please cite the original data papers and this repository if you reuse the code.

**License: All rights reserved.** This repository is released for academic review and reproducibility; reuse requires written permission from the corresponding author (see `LICENSE`).

---

*Analysis performed 2026-09. Contact: liqingqing_1999@163.com*
