# LEAP-InovAND Resource

This repository contains analysis code for our last paper: LEAP-InovAND a multiscale resource to explore genetics, brain imaging and clinical data in autism.

## Overview

Analysis is organized by paper section. Each analysis folder contains:
- `scripts/` — R and Python scripts
- `outputs/` — generated figures and tables
  - `figures/` — PDF figures for publication
  - `tables/` — CSV/TSV tables

## 1. Clustering Analysis

Clinical clustering based on IQ and SRS-2 dimensions, following the methodology described in the paper:

- **PCA** on a broad set of clinical measures (IQ subscales, SRS-2, Vineland, RBS-R, SSP) to identify primary variance components
- **Feature selection** rationale: IQ and SRS maximize sample size (n=1,026) and contribute most to PC1/PC2
- **Cluster validation** (k=2–10): NbClust (26 indices), AIC/BIC/ICL, Silhouette, VRS, Davies-Bouldin, Pseudo-F, Gap statistic
- **Final clustering**: k=3 clusters on IQ and SRS

### Scripts (run in order)

```bash
cd 1_clustering/scripts

# 1. PCA on clinical features
Rscript 01_pca_features.R

# 2. Feature selection rationale (sample sizes, PCA loadings)
Rscript 02_feature_selection_rationale.R

# 3. Cluster number validation (NbClust, Gap, etc.)
Rscript 03_cluster_validation.R

# 4. Run final clustering (k=3)
python 04_run_clustering.py
```

### Outputs

| Output | Description |
|-------|-------------|
| `figures/PCA_variance.pdf` | Scree plot |
| `figures/PCA_loadings.pdf` | Variable loadings (Supplementary Fig. 21 / S20) |
| `figures/PCA_contrib_PC1.pdf`, `PCA_contrib_PC2.pdf` | Variable contributions |
| `figures/sample_size_completeness.pdf` | Sample size by variable combination |
| `figures/gap_statistic.pdf` | Gap statistic for k selection |
| `figures/NbClust_consensus.pdf` | NbClust index consensus |
| `figures/silhouette_k3.pdf` | Silhouette plot for k=3 |
| `figures/clustering_scatter.pdf` | IQ × SRS scatter with clusters |
| `tables/PCA_variable_importance.csv` | PCA loadings and contributions |
| `tables/sample_sizes_by_variable_combination.csv` | Supplementary Table 12 |
| `tables/cluster_validation_indices.csv` | Supplementary Table 13 |
| `tables/cluster_assignments.csv` | Individual cluster labels |

### Data

Scripts expect `individuals_metrics.tsv` in the data directory. Set the path via:

```bash
export LEAP_INOVAND_DATA="/path/to/data/directory"
```

Default path (relative to repository): `../imaging2genet/0_input/dataframes/`

### R Dependencies

```r
install.packages(c("readr", "FactoMineR", "factoextra", "NbClust", "mclust", "cluster", "fpc", "ggplot2", "corrplot", "dplyr", "gt"))
```

### Python Dependencies

```
pandas numpy matplotlib seaborn scikit-learn
```

## Citation

If you use this resource, please cite:

> LEAP-InovAND: A multiscale resource to explore genetics, brain imaging, and clinical data in autism.  
> medRxiv preprint (2025). https://doi.org/10.1101/2025.11.24.25340858
