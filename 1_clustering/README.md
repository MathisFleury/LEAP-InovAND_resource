# 1. Clustering Analysis

Clinical clustering on IQ and SRS-2 dimensions for the LEAP-InovAND Nature Neuroscience paper.

## Methodology Summary

1. **Broad clinical measures**: IQ subscales (full-scale, verbal, non-verbal), SRS-2, Vineland (composite + domains), RBS-R, SSP
2. **PCA**: Exclude missing data and relatives; standardize features. First two components explain >75% variance; IQ and SRS contribute most; Vineland to PC1 but minimally to PC2 (Supplementary Fig. 21)
3. **Feature selection**: IQ and SRS maximize sample size (n=1,026 vs ~500 with VABS) — Supplementary Table 12
4. **Cluster validation** (k=2–10): NbClust (26 indices) → k=3 (10/26); secondary peak at k=8. Additional indices (AIC, BIC, ICL, Silhouette, VRS, Davies-Bouldin, Pseudo-F, Gap) favor k=3 — Supplementary Table 13
5. **Final clustering**: k=3 on IQ and SRS

## Scripts

| Script | Description |
|--------|-------------|
| `01_pca_features.R` | PCA on all clinical measures; scree plot, loadings, variable contributions |
| `02_feature_selection_rationale.R` | Sample size completeness, PCA loadings table |
| `03_cluster_validation.R` | NbClust + model-based (AIC/BIC/ICL) + partition (Silhouette, VRS, Davies-Bouldin, Pseudo-F) + Gap |
| `04_run_clustering.py` | K-means k=3 on IQ and SRS; scatter plot, cluster assignments |

## Run

```bash
# From 1_clustering/
./run_all.sh

# Or individually from 1_clustering/scripts/
Rscript 01_pca_features.R
Rscript 02_feature_selection_rationale.R
Rscript 03_cluster_validation.R
python 04_run_clustering.py
```

## Data

Requires `individuals_metrics.tsv` with columns: `ID`, `Relation_to_proposant`, `total_IQ`, `performance_IQ`, `SRS_tscore`, and optionally `ssp_total`, `RBS-R_total`, `vabsdscoresc_dss`, `vabsdscoresd_dss`, `vabsdscoress_dss`, `vabsabcabc_standard` for full PCA.

Set path: `export LEAP_INOVAND_DATA="/path/to/data/dir"`
