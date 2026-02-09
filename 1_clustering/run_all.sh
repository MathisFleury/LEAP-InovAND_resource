#!/bin/bash
# Run all clustering analysis scripts in order
set -e
cd "$(dirname "$0")/scripts"

echo "=== 1. PCA on clinical features ==="
Rscript 01_pca_features.R

echo "=== 2. Feature selection rationale ==="
Rscript 02_feature_selection_rationale.R

echo "=== 3. Cluster validation ==="
Rscript 03_cluster_validation.R

echo "=== 4. Run clustering (k=3) ==="
python 04_run_clustering.py

echo "=== Done. Outputs in 1_clustering/outputs/ ==="
