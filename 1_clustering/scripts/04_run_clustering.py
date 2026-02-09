#!/usr/bin/env python3
# =============================================================================
# 04 - Run Clustering (k=3)
# =============================================================================
# Performs k-means clustering on IQ and SRS (validated in script 03).
# Outputs: cluster scatter plot, dataframe with cluster labels.
#
# Paper: Final clustering on IQ and SRS dimensions, k=3 clusters.
# =============================================================================

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

# --- Configuration ---
_script_dir = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.environ.get(
    "LEAP_INOVAND_DATA",
    os.path.join(_script_dir, "..", "..", "..", "imaging2genet", "0_input", "dataframes")
)
INDIVIDUALS_METRICS = os.path.join(DATA_PATH, "individuals_metrics.tsv")
OUTPUT_BASE = os.path.normpath(os.path.join(_script_dir, "..", "outputs"))
FIGURES_DIR = os.path.join(OUTPUT_BASE, "figures")
TABLES_DIR = os.path.join(OUTPUT_BASE, "tables")
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

N_CLUSTERS = 3
RANDOM_STATE = 42

# --- Data preparation ---
print(f"Loading data from: {INDIVIDUALS_METRICS}")
df = pd.read_csv(INDIVIDUALS_METRICS, sep="\t", low_memory=False)
df = df.drop_duplicates(subset=["ID"])
df = df.replace(999, np.nan)
df = df.replace(998, np.nan)
df = df[df["Relation_to_proposant"] == "participant"]

# IQ: total_IQ with performance_IQ fallback
df["IQ"] = df["total_IQ"].fillna(df["performance_IQ"])

clinical_features = ["IQ", "SRS_tscore"]
df_clust = df[clinical_features + ["ID", "Population1", "PopulationS1"]].dropna()
X = df_clust[clinical_features].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"Sample size for clustering: {len(df_clust)}")

# --- Clustering ---
kmeans = KMeans(n_clusters=N_CLUSTERS, random_state=RANDOM_STATE, n_init=25)
cluster_labels = kmeans.fit_predict(X_scaled)
df_clust = df_clust.copy()
df_clust["Cluster"] = [f"C{c + 1}" for c in cluster_labels]

# Merge cluster labels back to full dataset
df_with_clusters = df.merge(
    df_clust[["ID", "Cluster"]],
    on="ID",
    how="left"
)

# --- Scatter plot ---
palette = {"C1": "#324095", "C2": "#5CAEE1", "C3": "#C1C2BC"}
fig, ax = plt.subplots(figsize=(10, 8))
sns.scatterplot(
    data=df_clust,
    x="SRS_tscore",
    y="IQ",
    hue="Cluster",
    palette=palette,
    alpha=0.7,
    s=50,
    ax=ax
)
ax.set_xlabel("SRS-2 t-score")
ax.set_ylabel("Full-scale IQ")
ax.set_title("Clinical Clustering: IQ and SRS (k=3)")
ax.legend(title="Cluster")
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "clustering_scatter.pdf"), dpi=300, bbox_inches="tight")
plt.close()
print(f"Figure saved: clustering_scatter.pdf")

# --- Save cluster assignments ---
df_clust_out = df_clust[["ID", "IQ", "SRS_tscore", "Cluster", "Population1", "PopulationS1"]]
df_clust_out.to_csv(os.path.join(TABLES_DIR, "cluster_assignments.csv"), index=False)
df_with_clusters.to_csv(os.path.join(TABLES_DIR, "individuals_metrics_with_clusters.csv"), index=False, sep="\t")
print(f"Cluster assignments saved: cluster_assignments.csv")
print(f"Full dataset with clusters saved: individuals_metrics_with_clusters.csv")

# --- Cluster summary ---
print("\nCluster sizes:")
print(df_clust["Cluster"].value_counts().sort_index())
print("\nDone. Outputs in:", OUTPUT_BASE)
