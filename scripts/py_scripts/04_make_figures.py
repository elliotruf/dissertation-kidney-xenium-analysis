# For visualising combined_results
# Includes:
# Heatmap of interaction strength over time
# Heatmap of source -> target communication
# Line plots of pathways of interest

# ============================================
# User Settings
# ============================================

EXPERIMENT = "Xen1"                               # Experiment name for naming files
GROUPBY = "seurat_clusters"                        # Broad / fine cell type / seurat clusters
TIME_POINT_ORDER = ["1wk", "2wk", "4wk", "12wk"]  # Order to display timepoints

ANN_DATA = "Xen1diet_prepared.h5ad"
DATA = "Xen1_seurat_clusters_combined_results_sender_receiver_targeted_allpairs.csv" # Name of combined results file to analyse

# ====================================
# Imports
# ====================================

import pandas as pd
import seaborn as sns
import scanpy as sc
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.stats import zscore

# ====================================
# Project paths
# ====================================

PROJECT_DIR = Path(__name__).resolve().parents[0]

DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
OBJECTS_DIR = RESULTS_DIR / "objects"
OUTPUTS_DIR = RESULTS_DIR / "liana_cell_cell_communication"
TABLES_DIR = OUTPUTS_DIR / "tables"
FIGURES_DIR = OUTPUTS_DIR / "figures"

adata = sc.read_h5ad(OBJECTS_DIR / ANN_DATA)

combined = pd.read_csv(TABLES_DIR / DATA)

# Create interaction label containing
# sender, receiver, ligand and receptor
cluster_lookup = (
    adata.obs[["seurat_clusters", "cell_type"]]
    .drop_duplicates()
)
combined["source"] = combined["source"].astype(str)
combined["target"] = combined["target"].astype(str)

cluster_lookup["seurat_clusters"] = (
    cluster_lookup["seurat_clusters"].astype(str)
)

cluster_to_celltype = dict(
    zip(
        cluster_lookup["seurat_clusters"],
        cluster_lookup["cell_type"]
    )
)

combined["source_celltype"] = (
    combined["source"].map(cluster_to_celltype)
)

combined["target_celltype"] = (
    combined["target"].map(cluster_to_celltype)
)

combined["interaction"] = (
    combined["source_celltype"]
    + " → "
    + combined["target_celltype"]
    + " : "
    + combined["ligand_complex"].astype(str)
    + " → "
    + combined["receptor_complex"].astype(str)
)

heatmap_df = (
    combined
    .groupby(
        ["interaction", "time_point"]
    )["lrscore"]
    .max()
    .reset_index()
)

heatmap_df = heatmap_df.pivot(
    index="interaction",
    columns="time_point",
    values="lrscore"
)

# Order time points
time_point_order = TIME_POINT_ORDER

heatmap_df = heatmap_df.reindex(columns=time_point_order)

# Keep only interactions with scores at every timepoint
heatmap_complete = heatmap_df.dropna(
    subset=TIME_POINT_ORDER
).copy()

print("Total interactions:", len(heatmap_df))
print("Complete interactions:", len(heatmap_complete))

# =================================
# Figure 1: RAW LR score heatmap
# =================================
plt.figure(figsize=(12,10))

sns.heatmap(
    heatmap_df,
    cmap="viridis",
    linewidths=0,
    cbar_kws={"label": "LIANA LR score"}
)

# Aesthetic adjustments
plt.xlabel("Developmental stage")
plt.ylabel("Ligand-receptor pair")
plt.title("Absolute ligand-receptor interaction dynamics across kidney development")

plt.tight_layout()

# Save raw heatmap
print("Saving raw LR score heatmap...")
plt.savefig(
    FIGURES_DIR / f"{EXPERIMENT}_{GROUPBY}_heatmap_raw_allpairs.pdf",
    dpi=600,
    bbox_inches="tight"
)

plt.close()

# ===================================
# Figure 2: Row-scaled heatmap
# ===================================

# Start with a copy of the raw interaction score matrix
heatmap_scaled = heatmap_df.copy()

# Remove infinite values
heatmap_scaled = heatmap_scaled.replace(
    [np.inf, -np.inf],
    np.nan
)

# Remove interactions containing missing values
heatmap_scaled = heatmap_scaled.dropna(
    axis=0,
    how="any"
)

print(
    "Interactions remaining before scaling:",
    heatmap_scaled.shape[0]
)

# Calculate row standard deviations
row_sd = heatmap_scaled.std(axis=1)

# Remove rows with no variation across conditions
heatmap_scaled = heatmap_scaled.loc[
    row_sd > 0
]

print(
    "Interactions remaining after removing zero-variance rows:",
    heatmap_scaled.shape[0]
)

# Row-scale interaction scores
heatmap_scaled = heatmap_scaled.sub(
    heatmap_scaled.mean(axis=1),
    axis=0
).div(
    heatmap_scaled.std(axis=1),
    axis=0
)

# Final safety check
if not np.isfinite(heatmap_scaled.to_numpy()).all():
    raise ValueError(
        "Scaled heatmap still contains non-finite values."
    )

lim = 2

# Cluster interactions by temporal patterns
cluster = sns.clustermap(
    heatmap_scaled,
    row_cluster=True,
    col_cluster=False,
    metric="correlation",
    method="average"
)

# Extract clustered row order
ordered = heatmap_scaled.index[
    cluster.dendrogram_row.reordered_ind
]

heatmap_scaled = heatmap_scaled.loc[ordered]

# Close the clustermap figure, as it is only being
# used to obtain the clustering order
plt.close(cluster.fig)

# Plot scaled heatmap
plt.figure(figsize=(12, 10))

sns.heatmap(
    heatmap_scaled,
    cmap="RdBu_r",
    center=0,
    vmin=-lim,
    vmax=lim,
    linewidths=0,
    cbar_kws={"label": "Row z-score"}
)

# Aesthetic adjustments
plt.title(
    "Relative ligand-receptor interaction dynamics across kidney development"
)
plt.xlabel("Developmental stage")
plt.ylabel("Ligand-receptor pair")

plt.tight_layout()

# Save scaled heatmap
print("Saving relative interaction heatmap...")

plt.savefig(
    FIGURES_DIR / f"{EXPERIMENT}_{GROUPBY}_heatmap_scaled_allpairs.pdf",
    dpi=600,
    bbox_inches="tight"
)

plt.close()

# Done!
print("Done!")