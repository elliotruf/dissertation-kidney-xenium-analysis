# For visualising combined_results
# Includes:
# Heatmap of interaction strength over time
# Heatmap of source -> target communication
# Line plots of pathways of interest

# ============================================
# User Settings
# ============================================

EXPERIMENT = "Xen1"                               # Experiment name for naming files
GROUPBY = "seurat_cluster"                        # Broad / fine cell type / seurat clusters
TIME_POINT_ORDER = ["1wk", "2wk", "4wk", "12wk"]  # Order to display timepoints

ANN_DATA = "Xen1diet_prepared.h5ad"
DATA = "Xen1_seurat_clusters_combined_results_sender_receiver_targeted_allpairs.csv" # Name of combined results file to analyse

# ====================================
# Imports
# ====================================

import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
from pathlib import Path
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
combined["source"] = combined["source"].astype(str)
combined["target"] = combined["target"].astype(str)

combined["interaction"] = (
    combined["ligand_complex"].astype(str)
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
# Relative pair interactions

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

# ============================================
# Figure 3: Broad cell-type LR pair distribution
# ============================================

broad_distribution = (
    combined
    .groupby(
        ["time_point", "source", "target"]
    )["interaction"]
    .nunique()
    .reset_index(name="n_lr_pairs")
)

print("\nBroad distribution:")
print(broad_distribution.head(20))

fig, axes = plt.subplots(
    2, 2,
    figsize=(16, 14)
)

axes = axes.flatten()

max_lr_pairs = broad_distribution["n_lr_pairs"].max()

for ax, timepoint in zip(
    axes,
    TIME_POINT_ORDER
):

    plot_df = (
        broad_distribution[
            broad_distribution["time_point"] == timepoint
        ]
        .pivot(
            index="source",
            columns="target",
            values="n_lr_pairs"
        )
        .fillna(0)
    )

    print(
        f"{timepoint}:",
        plot_df.shape,
        "non-NA values =",
        plot_df.notna().sum().sum()
    )

    sns.heatmap(
        plot_df,
        ax=ax,
        cmap="viridis",
        vmin=0,
        vmax=max_lr_pairs,
        linewidths=0,
        cbar=True,
        cbar_kws={
            "label": "Number of targeted LR pairs"
        }
    )

    ax.set_title(timepoint)
    ax.set_xlabel("Target Seurat cluster")
    ax.set_ylabel("Source Seurat cluster")

    ax.tick_params(
        axis="x",
        rotation=45
    )

    ax.tick_params(
        axis="y",
        rotation=0
    )

fig.suptitle(
    "Distribution of targeted ligand-receptor pairs between Seurat clusters across Xen1 kidney development",
    fontsize=16
)

plt.tight_layout()

print("Saving broad cell-type LR pair distribution heatmap...")

plt.savefig(
    FIGURES_DIR /
    f"{EXPERIMENT}_broad_celltype_LR_pair_distribution.pdf",
    dpi=600,
    bbox_inches="tight"
)

plt.close()

print("Broad cell-type LR pair distribution heatmap saved.")