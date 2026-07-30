# For visualising combined_results
# Includes:
# Heatmap of interaction strength over time
# Heatmap of source -> target communication
# Line plots of pathways of interest

# ============================================
# User Settings
# ============================================

EXPERIMENT = "Xen1"                               # Experiment name for naming files
GROUPBY = "broad_celltype"                        # Broad / fine cell type 
TIME_POINT_ORDER = ["1wk", "2wk", "4wk", "12wk"]  # Order to display timepoints

DATA = "Xen1_broad_celltype_combined_results.csv" # Name of combined results file to analyse

# ====================================
# Imports
# ====================================

import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.stats import zscore

# ====================================
# Project paths
# ====================================

PROJECT_DIR = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_DIR / "data"
PAIR_DIR = PROJECT_DIR / "info"
RESULTS_DIR = PROJECT_DIR / "results" / "liana_cell_cell_communication"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"

combined = pd.read_csv(TABLES_DIR / f"{EXPERIMENT}_{GROUPBY}_combined_results.csv")

heatmap_df = (
    combined
    .groupby(["interaction", "time_point"])["lrscore"]
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

# Mark NAs
heatmap_df = heatmap_df.fillna(0)

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
    FIGURES_DIR / f"{EXPERIMENT}_{GROUPBY}_heatmap_raw.pdf",
    dpi=600,
    bbox_inches="tight"
)

plt.close()

# ===================================
# Figure 2: Row-scaled heatmap
# ===================================

# Adjust scale
heatmap_scaled = heatmap_df.copy()

heatmap_scaled = heatmap_scaled.sub(
    heatmap_scaled.mean(axis=1),
    axis=0
)

heatmap_scaled = heatmap_scaled.div(
    heatmap_scaled.std(axis=1).replace(0,1),
    axis=0
)

lim  = abs(heatmap_scaled.values).max()

# Cluster by z-score
cluster = sns.clustermap(
    heatmap_scaled,
    row_cluster=True,
    col_cluster=False
)

# Plot scaled heatmap
ordered = heatmap_scaled.index[cluster.dendrogram_row.reordered_ind]

heatmap_scaled = heatmap_scaled.loc[ordered]

plt.figure(figsize = (12,10))

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
plt.title("Relative ligand-receptor interaction dynamics across kidney development")
plt.xlabel("Developmental stage")
plt.ylabel("Ligand-receptor pair")

plt.tight_layout()

# Save scaled heatmap
print("Saving relative interaction heatmap...")
plt.savefig(
    FIGURES_DIR / f"{EXPERIMENT}_{GROUPBY}_heatmap_scaled.pdf",
    dpi=600,
    bbox_inches="tight"
)

# Done!
print("Done!")