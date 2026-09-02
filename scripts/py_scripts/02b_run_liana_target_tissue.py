# =====================
# Run LIANA+
# =====================
# For running non-spatial cell-cell communication inference 
# on a prepared.h5ad file using LIANA+

# ===================
# User Settings
# ===================

EXPERIMENT = "Xen1"                 # Experiment name for naming files
DATA = "Xen1diet_prepared.h5ad"     # Prepared .h5ad file for analysis
TIMEPOINT = "12wk"                   # Time point for analysis
GROUPBY = "seurat_clusters"          # Fine/broad cell type annotations

# Analysis Parameters
SENDER_RECEIVER_FILE = "info/xen1_cell_chat_pairs_2.csv" # Targeted ligand-receptor pair list

# ====================================
# Imports
# ====================================

import scanpy as sc
import liana as li
from pathlib import Path
import pandas as pd

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

# ==================
# Load data
# ==================

print("Loading data...")

adata = sc.read_h5ad(OBJECTS_DIR / DATA)

# Subset by TIMEPOINT
print("Creating timepoint subset...")

adata_tp = adata[
    adata.obs["time_point"] == TIMEPOINT
].copy()

# Use normalised expression
adata.X = adata.layers["logcounts"].copy()

# =========================
# Run LIANA+
# =========================

print("Running LIANA+...")

li.mt.rank_aggregate(
    adata_tp,
    groupby=GROUPBY,
    resource_name="mouseconsensus",
    expr_prop=0.1,                      # Requires ligand/receptor to be expressed >10% of cells in a cell type
    use_raw=False,
    inplace=True
)

# Load results
results = adata_tp.uns["liana_res"]

# Load targeted sender/receiver pairs
pairs = pd.read_csv(
    PROJECT_DIR / SENDER_RECEIVER_FILE,
    dtype={"source": str, "target": str}
)

print("Targeted sender-receiver pairs:")
print(pairs)

# Select for targeted sender/receiver pairs
targeted_pairs = results.merge(
    pairs,
    on=["source", "target"],
    how="inner"
)

targeted_pairs["experiment"] = EXPERIMENT
targeted_pairs["time_point"] = TIMEPOINT
targeted_pairs["annotation"] = GROUPBY

# Save .CSV files
print("Saving results...")
results.to_csv(
    TABLES_DIR / f"{EXPERIMENT}_{TIMEPOINT}_{GROUPBY}_results_sender_receiver_all.csv",                   # All results
    index=False
)

targeted_pairs.to_csv(
    TABLES_DIR / f"{EXPERIMENT}_{TIMEPOINT}_{GROUPBY}_results_sender_receiver_targeted.csv",    # Results of supervisor pairs
        index=False
    )

# Done!
print("Done!")