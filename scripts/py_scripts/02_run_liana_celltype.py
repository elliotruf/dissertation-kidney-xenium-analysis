# Liana+ step one
# Xen1: Week 1 only
# Broad cell types
# Non spatial

# ===================
# User Settings
# ===================

EXPERIMENT = "Xen1"
TIMEPOINT = "1wk"
GROUPBY = "broad_celltype"

DATA = "liana_python/data/Xen1diet_prepared.h5ad"   # .h5ad file ready for analysis
PAIR_FILE = "info/Xenium_Probe_List_May2025_3_clean.csv" # curated ligand/receptor pair list 

OUTDIR = "liana_python/results/tables"

# ===================

import scanpy as sc
import liana as li
import pandas as pd

# Load DATA
adata = sc.read_h5ad(DATA)

# Subset by TIMEPOINT
adata_tp = adata[
    adata.obs["time_point"] == TIMEPOINT
].copy()

# Use normalised expression
adata.X = adata.layers["logcounts"].copy()

# Run LIANA+
li.mt.rank_aggregate(
    adata_tp,
    groupby=GROUPBY,
    resource_name="mouseconsensus",
    expr_prop=0.1,                      # Requires ligand/receptor to be expressed >10% of cells in a cell type
    use_raw=False,
    inplace=True
)

# Subset by PAIR_FILE
results = adata_tp.uns["liana_res"]

pairs = pd.read_csv(PAIR_FILE)

# Remove Excel row-number column
pairs = pairs.drop(columns="Unnamed: 0")

# Standardise names
pairs = pairs.rename(
    columns={
        "Ligand": "ligand_complex",
        "Receptor": "receptor_complex"
    }
)

# Add interaction label for interpretability
pairs["interaction"] = (
    pairs["ligand_complex"]
    + " / "
    + pairs["receptor_complex"]
)

# Select for supervisor pairs
supervisor_pairs = results.merge(
    pairs,
    on=["ligand_complex", "receptor_complex"],
    how="inner"
)

supervisor_pairs["experiment"] = EXPERIMENT
supervisor_pairs["time_point"] = TIMEPOINT
supervisor_pairs["annotation"] = GROUPBY

# Save .CSV files
results.to_csv(                                                             # All results
    f"{OUTDIR}/{EXPERIMENT}_{TIMEPOINT}_{GROUPBY}_results.csv",
    index=False
)

supervisor_pairs.to_csv(
    f"{OUTDIR}/{EXPERIMENT}_{TIMEPOINT}_{GROUPBY}_supervisor_pairs_results.csv",     # Results of supervisor pairs
        index=False
    )

# Done!
print("Done!")