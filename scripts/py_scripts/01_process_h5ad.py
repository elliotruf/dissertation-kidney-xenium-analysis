# ==================================
# Prepare h5ad.py file for LIANA+
# ==================================
# Loads h5ad.py file, extracts necessary info for LIANA+,
# and generates a new prepared.h5ad file.

# ===================================
# User Settings
# ===================================

EXPERIMENT = "Xen1diet"
DATA = "Xen1diet.h5ad"  # h5ad.py file to be prepared

# ====================================
# Imports
# ====================================

import scanpy as sc
import liana as li
from pathlib import Path

# ====================================
# Project paths
# ====================================

PROJECT_DIR = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
OBJECT_DIR = RESULTS_DIR / "objects"

# ====================================
# Load data
# ====================================

print("Loading data...")
adata = sc.read_h5ad(DATA_DIR / DATA)

print("Preparing...")
adata.obsm["spatial"] = adata.obs[
    ["sdimx", "sdimy"]
].to_numpy()

adata.layers["counts"] = adata.X.copy()

# Save
print("Saving prepared .h5ad file...")
adata.write_h5ad(OBJECT_DIR / f"{EXPERIMENT}_prepared.h5ad")

# Done!
print("Done!")