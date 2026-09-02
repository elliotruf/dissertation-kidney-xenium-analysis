# =======================================
# Combine LIANA+ Results .csv files
# =======================================
# For combining multiple targeted pair results .csvs (i.e all time points of an experiment)
# for easier interrogation

EXPERIMENT = "Xen1"
GROUPBY = "seurat_clusters"

# ====================================
# Imports
# ====================================

import pandas as pd
from pathlib import Path

# ====================================
# Project paths
# ====================================

PROJECT_DIR = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
OBJECTS_DIR = RESULTS_DIR / "objects"
OUTPUTS_DIR = RESULTS_DIR / "liana_cell_cell_communication"
TABLES_DIR = OUTPUTS_DIR / "tables"
FIGURES_DIR = OUTPUTS_DIR / "figures"

# ====================================

files = sorted(
    Path(TABLES_DIR).glob(
        f"{EXPERIMENT}_*_{GROUPBY}_results_sender_receiver_targeted_allpairs.csv"
    )
)

combined = pd.concat(
    [pd.read_csv(f) for f in files],
    ignore_index=True
)

# Save
print("Saving combined results...")
combined.to_csv(
    TABLES_DIR / 
    f"{EXPERIMENT}_{GROUPBY}_combined_results_sender_receiver_targeted_allpairs.csv",
    index=False
)

# Done!
print("Done!")