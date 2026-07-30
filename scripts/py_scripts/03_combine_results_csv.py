# =======================================
# Combine LIANA+ Results .csv files
# =======================================
# For combining multiple targeted pair results .csvs (i.e all time points of an experiment)
# for easier interrogation

EXPERIMENT = "Xen1"
GROUPBY = "broad_celltype"

# ===========================

import pandas as pd
from pathlib import Path

files = sorted(
    Path("liana_python/results/tables").glob(
        f"{EXPERIMENT}_*_{GROUPBY}_supervisor_pairs_results.csv"
    )
)

combined = pd.concat(
    [pd.read_csv(f) for f in files],
    ignore_index=True
)

combined.to_csv(
    f"liana_python/results/tables/{EXPERIMENT}_{GROUPBY}_combined_results.csv",
    index=False
)