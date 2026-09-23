
# Decoding Kidney Microenvironments using Xenium Spatial Transcriptomics

## Project Overview

This project analyses murine kidney Xenium spatial transcriptomic
datasets to investigate cellular and molecular changes throughout
development as well as injury-recovery following stromal GATA3 knockout.
The repository contains workflows for Region of Interest (ROI)
extraction, pseudobulk differential expression, whole-kidney cell-type
differential expression, pathway enrichment and visualisation.

## Repository Structure

``` text
scripts/
├── py_scripts/
|   ├── 01_process_h5ad.py
|   ├── 02_run_liana.py
|   ├── 03_combine_results.py
|   └── 04_make_figures.py
└── R_scripts/
    ├── helpers/
    ├── 01a_ROI_extraction.R
    ├── 01b_sex_pool.R
    ├── 02_ROI_pseudobulk.R
    ├── 03_ROI_broad_celltype_pseudobulk.R
    ├── 04_wholekidney_celltype.R
    ├── 05_cell_composition.R
    └── 06_convert_seurat_to_h5ad.R

results/
├── objects/
├── ROI_pseudobulk/
├── ROI_broad_celltype_pseudobulk/
├── ROI_singlecell_mixed_model/
├── whole_kidney_celltype/
├── cell_composition/
└── liana_cell_cell_communication/
```

`scripts/` contains all R and Python scripts required for analyses and
visualisations. Helper R functions are stored under `scripts/helpers/`
and grouped by purpose.

`results/` is the output directory for all analysis outputs, including
tables, tibbles, and visualisations. `.Rds` objects can be found in
`objects/`, while analysis outputs for each analysis pipeline can be
found in their corresponding subdirectory in `results/`.

## Data Requirements

The project expects the following directory format.

``` text
[Project]/
├── data/
|   ├── [dataset].rds
|   ├── [dataset].h5ad
|   ├── [dataset_group_tissue]/
|   |   ├── [time point]/
|   |   |   ├── [ROI_1_cells_stats].csv
|   |   |   ├── [ROI_2_cells_stats].csv
|   |   |   └── ...
|   |   ├── [time point]/
|   |   |   ├── [ROI_1_cells_stats].csv
|   |   |   ├── [ROI_2_cells_stats].csv
|   |   |   └── ...
|   |   └── ...
|   └── [targeted_ligand_receptor_pairs].csv
├── scripts/
├── results/
└── README.Rmd
```

### Directory descriptions

| Item | Description |
|----|----|
| **`[dataset].rds`** | Whole-dataset Seurat object (e.g. `Xen1diet.rds`). |
| **`[dataset].h5ad`** | AnnData version of the same dataset used for LIANA+ (generated automatically by `06_convert_seurat_to_h5ad.R`). |
| **`[dataset_group_tissue]/`** | Folder containing ROI exports for a particular tissue, region, or experimental group (e.g. `Xen1_KO_Vessels`). |
| **`[time_point]`** | Subdirectory containing ROI exports for a single developmental or experimental time point (e.g. `wk1`, `wk4`, `Sham`). Folder names are user-defined but must match values supplied in the script USER SETTINGS. |
| **`[ROI_X_cells_stats.csv]`** | `cells_stats.csv` files exported directly from Xenium Explorer. File names are not important; each file represents one ROI. |
| **`[targeted_ligand_receptor_pairs.csv]`** | Optional list of ligand-receptor pairs used for targeted LIANA+ analyses. |

### Seurat Object requirements

The Seurat object should contain:

- Cell-level Xenium transcript counts stored in the default assay
- Cell metadata including:
  - developmental or experimental time point,
  - biological replicate/sample identifier,
  - fine cell-type annotation,
  - broad cell-type annotation
  - any additional metadata required for grouping or downstream
    analyses.
    - examples include sex (M/F) or genotype stratification (wild type /
      knockout).

For LIANA+ analysis, the Seurat object must first be converted to
AnnData (.h5ad) format using `06_convert_seurat_to_h5ad.R`.

### Dependencies

#### R

Required CRAN packages:

``` r
# CRAN packages
cran_packages <- c(
  "Seurat",
  "tidyverse",
  "readxl",
  "fs",
  "ggplot2",
  "scales",
  "Matrix",
  "ggrepel",
  "pheatmap",
  "stringr",
  "presto"
)

# Install CRAN packages if needed
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
```

Required Bioconductor packages:

``` r
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Required Bioconductor packages
bioc_packages <- c(
  "edgeR",
  "clusterProfiler",
  "org.Mm.eg.db",
  "enrichplot",
  "GO.db"
)

# Install Bioconductor packages if needed
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}
```

#### Python

LIANA+ analyses require a Conda environment.

``` bash
conda activate liana
```

Required Python packages include: - scanpy - anndata - liana -
decoupler - pandas - seaborn - matplotlib - path

## Analysis Workflow

The overall analysis proceeds through five major stages:

1.  ROI extraction from ROIs exported from Xenium Explorer.
2.  ROI-based differential expression using `edgeR` and enrichment
    analyses using `clusterProfiler`.
3.  Whole-kidney cell-type differential expression using `Seurat` and
    enrichment analyses using `clusterProfiler`.
4.  Cell composition analysis.
5.  Cell-cell communication inference using `LIANA+`.

\[ DIAGRAM \]

## Pipeline Components

`01a_ROI_extraction.R` – Extracts cell IDs from manually selected
Regions of Interest (ROIs) in \``_cells_stats.csv` format exported from
Xenium Explorer, matches them to a whole-dataset Seurat object, and
generates a list of ROI-specific Seurat objects for downstream analyses.

`01b_sex_pool.R` – Combines separate “Male” and “Female” Seurat objects
into a single pooled “ALL” object. This is particularly useful for
analyses of vascular ROIs, where individual sexes may contain
insufficient cell numbers for downstream analyses.

`02_ROI_pseudobulk.R` – Performs pseudobulk differential expression
analysis and Gene Ontology enrichment analysis for each ROI,
enabling comparison of transcriptional changes within ROIs across
developmental or experimental time points.

`03_ROI_broad_celltype_pseudobulk.R` – Performs pseudobulk differential
expression analysis and Gene Ontology enrichment analysis on a specified
broad cell type (e.g. “Stroma”, “Proximal Tubule”) within each ROI,
enabling comparison of cell-type-specific transcriptional changes across
time points.

`04_wholekidney_celltype_analysis.R` – Performs whole-kidney
differential expression analysis and Gene Ontology enrichment analysis
for each broad cell type (e.g. all Stromal cells or all Proximal Tubule
cells), allowing assessment of transcriptional changes across
developmental or experimental time points.

`05_cell_composition.R` – Calculates and visualises changes in cell-type
composition across developmental or experimental time points using
stacked bar charts and summary tables.

`convert_seurat_to_h5ad.R` – Converts a Seurat object into the `.h5ad`
format required for the Python-based LIANA+ cell-cell communication
inference workflow.

### User Settings

All scripts contain a USER SETTINGS section near the beginning of the
file, allowing key analysis parameters to be modified without editing
the main body of the script. Typical options include:

- dataset selection
- developmental or experimental time point,
- grouping variables,
- statistical thresholds (e.g. adjusted P-value and log2 fold-change),
- Gene Ontology ontology selection,
- output locations.

Available options are tailored to the requirements of each individual
script. Any changes should be saved before sourcing or running the
script.

### Step 1. ROI Extraction

Following ROI selection in Xenium Explorer, export the ROI files in the
expected format and run the ROI extraction workflow:

``` r
source("scripts/R_scripts/01a_ROI_extraction.R")
```

For analyses where individual datasets contain relatively few cells (for
example, vascular ROIs), male and female datasets can optionally be
pooled to increase statistical power:

``` r
source("scripts/R_scripts/01b_sex_pool.R")
```

### Step 2. ROI Analysis

Following ROI extraction, analysis can be performed.

#### ROI pseudobulk differential expression:

``` r
source("scripts/R_scripts/02_ROI_pseudobulk.R")
```

Each ROI is treated as an independent biological replicate. Gene counts
are aggregated into pseudobulk profiles before differential expression
and Gene Ontology enrichment analyses are performed. Pseudobulk DE is
performed using EdgeR (Chen et al. 2025). The workflow automatically
generates:

- differential expression tables,
- volcano plots,
- Gene Ontology enrichment anlayses,
- GO dot plot visualisations.

#### ROI broad cell-type pseudobulk analyses:

``` r
source("scripts/R_scripts/03_ROI_broad_celltype_pseudobulk.R")
```

This workflow first divides each ROI according to broad cell type (for
example, Stroma, Proximal Tubule, Immune, Endothelium, etc.) before
performing pseudobulk differential expression independently for each
cell type. Pseudobulk DE is performed using EdgeR (Chen et al. 2025).

Analyses are skipped automatically if:

- fewer than four (4) ROI pseudobulks are available,
- fewer than two (2) biological replicates are present,
- too few genes remain following filtering.

The workflow automatically generates:

- differential expression tables,
- volcano plots,
- Gene Ontology enrichment analyses,
- GO dot plot visualisations.

### Step 3. Whole-kidney cell-type Analysis

``` r
source("scripts/R_scripts/04_wholekidney_celltype_analysis.R")
```

Rather than analysing individual ROIs, this workflow analyses every cell
of a specified broad cell type across the entire kidney sample. Outputs
include:

- differential expression tables,
- stacked bar chart summaries,
- heatmaps of differentially expressed genes,
- Gene Ontology enrichment analyses,
- GO dot plot visualisations.

### Step 4. Cell Composition

``` r
source("scripts/R_scripts/05_cell_composition.R")
```

This workflow quantifies changes in cell-type composition across
developmental or experimental time points and generates:

- summary tables,
- stacked bar chart visualisations.

### Step 5. Cell-Cell Communication Inference using LIANA+

Cell-cell communication analysis is performed in Python using LIANA+
(Dimitrov et al. 2024).

First convert the dataset’s Seurat object to AnnData (.h5ad) format:

``` r
source("scripts/R_scripts/06_convert_Seurat_to_h5ad.R")
```

In your machine’s command terminal, navigate to the Project Directory.
Then activate the LIANA Conda environment.

``` bash
conda activate liana
code .
```

Run the pipeline in order:

``` bash
python scripts/py_scripts/01_process_h5ad.py
python scripts/py_scripts/02_run_liana.py
python scripts/py_scripts/03_combine_results.py
python scripts/py_scripts/04_make figures.py
```

This workflow:

1.  Prepares the AnnData (.h5ad) object for LIANA+,
2.  infers ligand-receptor interactions independently for each
    developmental or experimental timepoint,
3.  combines results across time points
4.  generates publication-quality heatmaps and summary figures

### Helper Functions

Helper functions used throughout the R workflows are stored in
`scripts/R_scripts/helpers/`

These functions provide reusable code for:

- pseudobulk processing,
- differential expression with edgeR,
- differential expression with Seurat,
- Gene Ontology analysis,
- visualisation,
- file management.

## Outputs

Each workflow writes its outputs to the corresponding subdirectory
within `results/`.

Typical outputs include:

- intermediate objects, (`.rds`)
- differential expression tables (`.csv`),
- Gene Ontology enrichment tables, (`.csv`),
- LIANA+ interaction tables (`.csv`),
- publication-quality figures (`.pdf`)

## Reproducibility

This project was developed in RStudio (v2026.01.1+403 “Apple Blossom”)
using R (v4.6.1). Cell-cell communication analyses were performed using
Python (v3.11) within a dedicated Conda environment containing LIANA+
and its dependencies.

Relative file paths are managed automatically within each script using
`pathlib` (Python) and project-based paths (R), allowing the repository
to be relocated without modifying hard-coded directories.

All intermediate analysis objects are saved as `.rds` files within
`results/objects`, enabling analyses to be resumed without recomputing
earlier steps.
