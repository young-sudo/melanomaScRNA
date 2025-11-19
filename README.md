# Melanoma scRNA-seq reanalysis and trajectory modeling

_by Younginn Park_

Reproducible melanoma scRNA-seq pipeline combining Tirosh et al. reanalysis with extended trajectory and pseudotime modeling using Seurat and Monocle3.

The dataset used in this project can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056) with download links ([https](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72056&format=file&file=GSE72056%5Fmelanoma%5Fsingle%5Fcell%5Frevised%5Fv2%2Etxt%2Egz)/[ftp](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/suppl/GSE72056%5Fmelanoma%5Fsingle%5Fcell%5Frevised%5Fv2.txt.gz)).

Reference article: Tirosh et al., 2016 - Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq ([PubMed](https://pubmed.ncbi.nlm.nih.gov/27124452/))

<p align="center">
  <img src="figures/plot_m-tsne-novf-tumor.png" alt="malignant-tsne" width="300"/>
  <img src="figures/plot_nm-tsne-vf-cell.png" alt="malignant-tsne" width="300"/>
</p>

---

# Project Overview

This repository contains a modular, reproducible workflow for analyzing melanoma scRNA-seq data, including both reanalysis of published results and extended trajectory modeling.

### Step-by-step pipeline

1. **Data acquisition and reading** (`00-download_and_read.R`):

   * Download and load the raw Tirosh melanoma single-cell data into R.

2. **Expression and metadata preparation** (`01-prepare_expression_and_metadata.R`):

   * Clean the raw data.
   * Construct gene expression matrix and cell metadata.
   * Handle duplicate genes and remove inconsistent entries.

3. **Seurat object creation and basic QC** (`02-create_seurat_object_and_basic_qc.R`):

   * Create a Seurat object.
   * Add mitochondrial percentage.
   * Initial QC plots and filtering for feature/cell counts.

4. **Malignant and non-malignant cell separation** (`03-filter_and_split_malignant.R`):

   * Subset cells into malignant vs. non-malignant populations.

5. **Non-malignant cells preprocessing and dimensionality reduction** (`04-nonmalignant_preprocess_and_dimreduce.R`):

   * Normalize, scale, identify variable genes.
   * PCA, t-SNE, and UMAP dimensionality reduction.

6. **Non-malignant clustering and visualization** (`05-nonmalignant_clustering_and_qc_plots.R`):

   * Construct SNN graph and find clusters.
   * Visualize clusters and tumor ID distributions.
   * Save diagnostic plots to `figures/`.

7. **Malignant cells preprocessing and dimensionality reduction** (`06-malignant_preprocess_and_dimreduce.R`):

   * Normalize, scale, PCA, t-SNE, UMAP for malignant population.

8. **Trajectory analysis with Monocle3** (`07-trajectory_monocle3.R`):

   * Convert Seurat malignant subset to Monocle3 CellDataSet.
   * Preprocess, dimensionality reduction, cluster, and learn trajectory graph.
   * Optional pseudotime ordering.

---

# Usage

* **Run all scripts sequentially:**

```bash
Rscript run_all.R
```

* **Run individual steps:**

  * Scripts can be run independently in the order: `00` → `01` → `02` → `03` → `04` → `05` → `06` → `07`.

* Ensure that required packages (`Seurat`, `SeuratWrappers`, `monocle3`, `dplyr`, `data.table`, etc.) are installed and available.

