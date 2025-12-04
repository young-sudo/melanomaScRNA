#!/usr/bin/env Rscript

# 05-nonmalignant_clustering_and_qc_plots.R
library(Seurat)
library(optparse)

non_malignant <- readRDS("results/non-malignant/dimreduce/non_malignant_dimreduced.rds")
out_dir <- "results/non_malignant/cluster"

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=non_malignant, help="Input file"),
  make_option(c("-o", "--outdir"), type="character", default=out_dir, help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
non_malignant <- opt$input
out_dir <- opt$outdir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# Determine PCs, neighbors and clusters (parameters tunable)
# DimHeatmap(non_malignant, dims = 1:15, cells = 500, balanced = TRUE)
# ElbowPlot(non_malignant)
non_malignant <- FindNeighbors(non_malignant, dims = 1:15)
non_malignant <- FindClusters(non_malignant, resolution = 0.1, random.seed = 40)

# Save cluster annotations
outfile <- file.path(out_dir, "non_malignant_clustered.rds")
saveRDS(non_malignant, file = outfile)

# Create and save some diagnostic plots programmatically (png files)
plots_dir <- file.path(out_dir, "figures")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

png(file.path(plots_dir, "nonmalignant_umap_celltype.png"), width = 800, height = 800)
print(DimPlot(non_malignant, reduction = "umap", group.by = "cell_type", pt.size = 1.2))
dev.off()

png(file.path(plots_dir, "nonmalignant_umap_tumorid.png"), width = 800, height = 800)
print(DimPlot(non_malignant, reduction = "umap", group.by = "tumor_id", pt.size = 1.2))
dev.off()

message("Saved clustering results and UMAP plots in ", out_dir)
