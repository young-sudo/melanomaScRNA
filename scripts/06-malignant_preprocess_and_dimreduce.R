#!/usr/bin/env Rscript

# 06-malignant_preprocess_and_dimreduce.R
library(Seurat)
library(optparse)

objs <- readRDS("results/malignant/seurat/seurat_split.rds")
out_dir <- "results/malignant/dimreduce"

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=malignant, help="Input file"),
  make_option(c("-o", "--outdir"), type="character", default=out_dir, help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
malignant <- opt$input
out_dir <- opt$outdir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

malignant <- objs$malignant

malignant <- NormalizeData(malignant)
all.genes <- rownames(malignant)
malignant <- ScaleData(malignant, features = all.genes)
malignant <- RunPCA(malignant, features = all.genes)

# tSNE/UMAP
set.seed(40)
malignant <- RunTSNE(malignant, dims = 1:15, seed.use = 40)
set.seed(1)
malignant <- RunUMAP(malignant, dims = 1:15, seed.use = 1)

outfile <- file.path(out_dir, "malignant_dimreduced.rds")
saveRDS(malignant, file = outfile)
message("Saved malignant processed object to ", outfile)


fig_dir <- file.path(out_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# 1. t-SNE colored by tumor ID
tsne_tumor <- DimPlot(
  malignant,
  reduction = "tsne",
  group.by = "tumor_id",
  pt.size = 1.5,
  label = FALSE
) + ggtitle(NULL)

ggsave(
  filename = file.path(fig_dir, "malignant_tsne_tumor_id.png"),
  plot = tsne_tumor,
  width = 6, height = 5, dpi = 300
)

# 2. UMAP colored by tumor ID
umap_tumor <- DimPlot(
  malignant,
  reduction = "umap",
  group.by = "tumor_id",
  pt.size = 1.5,
  label = FALSE
) + ggtitle(NULL)

ggsave(
  filename = file.path(fig_dir, "malignant_umap_tumor_id.png"),
  plot = umap_tumor,
  width = 6, height = 5, dpi = 300
)
