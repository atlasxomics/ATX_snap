#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: rebuild_gene_h5ad.R <Seurat RDS> <output H5AD>")
}

rds_path <- normalizePath(args[[1]], mustWork = TRUE)
output_path <- normalizePath(args[[2]], mustWork = FALSE)
output_dir <- dirname(output_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source("/root/wf/R/seurat.R")

message("Rebuilding gene H5AD from Seurat RDS: ", rds_path)
obj <- readRDS(rds_path)

# Gene RDS objects are intentionally saved before the H5AD-specific cell rename.
# Restore the canonical sample#barcode-1 names used by the Snap observation,
# UMAP, and spatial tables before converting the recovered object.
obj <- rename_cells(list(obj))[[1]]

# The legacy RDS may contain sparse counts. Converting one sample at a time is
# safely below the sparse nonzero limit and avoids repeating ArchR imputation.
counts <- Seurat::GetAssayData(
  object = obj,
  assay = "Spatial",
  slot = "counts"
)
if (inherits(counts, "sparseMatrix")) {
  counts <- as.matrix(counts)
}
obj <- Seurat::SetAssayData(
  object = obj,
  assay = "Spatial",
  slot = "counts",
  new.data = counts
)
rm(counts)
gc(verbose = FALSE, full = TRUE)

temporary_h5seurat <- file.path(output_dir, "rebuild_gene_temp.h5Seurat")
temporary_h5ad <- file.path(output_dir, "rebuild_gene_temp.h5ad")
unlink(c(temporary_h5seurat, temporary_h5ad), force = TRUE)

SeuratDisk::SaveH5Seurat(
  obj,
  filename = temporary_h5seurat,
  overwrite = TRUE
)
SeuratDisk::Convert(
  temporary_h5seurat,
  dest = "h5ad",
  overwrite = TRUE
)

if (!file.exists(temporary_h5ad)) {
  stop("SeuratDisk did not create the expected temporary H5AD.")
}
if (!file.rename(temporary_h5ad, output_path)) {
  stop("Unable to move rebuilt H5AD to: ", output_path)
}
unlink(temporary_h5seurat, force = TRUE)
message("Rebuilt gene H5AD: ", output_path)
