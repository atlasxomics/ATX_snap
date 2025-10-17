library("ArchR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Rnorvegicus.UCSC.rn6")
library("chromVARmotifs")
library("circlize")
library("ComplexHeatmap")
library("data.table")
library("dplyr")
library("GenomicRanges")
library("gridExtra")
library("harmony")
library("plyr")
library("qdap")
library("readr")
library("Seurat")
library("seqLogo")
library("ShinyCell")
library("tidyverse")

source("/root/wf/R/archr.R")
source("/root/wf/R/seurat.R")
source("/root/wf/R/utils.R")
source("/root/wf/R/getDeviation_ArchR.R")

# Globals ---------------------------------------------------------------------
set.seed(42)  # Set seed for reproducibility

args <- commandArgs(trailingOnly = TRUE)
print(args)

project_name <- args[1]
genome <- args[2]
metadata_path <- args[3]
embedding_path <- args[4]
tile_size <- 5000
min_tss <- 0  # Use filtering from SnapATAC2
min_frags <- 0  # Use filtering from SnapATAC2
num_threads <- 50

runs <- strsplit(args[5:length(args)], ",")
runs

inputs <- c()  # Inputs for ArrowFiles (run_id : fragment_file path)
for (run in runs) {
  inputs[run[1]] <- run[2]
}
inputs

out_dir <- paste0(project_name, "_ArchRProject")

# Set genome size for peak calling ----
genome_sizes <- list("hg38" = 3.3e+09, "mm10" = 3.0e+09, "rnor6" = 2.9e+09)
genome_size <- genome_sizes[[genome]]

tempdir <- "/root"
archrproj_dir <- paste0(project_name, "_ArchRProject")

# Create ArchRProject ---------------------------------------------------------
addArchRThreads(threads = num_threads)

proj <- create_archrproject( # from archr.R
  inputs, genome, min_tss, min_frags, tile_size, out_dir
)

# Add Conditions to CellColData ----
for (run in runs) {
  proj$Condition[proj$Sample == run[1]] <- run[3]
}
saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Copy over filtering from SnapATAC2
tryCatch({
  obs <- read.csv(metadata_path)
}, error = function(e) {
  message("Could not read obs.csv: ", e$message)
})

tryCatch({
  cells <- obs[["X"]]
}, error = function(e) {
  cells <- proj$cellName
  message(
    "Expected column 'cluster' not found in obs.csv; skipping filtering...",
    e$message
  )
})

tryCatch({
  clusters <- obs[["cluster"]]
}, error = function(e) {
  clusters <- rep("C0", length(cells))
  message(
    "No cluster information found in obs; assinging cluster 0",
    e$message
  )
})

proj <- proj[proj$cellNames %in% cells]

# Add clusters from SnapATAC2
proj <- addCellColData(
  ArchRProj = proj,
  data = clusters,
  cells = cells,
  name = "Clusters",
  force = TRUE
)

# Parse conditions into 'treatments', add as columns to CellColData ----
conds <- strsplit(proj$Condition, split = "\\s|-")

for (i in seq_along(conds[[1]])) {
  proj <- ArchR::addCellColData(
    proj,
    data = extract_nth_ele(conds, i), # from utils.R
    name = paste0("condition_", i),
    cells = proj$cellNames,
    force = TRUE
  )
}
treatment <- names(ArchR::getCellColData(proj))[
  grep("condition_", names(ArchR::getCellColData(proj)))
]

print(paste("Treatments:", treatment))

# Set global values for project data ----
n_samples <- length(unique(proj$Sample))
n_cond <- length(unique(proj$Condition))
n_cells <- length(proj$cellNames)

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Copy reduced dims from Snap or compute new ----------------------------------

if (file.exists(embedding_path)) {

  message("Using previously computed embeddings from Snap...")
  embedding <- as.matrix(
    read.csv(embedding_path, row.names = 1, header = FALSE)
  )
  archr_cells <- proj@cellColData@rownames
  common <- intersect(archr_cells, rownames(embedding))

  if (length(common) == 0) {
    stop("No overlapping cell names between embedding and ArchR project!")
  }

  # Reorder embedding to match ArchR cell order
  embedding <- embedding[common, , drop = FALSE]
  embedding <- embedding[
    match(archr_cells, rownames(embedding)), , drop = FALSE
  ]

  missing <- sum(is.na(rownames(embedding)))
  if (missing > 0) {
    warning(missing, " ArchR cells have no embedding; filled with NA rows.")
  }

  rd <- SimpleList(
    matDR = embedding, date = Sys.time(), scaleDims = NA, corToDepth = NA
  )
  proj@reducedDims[["Spectral"]] <- rd
  rd_name <- "Spectral"

} else {
  message("No precomputed embedding found. Running IterativeLSI...")
  proj <- add_lsi(proj, 2, 0.7, 25000)
  rd_name <- "IterativeLSI"
}

# Impute weights --------------------------------------------------------------
proj <- ArchR::addImputeWeights(proj, reducedDims = rd_name)

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Gene Expression Analysis ----------------------------------------------------

# Create SeuratObjects for gene matrix ----
# Extract metadata for Seurat object --
metadata <- getCellColData(ArchRProj = proj)

# Set metadata rownames to barcodes
rownames(metadata) <- str_split_fixed(
  str_split_fixed(
    row.names(metadata),
    "#",
    2
  )[, 2],
  "-",
  2
)[, 1]

# Create col for log10 of fragment counts
metadata["log10_nFrags"] <- log(metadata$nFrags)

# Extract gene matrix for SeuratObject --
# Get imputation weights once
impute_weights <- getImputeWeights(proj)

gene_matrix <- ArchR::getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  asMatrix = TRUE
)
gene_row_names <- gene_matrix@elementMetadata$name

# Identify empty features for filtering volcano plots --
print("Identifying empty features...")
empty_feat_idx <- which(Matrix::rowSums(
  SummarizedExperiment::assay(gene_matrix, "GeneScoreMatrix")
) == 0)
empty_feat <- gene_row_names[empty_feat_idx]
print(paste("Found", length(empty_feat), "empty features"))


matrix <- ArchR::imputeMatrix(
  mat = SummarizedExperiment::assay(gene_matrix),
  imputeWeights = impute_weights
)

rm(gene_matrix, impute_weights)
gc()

rownames(matrix) <- gene_row_names

# Create and save SeuratObjects --
print("Creating SeuratObjects...")

seurat_objs <- c()
for (run in runs) {

  obj <- build_atlas_seurat_object( # from seurat.R
    run_id = run[1],
    matrix = matrix,
    metadata = metadata,
    spatial_path = run[5]
  )

  saveRDS(obj, file = paste0(run[1], "_SeuratObj.rds"))
  seurat_objs <- c(seurat_objs, obj)
}

print("Available SeuratObjects:")
seurat_objs

all <- rename_cells(seurat_objs)  # from seurat.R
rm(seurat_objs)
gc()

# Convert Seurat to h5ad and save ----
for (obj in all) {
  seurat_to_h5ad(obj, FALSE, paste0(unique(obj$Sample), "_g"))  # from utils.R
}
rm(all)
gc()

# Identify marker genes ----
# Marker genes per cluster, save marker gene rds, csv, heatmap.csv --
n_clust <- length(unique(proj$Clusters))

cluster_marker_genes <- get_marker_genes( # from archr.R
  proj,
  group_by = "Clusters",
  markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
  heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10"
)

write.csv(
  cluster_marker_genes$marker_list,
  "ranked_genes_per_cluster.csv",
  row.names = FALSE
)
write.csv(cluster_marker_genes$heatmap_gs, "genes_per_cluster_hm.csv")

# Marker genes per sample, save marker gene rds, csv, heatmap.csv --
if (n_samples > 1) {

  sample_marker_genes <- get_marker_genes(
    proj,
    group_by = "Sample",
    markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
    heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10"
  )

  write.csv(
    sample_marker_genes$marker_list,
    "ranked_genes_per_sample.csv",
    row.names = FALSE
  )
  write.csv(sample_marker_genes$heatmap_gs, "genes_per_sample_hm.csv")

}

# Marker genes per treatment, save marker gene rds, csv, heatmap.csv --
if (n_cond > 1) {

  for (i in seq_along(treatment)) {

    treatment_marker_genes <- get_marker_genes(
      proj,
      group_by = treatment[i],
      markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
      heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10"
    )

    write.csv(
      treatment_marker_genes$marker_list,
      paste0("ranked_genes_per_condition_", i, ".csv"),
      row.names = FALSE
    )
    write.csv(
      treatment_marker_genes$heatmap_gs,
      paste0("genes_per_condition_", i, "_hm.csv")
    )
  }
}

# Volcano plots for genes ----
if (n_cond > 1) {
  for (j in seq_along(treatment)) {

    # Get gene markers df for all clusters together --
    marker_genes_df <- get_marker_df(
      proj = proj,
      group_by = treatment[j],
      matrix = "GeneScoreMatrix",
      seq_names = NULL,
      max_cells = n_cells,  # Equals total cells in project
      test_method = "ttest",
      diff_metric = "Log2FC"
    )

    # Create a merged marker genes df for clusters for which no condition is
    # >90% of all cells --
    req_clusters <- get_required_clusters(proj, treatment[j])
    marker_genes_by_cluster_df <- get_marker_df_clusters(
      proj = proj,
      clusters = req_clusters,
      group_by = treatment[j],
      seq_names = "z",
      matrix = "GeneScoreMatrix",
      test_method = "ttest",
      diff_metric = "Log2FC"
    )

    # Per condition, merge dfs and cleanup data --
    conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
    for (cond in conditions) {

      volcano_table <- get_volcano_table( # from archr.R
        marker_genes_df,
        marker_genes_by_cluster_df,
        cond,
        "gene",
        empty_feat,
        fc_col = "Log2FC"
      )

      write.table(
        volcano_table,
        paste0(
          "volcanoMarkers_genes_", j, "_", cond, ".csv"
        ),
        sep = ",",
        quote = FALSE,
        row.names = FALSE
      )
      print(paste0("volcanoMarkers_genes_", j, "_", cond, ".csv is done!"))

      features <- unique(volcano_table$cluster)
      others <- paste(conditions[conditions != cond], collapse = "|")
      volcano_plots <- list()
      for (i in seq_along(features)) {
        volcano_plots[[i]] <- scvolcano(
          volcano_table, cond, others, features[[i]], fc_col = "Log2FC"
        )
      }

      pdf(paste0("volcano_plots_", cond, ".pdf"))
      for (plot in volcano_plots) {
        print(plot)
      }
      dev.off()
    }
  }
} else {
  print("There are not enough conditions to be compared with!")
}

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)
