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
gene_artifacts_dir <- if (length(args) >= 5) args[5] else ""
resume_from_gene_artifacts <- (
  !is.na(gene_artifacts_dir) &&
    nzchar(gene_artifacts_dir) &&
    tolower(gene_artifacts_dir) != "none"
)
tile_size <- 5000
min_tss <- 0  # Use filtering from SnapATAC2
min_frags <- 0  # Use filtering from SnapATAC2
num_threads <- 50

run_args <- if (resume_from_gene_artifacts) {
  if (length(args) < 6) character(0) else args[6:length(args)]
} else {
  if (length(args) < 5) character(0) else args[5:length(args)]
}

if (length(run_args) == 0) {
  stop("No runs were supplied to archr_genes.R.")
}

runs <- strsplit(run_args, ",")
runs

inputs <- c()  # Inputs for ArrowFiles (run_id : fragment_file path)
for (run in runs) {
  inputs[run[1]] <- run[2]
}
inputs

# Map run_id -> optional sample_name (if provided)
sample_name_map <- c()
for (run in runs) {
  run_id <- run[1]
  sample_name <- if (length(run) >= 6) trimws(run[6]) else ""
  if (!is.na(sample_name) && nzchar(sample_name) &&
      tolower(sample_name) != "none") {
    sample_name_map[run_id] <- sample_name
  }
}

remap_sample_ids <- function(values, sample_name_map) {
  remapped <- values
  for (i in seq_along(values)) {
    value <- as.character(values[[i]])
    if (value %in% names(sample_name_map)) {
      remapped[[i]] <- sample_name_map[[value]]
    }
  }
  return(remapped)
}

copy_gene_artifact_files <- function(artifact_dir, pattern) {
  files <- list.files(
    artifact_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(files) == 0) {
    return(character(0))
  }

  copied <- character(0)
  for (file in files) {
    destination <- basename(file)
    if (!file.exists(destination)) {
      file.copy(file, destination, overwrite = FALSE)
    }
    copied <- c(copied, destination)
  }

  return(unique(copied))
}

find_archr_project_dir <- function(artifact_dir, project_dir_name) {
  candidates <- c(
    file.path(artifact_dir, project_dir_name),
    artifact_dir
  )

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "Save-ArchR-Project.rds"))) {
      return(candidate)
    }
  }

  recursive_matches <- list.files(
    artifact_dir,
    pattern = "^Save-ArchR-Project\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(recursive_matches) > 0) {
    return(dirname(recursive_matches[[1]]))
  }

  stop(
    "Could not find Save-ArchR-Project.rds in gene artifacts directory: ",
    artifact_dir
  )
}

prepare_resume_seurat_object <- function(obj) {
  if (any(grepl("#", colnames(obj), fixed = TRUE))) {
    return(obj)
  }

  rename_cells(list(obj))[[1]]
}

# Set genome size for peak calling ----
genome_sizes <- list("hg38" = 3.3e+09, "mm10" = 3.0e+09, "rnor6" = 2.9e+09)
genome_size <- genome_sizes[[genome]]

archrproj_dir <- paste0(project_name, "_ArchRProject")

# Create ArchRProject ---------------------------------------------------------
addArchRThreads(threads = num_threads)

if (resume_from_gene_artifacts) {

  message("Resuming gene task from artifacts directory: ", gene_artifacts_dir)
  gene_artifacts_dir <- normalizePath(gene_artifacts_dir, mustWork = TRUE)

  archrproj_source_dir <- find_archr_project_dir(
    gene_artifacts_dir,
    archrproj_dir
  )
  message("Loading ArchRProject from: ", archrproj_source_dir)
  proj <- ArchR::loadArchRProject(
    path = archrproj_source_dir,
    force = TRUE,
    showLogo = FALSE
  )

  # Refresh run metadata in case the resumed project came from a prior launch.
  for (run in runs) {
    proj$Condition[proj$Sample == run[1]] <- run[3]
    proj$sample_name[proj$Sample == run[1]] <- run[6]
  }

  conds <- strsplit(proj$Condition, split = "\\s|-")
  for (i in seq_along(conds[[1]])) {
    proj <- ArchR::addCellColData(
      proj,
      data = extract_nth_ele(conds, i),
      name = paste0("condition_", i),
      cells = proj$cellNames,
      force = TRUE
    )
  }
  treatment <- names(ArchR::getCellColData(proj))[
    grep("condition_", names(ArchR::getCellColData(proj)))
  ]

  print(paste("Treatments:", treatment))

  n_samples <- length(unique(proj$Sample))
  n_cond <- length(unique(proj$Condition))
  n_cells <- length(proj$cellNames)
  empty_feat <- character(0)

  copied_h5ads <- copy_gene_artifact_files(
    gene_artifacts_dir,
    "_g_converted\\.h5ad$"
  )
  copied_rds <- copy_gene_artifact_files(
    gene_artifacts_dir,
    "_SeuratObj\\.rds$"
  )

  message("Copied ", length(copied_h5ads), " existing gene h5ad file(s).")
  message("Copied ", length(copied_rds), " Seurat object RDS file(s).")

  for (run in runs) {
    h5ad_file <- paste0(run[1], "_g_converted.h5ad")
    if (file.exists(h5ad_file)) {
      message("Found existing ", h5ad_file, "; skipping conversion.")
      next
    }

    rds_file <- paste0(run[1], "_SeuratObj.rds")
    if (!file.exists(rds_file)) {
      stop(
        "Missing required Seurat object for run '",
        run[1],
        "': ",
        rds_file
      )
    }

    message("Converting ", rds_file, " to h5ad...")
    obj <- readRDS(rds_file)
    obj <- prepare_resume_seurat_object(obj)
    seurat_to_h5ad(obj, FALSE, paste0(unique(obj$Sample), "_g"))
    rm(obj)
    gc(verbose = FALSE, full = TRUE)
  }

} else {

proj <- create_archrproject( # from archr.R
  inputs, genome, min_tss, min_frags, tile_size, archrproj_dir
)

# Add Conditions to CellColData ----
for (run in runs) {
  proj$Condition[proj$Sample == run[1]] <- run[3]
  proj$sample_name[proj$Sample == run[1]] <- run[6]
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
  clusters <- as.character(obs[["cluster"]])
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

# Chunked imputation parameters
chunk_size <- 2000
total_features <- nrow(SummarizedExperiment::assay(
  gene_matrix,
  "GeneScoreMatrix"
))
n_chunks <- ceiling(total_features / chunk_size)

print(paste(Sys.time(), "Starting chunked imputation..."))
print(paste(
  "Processing", total_features, "features in", n_chunks,
  "chunks of size", chunk_size
))

imputed_chunks <- vector("list", n_chunks)

for (i in seq_len(n_chunks)) {

  start_idx <- (i - 1) * chunk_size + 1
  end_idx <- min(i * chunk_size, total_features)

  print(paste(Sys.time(), "Processing chunk", i, "of", n_chunks))

  mat_chunk <- SummarizedExperiment::assay(
    gene_matrix,
    "GeneScoreMatrix"
  )[start_idx:end_idx, , drop = FALSE]

  imputed_chunk <- ArchR::imputeMatrix(
    mat = mat_chunk,
    imputeWeights = impute_weights
  )

  imputed_chunks[[i]] <- imputed_chunk
  rm(mat_chunk, imputed_chunk)

  print(paste(
    "Completed", i, "Mem:", format(object.size(imputed_chunks), units = "MB")
  ))
}

print(paste(Sys.time(), "Combining chunks..."))
matrix <- do.call(rbind, imputed_chunks)

rm(imputed_chunks)
gc(verbose = FALSE)
print(paste(Sys.time(), "Imputation complete!"))

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
rm(all, matrix, metadata, gene_row_names)
gc(verbose = FALSE, full = TRUE)

}

# Identify marker genes ----
# Marker genes per cluster, save marker gene rds, csv, heatmap.csv --
n_clust <- length(unique(proj$Clusters))
marker_max_cells <- 500
marker_threads <- 1
print(paste(
  "Using maxCells =", marker_max_cells,
  "and threads =", marker_threads,
  "for ArchR gene marker tests to limit memory usage."
))

cluster_marker_genes <- get_marker_genes( # from archr.R
  proj,
  group_by = "Clusters",
  markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
  heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
  max_cells = marker_max_cells,
  threads = marker_threads
)

write.csv(
  cluster_marker_genes$marker_list,
  "ranked_genes_per_cluster.csv",
  row.names = FALSE
)
write.csv(cluster_marker_genes$heatmap_gs, "genes_per_cluster_hm.csv")

# Recompute gene heatmap for plotting (transpose, plotLog2FC different) --
cut_off <- "Pval <= 0.05 & Log2FC >= 0.1"

if (
  is.null(cluster_marker_genes$markers_gs) ||
    nrow(cluster_marker_genes$heatmap_gs) == 0 ||
    ncol(cluster_marker_genes$heatmap_gs) == 0
) {
  message("Skipping gene heatmap plot because no cluster marker genes were available.")
  empty_pdf("heatmap_genes.pdf", "No cluster marker genes available")
} else {
  heatmap_gs_plotting <- plotMarkerHeatmap(
    seMarker = cluster_marker_genes$markers_gs,
    cutOff = cut_off,
    transpose = TRUE
  )

  gene_hm <- ComplexHeatmap::draw(
    heatmap_gs_plotting,
    heatmap_legend_side = "bot",
    annotation_legend_side = "bot",
    column_title = paste0("Marker genes (", cut_off, ")"),
    column_title_gp = gpar(fontsize = 12)
  )

  print("Saving gene heatmap...")
  pdf("heatmap_genes.pdf")
  print(gene_hm)
  dev.off()
}
rm(list = intersect(c("cluster_marker_genes", "heatmap_gs_plotting", "gene_hm"), ls()))
gc(verbose = FALSE, full = TRUE)

# Marker genes per sample, save marker gene rds, csv, heatmap.csv --
if (n_samples > 1) {

  sample_marker_genes <- get_marker_genes(
    proj,
    group_by = "Sample",
    markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
    heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
    max_cells = marker_max_cells,
    threads = marker_threads
  )

  write.csv(
    sample_marker_genes$marker_list,
    "ranked_genes_per_sample.csv",
    row.names = FALSE
  )
  write.csv(sample_marker_genes$heatmap_gs, "genes_per_sample_hm.csv")

  if (length(sample_name_map) > 0) {
    sample_marker_list_named <- sample_marker_genes$marker_list
    marker_groups <- names(sample_marker_list_named)
    if (!is.data.frame(sample_marker_list_named) && !is.null(marker_groups)) {
      names(sample_marker_list_named) <- remap_sample_ids(
        marker_groups, sample_name_map
      )
    }
    write.csv(
      sample_marker_list_named,
      "ranked_genes_per_sample_name.csv",
      row.names = FALSE
    )

    sample_heatmap_named <- sample_marker_genes$heatmap_gs
    if (!is.null(colnames(sample_heatmap_named))) {
      colnames(sample_heatmap_named) <- remap_sample_ids(
        colnames(sample_heatmap_named), sample_name_map
      )
    }
    if (!is.null(rownames(sample_heatmap_named))) {
      rownames(sample_heatmap_named) <- remap_sample_ids(
        rownames(sample_heatmap_named), sample_name_map
      )
    }
    write.csv(
      sample_heatmap_named,
      "genes_per_sample_name_hm.csv"
    )
  }

  rm(list = intersect(
    c("sample_marker_genes", "sample_marker_list_named", "sample_heatmap_named"),
    ls()
  ))
  gc(verbose = FALSE, full = TRUE)
}

# Marker genes per treatment, save marker gene rds, csv, heatmap.csv --
if (n_cond > 1) {

  for (i in seq_along(treatment)) {

    treatment_marker_genes <- get_marker_genes(
      proj,
      group_by = treatment[i],
      markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
      heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
      max_cells = marker_max_cells,
      threads = marker_threads
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
    rm(treatment_marker_genes)
    gc(verbose = FALSE, full = TRUE)
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
      max_cells = marker_max_cells,
      test_method = "ttest",
      diff_metric = "Log2FC",
      threads = marker_threads
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
      max_cells = marker_max_cells,
      test_method = "ttest",
      diff_metric = "Log2FC",
      threads = marker_threads
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
      rm(volcano_table, volcano_plots)
      gc(verbose = FALSE, full = TRUE)
    }
    rm(marker_genes_df, marker_genes_by_cluster_df)
    gc(verbose = FALSE, full = TRUE)
  }
} else {
  print("There are not enough conditions to be compared with!")
}

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)
