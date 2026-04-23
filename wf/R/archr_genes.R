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

get_gene_feature_names <- function(gene_matrix) {
  feature_names <- gene_matrix@elementMetadata$name
  if (is.null(feature_names)) {
    feature_names <- SummarizedExperiment::rowData(gene_matrix)$name
  }
  feature_names
}

get_impute_weight_files <- function(impute_weights) {
  if (is.null(impute_weights) || length(impute_weights) == 0) {
    return(character(0))
  }

  weight_files <- unlist(impute_weights$Weights, use.names = FALSE)
  if (!is.character(weight_files)) {
    return(character(0))
  }

  weight_files[!is.na(weight_files) & nzchar(weight_files)]
}

repair_impute_weight_paths <- function(proj) {
  impute_weights <- ArchR::getImputeWeights(proj)
  weight_files <- get_impute_weight_files(impute_weights)

  if (length(weight_files) == 0 || all(file.exists(weight_files))) {
    return(proj)
  }

  weights <- impute_weights$Weights
  weight_dir <- file.path(ArchR::getOutputDirectory(proj), "ImputeWeights")
  repaired <- FALSE

  for (i in seq_along(weights)) {
    weight_file <- weights[[i]]
    if (!is.character(weight_file) || file.exists(weight_file)) {
      next
    }

    candidate <- file.path(weight_dir, basename(weight_file))
    if (file.exists(candidate)) {
      message("Repointing impute weight file to: ", candidate)
      weights[[i]] <- candidate
      repaired <- TRUE
    }
  }

  if (repaired) {
    impute_weights$Weights <- weights
    proj@imputeWeights <- impute_weights
  }

  proj
}

select_impute_reduced_dims <- function(proj) {
  reduced_dims <- names(proj@reducedDims)
  preferred <- c("Spectral", "IterativeLSI", "Harmony")

  for (rd_name in preferred) {
    if (rd_name %in% reduced_dims) {
      return(rd_name)
    }
  }

  if (length(reduced_dims) > 0) {
    return(reduced_dims[[1]])
  }

  stop("No reducedDims found in ArchRProject; cannot add impute weights.")
}

ensure_valid_impute_weights <- function(proj) {
  proj <- repair_impute_weight_paths(proj)

  impute_weights <- ArchR::getImputeWeights(proj)
  weight_files <- get_impute_weight_files(impute_weights)
  needs_rebuild <- (
    is.null(impute_weights) ||
      length(impute_weights) == 0 ||
      (length(weight_files) > 0 && any(!file.exists(weight_files)))
  )

  if (!needs_rebuild) {
    return(proj)
  }

  rd_name <- select_impute_reduced_dims(proj)
  message(
    "Impute weights are missing or invalid; rebuilding with reducedDims = ",
    rd_name
  )
  proj <- ArchR::addImputeWeights(proj, reducedDims = rd_name)
  gc(verbose = FALSE, full = TRUE)

  proj
}

create_missing_seurat_object <- function(
  proj,
  run,
  spatial_path,
  chunk_size = 2000
) {
  run_id <- run[1]
  rds_file <- paste0(run_id, "_SeuratObj.rds")

  message("Recreating missing Seurat object for run: ", run_id)

  metadata <- ArchR::getCellColData(ArchRProj = proj)
  rownames(metadata) <- str_split_fixed(
    str_split_fixed(
      row.names(metadata),
      "#",
      2
    )[, 2],
    "-",
    2
  )[, 1]
  metadata["log10_nFrags"] <- log(metadata$nFrags)

  impute_weights <- ArchR::getImputeWeights(proj)
  gene_matrix <- ArchR::getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    asMatrix = TRUE
  )
  gene_row_names <- get_gene_feature_names(gene_matrix)

  gene_assay <- SummarizedExperiment::assay(
    gene_matrix,
    "GeneScoreMatrix"
  )
  run_cols <- grep(pattern = run_id, colnames(gene_assay))
  if (length(run_cols) == 0) {
    stop("No GeneScoreMatrix columns found for missing run: ", run_id)
  }

  n_chunks <- ceiling(nrow(gene_assay) / chunk_size)
  imputed_chunks <- vector("list", n_chunks)

  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, nrow(gene_assay))
    message(
      "Imputing missing Seurat object chunk ",
      i,
      " of ",
      n_chunks,
      " for ",
      run_id
    )

    mat_chunk <- gene_assay[start_idx:end_idx, , drop = FALSE]
    imputed_chunk <- ArchR::imputeMatrix(
      mat = mat_chunk,
      imputeWeights = impute_weights
    )
    imputed_chunks[[i]] <- imputed_chunk[, run_cols, drop = FALSE]
    rm(mat_chunk, imputed_chunk)
    gc(verbose = FALSE, full = TRUE)
  }

  matrix <- do.call(rbind, imputed_chunks)
  rownames(matrix) <- gene_row_names
  rm(imputed_chunks, gene_matrix, gene_assay, impute_weights)
  gc(verbose = FALSE, full = TRUE)

  obj <- build_atlas_seurat_object(
    run_id = run_id,
    matrix = matrix,
    metadata = metadata,
    spatial_path = spatial_path
  )
  saveRDS(obj, file = rds_file)
  rm(matrix, metadata)
  gc(verbose = FALSE, full = TRUE)

  obj
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

  missing_seurat_runs <- c()
  for (run in runs) {
    h5ad_file <- paste0(run[1], "_g_converted.h5ad")
    rds_file <- paste0(run[1], "_SeuratObj.rds")
    if (!file.exists(h5ad_file) && !file.exists(rds_file)) {
      missing_seurat_runs <- c(missing_seurat_runs, run[1])
    }
  }

  if (length(missing_seurat_runs) > 0) {
    message(
      "Missing Seurat object(s) require imputation for run(s): ",
      paste(missing_seurat_runs, collapse = ", ")
    )
    proj <- ensure_valid_impute_weights(proj)
  }

  for (run in runs) {
    h5ad_file <- paste0(run[1], "_g_converted.h5ad")
    if (file.exists(h5ad_file)) {
      message("Found existing ", h5ad_file, "; skipping conversion.")
      next
    }

    rds_file <- paste0(run[1], "_SeuratObj.rds")
    if (!file.exists(rds_file)) {
      obj <- create_missing_seurat_object(
        proj = proj,
        run = run,
        spatial_path = run[5]
      )
    } else {
      message("Converting ", rds_file, " to h5ad...")
      obj <- readRDS(rds_file)
    }

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

# Gene differential statistics -------------------------------------------------
message("Skipping gene differential statistics for this run.")
write.csv(empty_result_df(), "ranked_genes_per_cluster.csv", row.names = FALSE)
write.csv(empty_result_df(), "genes_per_cluster_hm.csv", row.names = FALSE)
empty_pdf("heatmap_genes.pdf", "Gene differential statistics skipped")

if (n_samples > 1) {
  write.csv(empty_result_df(), "ranked_genes_per_sample.csv", row.names = FALSE)
  write.csv(empty_result_df(), "genes_per_sample_hm.csv", row.names = FALSE)

  if (length(sample_name_map) > 0) {
    write.csv(
      empty_result_df(),
      "ranked_genes_per_sample_name.csv",
      row.names = FALSE
    )
    write.csv(empty_result_df(), "genes_per_sample_name_hm.csv", row.names = FALSE)
  }
}

if (n_cond > 1) {
  for (i in seq_along(treatment)) {
    write.csv(
      empty_result_df(),
      paste0("ranked_genes_per_condition_", i, ".csv"),
      row.names = FALSE
    )
    write.csv(
      empty_result_df(),
      paste0("genes_per_condition_", i, "_hm.csv"),
      row.names = FALSE
    )
  }
}

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)
