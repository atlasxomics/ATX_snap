library("ArchR")
source("/root/wf/R/load_genomes.R")
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
library("TxDb.Mmusculus.UCSC.mm39.knownGene")
library("org.Mm.eg.db")

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

is_bool_arg <- function(value) {
  if (length(value) == 0 || is.na(value)) {
    return(FALSE)
  }
  tolower(trimws(as.character(value))) %in%
    c("true", "t", "1", "yes", "y", "false", "f", "0", "no", "n")
}

remaining_args <- function(start_idx) {
  if (length(args) < start_idx) {
    return(character(0))
  }
  args[start_idx:length(args)]
}

next_arg_idx <- 5
include_y_chromosome <- FALSE
if (length(args) >= next_arg_idx && is_bool_arg(args[next_arg_idx])) {
  include_y_chromosome <- parse_bool_arg(args[next_arg_idx])
  next_arg_idx <- next_arg_idx + 1
}

stage_mode <- "all"
if (
  length(args) >= next_arg_idx &&
    tolower(args[next_arg_idx]) %in% c("all", "prepare", "export")
) {
  stage_mode <- tolower(args[next_arg_idx])
  next_arg_idx <- next_arg_idx + 1
}

gene_artifacts_dir <- if (length(args) >= next_arg_idx) args[next_arg_idx] else ""
resume_from_gene_artifacts <- (
  !is.na(gene_artifacts_dir) &&
    nzchar(gene_artifacts_dir) &&
    tolower(gene_artifacts_dir) != "none" &&
    !grepl(",", gene_artifacts_dir, fixed = TRUE) &&
    dir.exists(gene_artifacts_dir)
)
tile_size <- 5000
min_tss <- 0  # Use filtering from SnapATAC2
min_frags <- 0  # Use filtering from SnapATAC2
num_threads <- 24

run_args <- if (resume_from_gene_artifacts) {
  remaining_args(next_arg_idx + 1)
} else {
  remaining_args(next_arg_idx)
}

if (length(run_args) == 0) {
  stop("No runs were supplied to archr_genes.R.")
}

if (stage_mode == "export" && !resume_from_gene_artifacts) {
  stop("Export mode requires a valid checkpointed gene artifacts directory.")
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

add_project_impute_weights <- function(proj, rd_name) {
  n_cells <- length(proj$cellNames)
  if (n_cells < 2) {
    stop(
      "Cannot add ArchR impute weights with fewer than 2 cells after filtering."
    )
  }

  impute_k <- min(15, n_cells - 1)
  if (impute_k < 15) {
    message(
      "Using addImputeWeights(k = ",
      impute_k,
      ") because only ",
      n_cells,
      " cells are available after filtering."
    )
  }

  ArchR::addImputeWeights(
    proj,
    reducedDims = rd_name,
    k = impute_k
  )
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
  proj <- add_project_impute_weights(proj, rd_name)
  gc(verbose = FALSE, full = TRUE)

  proj
}

write_gene_score_chunks <- function(
  proj,
  runs,
  chunk_root,
  feature_chunk_size = 500
) {
  dir.create(chunk_root, recursive = TRUE, showWarnings = FALSE)

  run_ids <- vapply(runs, function(run) run[[1]], character(1))
  run_chunk_dirs <- file.path(
    chunk_root,
    sprintf("run_%03d", seq_along(run_ids))
  )
  names(run_chunk_dirs) <- run_ids
  for (run_chunk_dir in run_chunk_dirs) {
    dir.create(run_chunk_dir, recursive = TRUE, showWarnings = FALSE)
  }

  run_cells <- lapply(run_ids, function(run_id) {
    proj$cellNames[as.character(proj$Sample) == run_id]
  })
  names(run_cells) <- run_ids

  empty_runs <- names(run_cells)[lengths(run_cells) == 0]
  if (length(empty_runs) > 0) {
    stop(
      "No ArchR cells found for run(s): ",
      paste(empty_runs, collapse = ", ")
    )
  }

  impute_weights <- ArchR::getImputeWeights(proj)
  if (is.null(impute_weights) || length(impute_weights) == 0) {
    stop("No imputation weights found in the checkpointed ArchR project.")
  }

  matrix_seqnames <- ArchR::getSeqnames(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix"
  )
  if (length(matrix_seqnames) == 0) {
    stop("No GeneScoreMatrix seqnames found in the ArchR project.")
  }

  chunk_index <- 0
  for (seq_index in seq_along(matrix_seqnames)) {
    seqname <- matrix_seqnames[[seq_index]]
    message(
      "Reading GeneScoreMatrix seqname ",
      seq_index,
      " of ",
      length(matrix_seqnames),
      ": ",
      seqname
    )

    # The pinned ArchR fork densifies each Arrow matrix when asMatrix is TRUE.
    # Restricting each call to one seqname bounds that allocation instead of
    # materializing every gene across all cells at once.
    gene_matrix <- ArchR::getMatrixFromProject(
      ArchRProj = proj,
      useMatrix = "GeneScoreMatrix",
      useSeqnames = seqname,
      threads = num_threads,
      asMatrix = TRUE
    )
    gene_assay <- SummarizedExperiment::assay(
      gene_matrix,
      "GeneScoreMatrix"
    )
    feature_names <- get_gene_feature_names(gene_matrix)

    if (length(feature_names) != nrow(gene_assay)) {
      stop("Gene feature names do not match matrix rows for seqname ", seqname)
    }

    run_col_indices <- lapply(run_cells, function(cells) {
      match(cells, colnames(gene_assay))
    })
    missing_runs <- names(run_col_indices)[vapply(
      run_col_indices,
      anyNA,
      logical(1)
    )]
    if (length(missing_runs) > 0) {
      stop(
        "GeneScoreMatrix is missing cells for run(s): ",
        paste(missing_runs, collapse = ", ")
      )
    }

    n_feature_chunks <- ceiling(nrow(gene_assay) / feature_chunk_size)
    for (feature_index in seq_len(n_feature_chunks)) {
      chunk_index <- chunk_index + 1
      start_idx <- (feature_index - 1) * feature_chunk_size + 1
      end_idx <- min(feature_index * feature_chunk_size, nrow(gene_assay))
      feature_idx <- start_idx:end_idx

      message(
        "Imputing feature chunk ",
        feature_index,
        " of ",
        n_feature_chunks,
        " for ",
        seqname,
        " (global chunk ",
        chunk_index,
        ")"
      )

      mat_chunk <- gene_assay[feature_idx, , drop = FALSE]
      imputed_chunk <- ArchR::imputeMatrix(
        mat = mat_chunk,
        imputeWeights = impute_weights,
        threads = num_threads
      )

      # The pinned ArchR imputeMatrix implementation subsets its final result
      # without drop = FALSE. A one-feature chunk is therefore returned as a
      # dimensionless vector. Restore the feature-by-cell shape explicitly and
      # validate all other chunks before persisting them.
      expected_dim <- dim(mat_chunk)
      if (is.null(dim(imputed_chunk))) {
        if (length(imputed_chunk) != prod(expected_dim)) {
          stop(
            "Dimensionless imputed chunk has ",
            length(imputed_chunk),
            " values; expected ",
            prod(expected_dim)
          )
        }
        imputed_chunk <- matrix(
          imputed_chunk,
          nrow = expected_dim[[1]],
          ncol = expected_dim[[2]]
        )
      }
      if (!identical(as.integer(dim(imputed_chunk)), as.integer(expected_dim))) {
        stop(
          "Imputed chunk dimensions ",
          paste(dim(imputed_chunk), collapse = " x "),
          " do not match expected dimensions ",
          paste(expected_dim, collapse = " x ")
        )
      }
      rownames(imputed_chunk) <- feature_names[feature_idx]
      colnames(imputed_chunk) <- colnames(mat_chunk)

      for (run_id in run_ids) {
        # Imputed gene scores are effectively dense. Keep each bounded,
        # per-sample chunk dense instead of paying the larger dgCMatrix
        # overhead or introducing its nonzero-element limit.
        run_chunk <- imputed_chunk[
          , run_col_indices[[run_id]], drop = FALSE
        ]
        if (!is.matrix(run_chunk) || inherits(run_chunk, "sparseMatrix")) {
          run_chunk <- as.matrix(run_chunk)
        }
        chunk_path <- file.path(
          run_chunk_dirs[[run_id]],
          sprintf("chunk_%05d.rds", chunk_index)
        )
        saveRDS(run_chunk, file = chunk_path, compress = FALSE)
        rm(run_chunk)
      }

      rm(mat_chunk, imputed_chunk)
      gc(verbose = FALSE, full = TRUE)
    }

    rm(gene_matrix, gene_assay, feature_names, run_col_indices)
    gc(verbose = FALSE, full = TRUE)
  }

  rm(impute_weights)
  gc(verbose = FALSE, full = TRUE)
  run_chunk_dirs
}

export_gene_score_objects <- function(runs, metadata, run_chunk_dirs) {
  for (run_index in seq_along(runs)) {
    run <- runs[[run_index]]
    run_id <- run[[1]]
    chunk_files <- sort(list.files(
      run_chunk_dirs[[run_id]],
      pattern = "^chunk_[0-9]+\\.rds$",
      full.names = TRUE
    ))

    if (length(chunk_files) == 0) {
      stop("No persisted gene chunks found for run: ", run_id)
    }

    message(
      "Assembling ",
      length(chunk_files),
      " persisted gene chunks for run ",
      run_id
    )
    run_chunks <- lapply(chunk_files, readRDS)
    matrix <- do.call(rbind, run_chunks)
    if (!is.matrix(matrix) || inherits(matrix, "sparseMatrix")) {
      matrix <- as.matrix(matrix)
    }
    rm(run_chunks)
    gc(verbose = FALSE, full = TRUE)

    obj <- build_atlas_seurat_object(
      run_id = run_id,
      matrix = matrix,
      metadata = metadata,
      spatial_path = run[[5]],
      # Seurat v5 stores assay counts sparsely. The task converts and verifies
      # the exported per-run H5AD as dense immediately after this R script.
      sparse_counts = TRUE
    )
    saveRDS(obj, file = paste0(run_id, "_SeuratObj.rds"))

    h5ad_obj <- rename_cells(list(obj))[[1]]
    seurat_to_h5ad(
      h5ad_obj,
      FALSE,
      paste0(unique(h5ad_obj$Sample), "_g")
    )

    rm(matrix, obj, h5ad_obj)
    unlink(run_chunk_dirs[[run_id]], recursive = TRUE, force = TRUE)
    gc(verbose = FALSE, full = TRUE)
  }
}

# Set genome size for peak calling ----
genome_sizes <- list(
  "hg38" = 3.3e+09,
  "mm10" = 3.0e+09,
  "mm39" = 2.7e+09,
  "rnor6" = 2.9e+09
)
genome_size <- genome_sizes[[genome]]
if (is.null(genome_size)) {
  stop("No genome size configured for genome: ", genome)
}

archrproj_dir <- paste0(project_name, "_ArchRProject")

# Create ArchRProject ---------------------------------------------------------
addArchRThreads(threads = num_threads)

if (resume_from_gene_artifacts) {

  message("Loading checkpointed gene project from: ", gene_artifacts_dir)
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
  proj <- ensure_valid_impute_weights(proj)

} else {

  proj <- create_archrproject( # from archr.R
    inputs,
    genome,
    min_tss,
    min_frags,
    tile_size,
    archrproj_dir,
    include_y_chromosome = include_y_chromosome
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
      "No cluster information found in obs; assigning cluster 0",
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

  # Copy reduced dims from Snap or compute new -------------------------------
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

  # Impute weights ------------------------------------------------------------
  proj <- add_project_impute_weights(proj, rd_name)
  saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)
}

if (stage_mode == "prepare") {
  message(
    "ArchR gene-project checkpoint complete; skipping GeneScoreMatrix export."
  )
  quit(save = "no", status = 0, runLast = FALSE)
}

# Gene Expression Analysis ----------------------------------------------------
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

chunk_root <- tempfile(pattern = "gene_score_chunks_", tmpdir = "/root")
run_chunk_dirs <- write_gene_score_chunks(
  proj = proj,
  runs = runs,
  chunk_root = chunk_root,
  feature_chunk_size = 500
)
export_gene_score_objects(
  runs = runs,
  metadata = metadata,
  run_chunk_dirs = run_chunk_dirs
)
unlink(chunk_root, recursive = TRUE, force = TRUE)
rm(metadata, run_chunk_dirs)
gc(verbose = FALSE, full = TRUE)

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
