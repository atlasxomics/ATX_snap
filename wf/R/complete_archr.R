# Shared by the full workflow and recovery of previously returned projects.

validate_arrow_matrices <- function(arrow_files, matrices) {
  for (matrix in matrices) {
    reference <- NULL
    for (arrow in arrow_files) {
      path <- paste0(matrix, "/Info/FeatureDF")
      features <- tryCatch(
        rhdf5::h5read(arrow, path),
        error = function(e) stop("Incomplete ", matrix, " in ", arrow, ": ", e$message)
      )
      if (NROW(features) == 0) stop("Empty ", path, " in ", arrow)
      if (!is.null(reference) && !identical(features, reference)) {
        stop("Inconsistent ", matrix, " feature definitions across Arrow files: ", arrow)
      }
      reference <- features
      if (matrix == "MotifMatrix") {
        completed <- rhdf5::h5read(arrow, "MotifMatrix/Info/Completed")
        if (!identical(as.character(completed), "Finished")) {
          stop("Incomplete motif deviation calculation in ", arrow)
        }
        if (!all(c("z", "deviations") %in% features$seqnames)) {
          stop("MotifMatrix must contain both z scores and deviations in ", arrow)
        }
      }
    }
  }
  invisible(TRUE)
}

rebuild_final_motifs <- function(proj, genome) {
  if (length(proj@peakSet) == 0) stop("Cannot compute motifs without a peak set.")
  proj <- add_motif_annotations(proj, genome)
  proj <- ArchR::addBgdPeaks(proj, force = TRUE)
  proj <- ArchR::addDeviationsMatrix(
    proj, peakAnnotation = "Motif", matrixName = "MotifMatrix",
    out = c("z", "deviations"), force = TRUE
  )
  validate_arrow_matrices(ArchR::getArrowFiles(proj), "MotifMatrix")
  proj
}

validate_complete_archr_project <- function(proj) {
  matrices <- c("GeneScoreMatrix", "PeakMatrix", "MotifMatrix")
  validate_arrow_matrices(ArchR::getArrowFiles(proj), matrices)
  metadata <- ArchR::getCellColData(proj)
  arrows <- ArchR::getArrowFiles(proj)
  for (sample in names(arrows)) {
    expected <- sub("^[^#]+#", "", rownames(metadata)[metadata$Sample == sample])
    for (matrix in matrices) {
      available <- rhdf5::h5read(arrows[[sample]], paste0(matrix, "/Info/CellNames"))
      if (!all(expected %in% available)) {
        stop(matrix, " is missing project cells in sample ", sample)
      }
    }
  }
  if (length(proj@peakSet) == 0) stop("Saved project has no peak set.")
  annotation <- ArchR::getPeakAnnotation(proj, "Motif")
  matches <- readRDS(annotation$Matches)
  same_ranges <- function(x, y) {
    identical(as.character(GenomicRanges::seqnames(x)), as.character(GenomicRanges::seqnames(y))) &&
      identical(BiocGenerics::start(x), BiocGenerics::start(y)) &&
      identical(BiocGenerics::end(x), BiocGenerics::end(y))
  }
  if (!same_ranges(SummarizedExperiment::rowRanges(matches), proj@peakSet)) {
    stop("Motif annotations do not match the active peak set.")
  }
  # Read actual matrix data for one cell per sample, bounding memory while
  # exercising every Arrow file through ArchR's public matrix reader.
  cells <- rownames(metadata)[!duplicated(metadata$Sample)]
  subset <- proj[cells, ]
  for (matrix in matrices) {
    mat <- ArchR::getMatrixFromProject(subset, useMatrix = matrix, threads = 1)
    if (nrow(mat) == 0 || !setequal(colnames(mat), cells)) {
      stop("Cannot read ", matrix, " for all samples in saved project.")
    }
    if (matrix == "PeakMatrix" &&
        !same_ranges(SummarizedExperiment::rowRanges(mat), proj@peakSet)) {
      stop("PeakMatrix does not match the active peak set.")
    }
    # Legacy motif scores use cluster peaks; later sample/condition calls may
    # change the active peak annotations. Validate readability, not equivalence
    # to a newly computed motif matrix on that different peak set.
  }
  invisible(TRUE)
}

save_complete_archr_project <- function(proj, output_dir) {
  # Saving to a different directory drops ArchR's imputation-weight metadata.
  # Preserve the weights and relocate their files alongside the other assets.
  weights <- ArchR::getImputeWeights(proj)
  source_dir <- ArchR::getOutputDirectory(proj)
  ArchR::saveArchRProject(proj, outputDirectory = output_dir, load = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  saved <- ArchR::loadArchRProject(output_dir, force = FALSE, showLogo = FALSE)
  saved <- rebase_group_coverage_paths(saved, output_dir)
  if (length(weights$Weights) > 0) {
    weight_dir <- file.path(output_dir, "ImputeWeights")
    dir.create(weight_dir, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_along(weights$Weights)) {
      source <- weights$Weights[[i]]
      if (!file.exists(source)) source <- file.path(source_dir, "ImputeWeights", basename(source))
      target <- file.path(weight_dir, basename(source))
      if (!file.exists(source)) stop("Missing imputation weight file: ", source)
      if (normalizePath(source) != normalizePath(target, mustWork = FALSE) &&
          !file.copy(source, target, overwrite = TRUE)) stop("Could not copy weight file: ", source)
      weights$Weights[[i]] <- target
    }
    saved@imputeWeights <- weights
  }
  saveRDS(saved, file.path(output_dir, "Save-ArchR-Project.rds"))
  saved <- ArchR::loadArchRProject(output_dir, force = FALSE, showLogo = FALSE)
  validate_complete_archr_project(saved)
  invisible(saved)
}
