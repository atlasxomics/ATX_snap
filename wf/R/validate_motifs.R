validate_motif_arrow_files <- function(arrows) {
  reference <- NULL
  for (arrow in arrows) {
    features <- tryCatch(
      rhdf5::h5read(arrow, "MotifMatrix/Info/FeatureDF"),
      error = function(e) stop("Incomplete MotifMatrix in ", arrow, ": ", e$message)
    )
    if (NROW(features) == 0 ||
        !all(c("z", "deviations") %in% features$seqnames)) {
      stop("Missing motif z scores or deviations in ", arrow)
    }
    if (!identical(as.character(rhdf5::h5read(arrow, "MotifMatrix/Info/Completed")), "Finished")) {
      stop("Unfinished MotifMatrix in ", arrow)
    }
    if (!is.null(reference) && !identical(reference, features)) {
      stop("Motif features differ across Arrow files: ", arrow)
    }
    reference <- features
  }
  invisible(TRUE)
}

read_recovery_matrices <- function(proj, matrices) {
  metadata <- ArchR::getCellColData(proj)
  cells <- rownames(metadata)[!duplicated(metadata$Sample)]
  subset <- proj[cells, ]
  result <- setNames(vector("list", length(matrices)), matrices)
  for (matrix in matrices) {
    mat <- ArchR::getMatrixFromProject(subset, useMatrix = matrix, threads = 1)
    if (nrow(mat) == 0 || !setequal(colnames(mat), cells)) {
      stop("Cannot read ", matrix, " for every sample.")
    }
    result[[matrix]] <- mat[, cells, drop = FALSE]
  }
  result
}
