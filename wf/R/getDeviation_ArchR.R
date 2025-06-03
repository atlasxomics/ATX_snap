library(ArchR)
library(Matrix)

getDeviation_ArchR <- function(
  ArchRProj,
  name,
  imputeWeights = getImputeWeights(ArchRProj),
  useAssay = "z",
  log2Norm = TRUE
) {
  motif_se <- getMatrixFromProject(ArchRProj, useMatrix = "MotifMatrix")
  mat <- assay(motif_se, useAssay)

  motif_names <- rowData(motif_se)$name
  idx <- which(tolower(motif_names) %in% tolower(name))
  if (length(idx) == 0) {
    stop("None of the provided motif names were found in the matrix.")
  }

  if (!is.null(imputeWeights)) {
    message("Applying imputation...")
    mat <- imputeMatrix(mat = as.matrix(mat), imputeWeights = imputeWeights)
  }

  result <- mat[idx, , drop = FALSE]
  rownames(result) <- motif_names[idx]

  if (log2Norm) {
    result <- log2(result + 1)
  }

  colnames(result) <- gsub(".*#", "", colnames(result))
  colnames(result) <- gsub("-.*", "", colnames(result))

  return(as.data.frame(as.matrix(result)))
}