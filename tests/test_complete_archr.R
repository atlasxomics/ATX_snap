# Run from the repository root with Rscript --vanilla tests/test_complete_archr.R
library(rhdf5)
source("wf/R/complete_archr.R")

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

root <- tempfile("matrix-validation-")
dir.create(root)
arrows <- file.path(root, c("one.arrow", "two.arrow"))
matrices <- c("GeneScoreMatrix", "PeakMatrix", "MotifMatrix")
features <- data.frame(seqnames = c("z", "deviations"), idx = c(1L, 1L), name = c("TF1", "TF1"))
for (arrow in arrows) {
  h5createFile(arrow)
  for (matrix in matrices) {
    h5createGroup(arrow, matrix)
    h5createGroup(arrow, paste0(matrix, "/Info"))
    h5write(features, arrow, paste0(matrix, "/Info/FeatureDF"))
  }
  h5write("Finished", arrow, "MotifMatrix/Info/Completed")
}
validate_arrow_matrices(arrows, matrices)

# Reproduce the reported error in the second Arrow: checking just the first
# file or checking the top-level group names would miss this.
h5delete(arrows[2], "MotifMatrix/Info/FeatureDF")
expect_error(validate_arrow_matrices(arrows, matrices), "Incomplete MotifMatrix.*two.arrow")
h5write(features, arrows[2], "MotifMatrix/Info/FeatureDF")

h5write("Started", arrows[2], "MotifMatrix/Info/Completed")
expect_error(validate_arrow_matrices(arrows, matrices), "Incomplete motif deviation")
h5write("Finished", arrows[2], "MotifMatrix/Info/Completed")

h5delete(arrows[2], "MotifMatrix/Info/FeatureDF")
other <- features
other$name <- c("TF2", "TF2")
h5write(other, arrows[2], "MotifMatrix/Info/FeatureDF")
expect_error(validate_arrow_matrices(arrows, matrices), "Inconsistent MotifMatrix")
h5delete(arrows[2], "MotifMatrix/Info/FeatureDF")
h5write(features, arrows[2], "MotifMatrix/Info/FeatureDF")

h5delete(arrows[2], "GeneScoreMatrix/Info/FeatureDF")
expect_error(validate_arrow_matrices(arrows, matrices), "Incomplete GeneScoreMatrix")
h5closeAll()
unlink(root, recursive = TRUE)
cat("HDF5 matrix validation regression checks passed.\n")
