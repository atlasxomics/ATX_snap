library("ArchR")
source("/root/wf/R/load_genomes.R")
source("/root/wf/R/archr.R")
source("/root/wf/R/validate_motifs.R")

set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: restore_motifs.R <copied_project> <genome> <threads>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
genome <- args[[2]]
if (!genome %in% c("hg38", "mm10", "mm39", "rnor6")) stop("Unsupported genome: ", genome)
ArchR::addArchRThreads(threads = as.integer(args[[3]]))
proj <- ArchR::loadArchRProject(project_dir, force = TRUE, showLogo = FALSE)
proj <- rebase_group_coverage_paths(proj, project_dir)
if (length(proj@peakSet) == 0) {
  stop("Input project has no peak set. Supply the project with restored peaks.")
}
original_peaks <- proj@peakSet
original_cells <- ArchR::getCellColData(proj)
original_dims <- proj@reducedDims
original_embeddings <- proj@embeddings
# Read real data from every sample before altering the copied Arrow files.
original_matrices <- read_recovery_matrices(proj, c("GeneScoreMatrix", "PeakMatrix"))
same_coordinates <- function(x, y) {
  identical(as.character(GenomicRanges::seqnames(x)), as.character(GenomicRanges::seqnames(y))) &&
    identical(BiocGenerics::start(x), BiocGenerics::start(y)) &&
    identical(BiocGenerics::end(x), BiocGenerics::end(y))
}
if (!same_coordinates(SummarizedExperiment::rowRanges(original_matrices$PeakMatrix), original_peaks)) {
  stop("Input PeakMatrix does not match the active peak set; repair peaks first.")
}

# Relocate retained imputation weights in the private project copy.
weights <- ArchR::getImputeWeights(proj)
if (length(weights$Weights) > 0) {
  for (i in seq_along(weights$Weights)) {
    target <- file.path(project_dir, "ImputeWeights", basename(weights$Weights[[i]]))
    if (!file.exists(target)) stop("Missing retained imputation weight file: ", target)
    weights$Weights[[i]] <- target
  }
  proj@imputeWeights <- weights
}

message("Rebuilding motifs on the existing ", length(proj@peakSet), " peaks.")
proj <- add_motif_annotations(proj, genome)
proj <- ArchR::addBgdPeaks(proj, force = TRUE)
proj <- ArchR::addDeviationsMatrix(
  proj, peakAnnotation = "Motif", matrixName = "MotifMatrix",
  out = c("z", "deviations"), force = TRUE
)
validate_motif_arrow_files(ArchR::getArrowFiles(proj))

# Save in place on the copied project so genes, peaks and weights are retained.
ArchR::saveArchRProject(proj, outputDirectory = project_dir, load = FALSE)
rebase_saved_archr_project(project_dir)
saved <- ArchR::loadArchRProject(project_dir, force = FALSE, showLogo = FALSE)
validate_motif_arrow_files(ArchR::getArrowFiles(saved))
stopifnot(
  identical(GenomicRanges::ranges(original_peaks), GenomicRanges::ranges(saved@peakSet)),
  identical(GenomicRanges::seqnames(original_peaks), GenomicRanges::seqnames(saved@peakSet)),
  identical(ArchR::getCellColData(saved), original_cells),
  identical(saved@reducedDims, original_dims),
  identical(saved@embeddings, original_embeddings)
)
matrices <- read_recovery_matrices(saved, c("GeneScoreMatrix", "PeakMatrix", "MotifMatrix"))
for (name in names(original_matrices)) {
  if (!identical(SummarizedExperiment::rowRanges(original_matrices[[name]]),
                 SummarizedExperiment::rowRanges(matrices[[name]])) ||
      !identical(SummarizedExperiment::assays(original_matrices[[name]]),
                 SummarizedExperiment::assays(matrices[[name]]))) {
    stop(name, " changed during motif recovery.")
  }
}
# Verify full cell coverage in every motif matrix, not just the sampled reads.
for (sample in names(ArchR::getArrowFiles(saved))) {
  arrow <- ArchR::getArrowFiles(saved)[[sample]]
  available <- rhdf5::h5read(arrow, "MotifMatrix/Info/CellNames")
  expected <- sub("^[^#]+#", "", rownames(original_cells)[original_cells$Sample == sample])
  if (!all(expected %in% available)) stop("MotifMatrix is missing cells for ", sample)
}
matches <- readRDS(ArchR::getPeakAnnotation(saved, "Motif")$Matches)
if (!same_coordinates(SummarizedExperiment::rowRanges(matches), saved@peakSet)) {
  stop("Restored motif annotations do not match the retained peaks.")
}
if (!setequal(SummarizedExperiment::rowData(matrices$MotifMatrix)$name, colnames(matches))) {
  stop("Saved motif scores do not match the restored annotations.")
}
write.csv(
  data.frame(genome = genome, peaks = length(saved@peakSet),
             motifs = nrow(matrices$MotifMatrix), cells = nrow(original_cells)),
  file.path(dirname(project_dir), "motif_recovery_summary.csv"), row.names = FALSE
)
message("Saved and validated motif recovery project: ", project_dir)
