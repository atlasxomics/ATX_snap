library("ArchR")
source("/root/wf/R/load_genomes.R")
source("/root/wf/R/archr.R")

set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: restore_peaks.R <project> <genome> <group_by> <include_y> <threads>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
genome <- args[[2]]
group_by <- args[[3]]
include_y <- parse_bool_arg(args[[4]])
ArchR::addArchRThreads(threads = as.integer(args[[5]]))
genome_sizes <- c(hg38 = 3.3e9, mm10 = 3.0e9, mm39 = 2.7e9, rnor6 = 2.9e9)
if (!genome %in% names(genome_sizes)) stop("Unsupported genome: ", genome)

select_final_peak_group <- function(metadata) {
  # Match the order in archr_motifs.R, including its metadata column order.
  treatments <- grep("condition_", colnames(metadata), value = TRUE)
  if ("Condition" %in% colnames(metadata) &&
      length(unique(metadata$Condition)) > 1 && length(treatments) > 0) {
    return(tail(treatments, 1))
  }
  if (length(unique(metadata$Sample)) > 1) return("Sample")
  "Clusters"
}

proj <- ArchR::loadArchRProject(project_dir, force = TRUE, showLogo = FALSE)
proj <- rebase_group_coverage_paths(proj, project_dir)
metadata <- ArchR::getCellColData(proj)
if (!nzchar(group_by)) group_by <- select_final_peak_group(metadata)
if (!group_by %in% colnames(metadata)) stop("Missing grouping column: ", group_by)
labels <- as.character(metadata[[group_by]])
if (anyNA(labels) || any(!nzchar(labels))) {
  stop("Grouping column contains missing or empty labels: ", group_by)
}
proj <- ArchR::addCellColData(
  proj, data = labels, cells = proj$cellNames, name = group_by, force = TRUE
)

message("Restoring the active peak set using group: ", group_by)
# Use the original coverage, MACS, PeakMatrix and motif annotation settings.
# This does not run motif deviations, enrichment or gene analysis.
proj <- get_annotated_peaks(
  proj, group_by, unname(genome_sizes[[genome]]), genome,
  include_y_chromosome = include_y
)
if (length(proj@peakSet) == 0) stop("Peak calling produced an empty peak set.")
output_root <- dirname(project_dir)
export_peak_beds_by_group(proj, group_by, output_root)
write.csv(
  data.frame(group_by = group_by, genome = genome, peaks = length(proj@peakSet)),
  file.path(output_root, "peak_recovery_summary.csv"), row.names = FALSE
)

# Save in place on the private copy, keeping all existing project files.
ArchR::saveArchRProject(proj, outputDirectory = project_dir, load = FALSE)
rebase_saved_archr_project(project_dir)
saved <- ArchR::loadArchRProject(project_dir, force = TRUE, showLogo = FALSE)
if (length(saved@peakSet) == 0 ||
    !"PeakMatrix" %in% ArchR::getAvailableMatrices(saved)) {
  stop("Saved recovery project is missing peaks or PeakMatrix.")
}
