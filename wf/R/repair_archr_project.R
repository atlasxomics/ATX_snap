# Run in the workflow image: Rscript /root/wf/R/repair_archr_project.R
# <input_project> <new_output_project> <genome> [group_by] [include_y] [threads]
library("ArchR")
source("/root/wf/R/load_genomes.R")
source("/root/wf/R/archr.R")
source("/root/wf/R/complete_archr.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 6) {
  stop("Usage: repair_archr_project.R <input> <new_output> <genome> [group_by] [include_y] [threads]")
}
input <- normalizePath(args[[1]], mustWork = TRUE)
output <- args[[2]]
genome <- args[[3]]
group_by <- if (length(args) >= 4) args[[4]] else ""
include_y <- if (length(args) >= 5) parse_bool_arg(args[[5]]) else FALSE
threads <- if (length(args) >= 6) as.integer(args[[6]]) else 8L
sizes <- c(hg38 = 3.3e9, mm10 = 3.0e9, mm39 = 2.7e9, rnor6 = 2.9e9)
if (!genome %in% names(sizes)) stop("Unsupported genome: ", genome)
if (file.exists(output)) stop("Output must be a new project directory.")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
output <- file.path(normalizePath(dirname(output)), basename(output))
if (startsWith(output, paste0(input, "/"))) stop("Output cannot be inside the input project.")
dir.create(output)
files <- list.files(input, full.names = TRUE, all.files = TRUE, no.. = TRUE)
if (!all(file.copy(files, output, recursive = TRUE))) stop("Failed to copy input project.")
set.seed(42)
ArchR::addArchRThreads(threads = threads)
proj <- ArchR::loadArchRProject(output, force = TRUE, showLogo = FALSE)
proj <- rebase_group_coverage_paths(proj, output)

gene_ok <- tryCatch({
  validate_arrow_matrices(ArchR::getArrowFiles(proj), "GeneScoreMatrix")
  # Exercise the data reader as well as feature metadata before preserving it.
  md <- ArchR::getCellColData(proj)
  cells <- rownames(md)[!duplicated(md$Sample)]
  ArchR::getMatrixFromProject(proj[cells, ], useMatrix = "GeneScoreMatrix", threads = 1)
  TRUE
}, error = function(e) {
  message("Rebuilding incomplete GeneScoreMatrix: ", e$message)
  FALSE
})
if (!gene_ok) {
  # Original createArrowFiles used the default gene-score model.
  proj <- ArchR::addGeneScoreMatrix(
    proj, excludeChr = archr_exclude_chroms(include_y), force = TRUE
  )
}

if (length(proj@peakSet) == 0 || nzchar(group_by)) {
  md <- ArchR::getCellColData(proj)
  if (!nzchar(group_by)) {
    treatments <- grep("condition_", colnames(md), value = TRUE)
    group_by <- if ("Condition" %in% colnames(md) &&
                    length(unique(md$Condition)) > 1 && length(treatments) > 0) {
      tail(treatments, 1)
    } else if (length(unique(md$Sample)) > 1) "Sample" else "Clusters"
  }
  if (!group_by %in% colnames(md)) stop("Missing grouping column: ", group_by)
  proj <- ArchR::addCellColData(
    proj, data = as.character(md[[group_by]]), cells = proj$cellNames,
    name = group_by, force = TRUE
  )
  proj <- get_annotated_peaks(proj, group_by, sizes[[genome]], genome, include_y)
} else {
  message("Preserving the existing active peak set; rebuilding its PeakMatrix.")
  proj <- ArchR::addPeakMatrix(proj, force = TRUE)
}
proj <- rebuild_final_motifs(proj, genome)
save_complete_archr_project(proj, output)
message("Validated complete project saved to: ", output)
