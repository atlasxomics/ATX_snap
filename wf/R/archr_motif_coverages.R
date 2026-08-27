library("ArchR")
source("/root/wf/R/load_genomes.R")

source("/root/wf/R/archr.R")

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    "Usage: archr_motif_coverages.R <input_project> <output_project> <threads>"
  )
}

input_project <- args[[1]]
output_project <- args[[2]]
num_threads <- max(1, as.integer(args[[3]]))

ArchR::addArchRThreads(threads = num_threads)
proj <- ArchR::loadArchRProject(input_project, force = TRUE, showLogo = FALSE)

# Keep labels such as "-1" as literal group names rather than numeric indices.
if (!is.character(proj$Clusters)) {
  proj <- ArchR::addCellColData(
    ArchRProj = proj,
    data = as.character(proj$Clusters),
    cells = proj$cellNames,
    name = "Clusters",
    force = TRUE
  )
}

message(
  "Generating cluster group coverages with ",
  num_threads,
  " ArchR worker(s)."
)
proj <- ArchR::addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Clusters",
  maxCells = 1500,
  force = TRUE
)

message("Saving motif coverage checkpoint to: ", output_project)
ArchR::saveArchRProject(
  ArchRProj = proj,
  outputDirectory = output_project,
  load = FALSE
)
rebase_saved_archr_project(output_project)
