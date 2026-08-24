library("ArchR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Mmusculus.UCSC.mm39")
library("BSgenome.Rnorvegicus.UCSC.rn6")
library("TxDb.Mmusculus.UCSC.mm39.knownGene")
library("org.Mm.eg.db")

source("/root/wf/R/archr.R")

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(
    paste(
      "Usage: archr_motif_peaks.R <input_project> <output_project>",
      "<genome> <include_y_chromosome> <genome_size> <threads>"
    )
  )
}

input_project <- args[[1]]
output_project <- args[[2]]
genome <- args[[3]]
include_y_chromosome <- parse_bool_arg(args[[4]])
genome_size <- as.numeric(args[[5]])
num_threads <- max(1, as.integer(args[[6]]))

ArchR::addArchRThreads(threads = num_threads)
proj <- ArchR::loadArchRProject(input_project, force = TRUE, showLogo = FALSE)
exclude_chr <- archr_exclude_chroms(include_y_chromosome)

message("Calling reproducible cluster peaks with ", num_threads, " worker(s).")
proj <- ArchR::addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Clusters",
  pathToMacs2 = ArchR::findMacs2(),
  genomeSize = genome_size,
  excludeChr = exclude_chr,
  maxPeaks = 300000,
  force = TRUE
)
proj <- ArchR::addPeakMatrix(ArchRProj = proj, force = TRUE)
proj <- add_motif_annotations(proj, genome)

message("Saving annotated motif peak checkpoint to: ", output_project)
ArchR::saveArchRProject(
  ArchRProj = proj,
  outputDirectory = output_project,
  load = FALSE
)
