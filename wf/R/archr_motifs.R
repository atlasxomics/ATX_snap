library("ArchR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Rnorvegicus.UCSC.rn6")
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
archrproj_path <- args[4]
num_threads <- 50

runs <- strsplit(args[5:length(args)], ",")
runs

inputs <- c()  # Inputs for ArrowFiles (run_id : fragment_file path)
for (run in runs) {
  inputs[run[1]] <- run[2]
}
inputs

out_dir <- paste0(project_name, "_ArchRProject")

# Set genome size for peak calling ----
genome_sizes <- list("hg38" = 3.3e+09, "mm10" = 3.0e+09, "rnor6" = 2.9e+09)
genome_size <- genome_sizes[[genome]]

tempdir <- "/root"
archrproj_dir <- paste0(project_name, "_ArchRProject")

# Create ArchRProject ---------------------------------------------------------
addArchRThreads(threads = 50)

proj <- loadArchRProject(archrproj_path)

rd_names <- names(proj@reducedDims@listData)
if (length(rd_names) == 0) {
  warning("No reduced dims found for ArchRProjects; adding LSI.")
  proj <- add_lsi(proj, 2, 0.7, 25000)
}
if ("Spectral" %in% names(proj@reducedDims@listData)) {
  rd_name <- "Spectral"
} else {
  rd_name <- "IterativeLSI"
}

# Add Conditions to CellColData ----
for (run in runs) {
  proj$Condition[proj$Sample == run[1]] <- run[3]
}

# Parse conditions into 'treatments', add as columns to CellColData ----
conds <- strsplit(proj$Condition, split = "\\s|-")
treatment <- names(ArchR::getCellColData(proj))[
  grep("condition_", names(ArchR::getCellColData(proj)))
]

print(paste("Treatments:", treatment))

# Set global values for project data ----
n_samples <- length(unique(proj$Sample))
n_cond <- length(unique(proj$Condition))
n_cells <- length(proj$cellNames)

# Extract metadata for Seurat object --
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

# Peak and motif analysis -----------------------------------------------------

# Due to mcapply bug, too many threads in peak calling can lead to OOM error;
# decrease theads here, per specification from input parameters ----
addArchRThreads(threads = num_threads)

# Peak calling and motif enrichment for clusters ----
proj <- get_annotated_peaks(proj, "Clusters", genome_size, genome)

saveArchRProject(ArchRProj = proj, archrproj_dir)

# Save run metrics in medians.csv ----
medians <- get_proj_medians(proj)
write.csv(medians, file = "medians.csv", row.names = FALSE)

# Get marker peaks for clusters, samples, treatments; save as csv ----

# Initialize base data frame, set significance cutoff --
peak_data <- data.frame(proj@peakSet@ranges, proj@peakSet@elementMetadata)
cut_off <- "Pval <= 0.05 & Log2FC >= 0.1"

# Marker peaks per clusters --
marker_peaks_c <- get_marker_peaks(proj, "Clusters", peak_data, cut_off)

write.csv(
  marker_peaks_c$marker_peak_list,
  "marker_peaks_per_cluster.csv",
  row.names = FALSE
)

write.csv(
  marker_peaks_c$total_peaks,
  "complete_peak_list_cluster.csv",
  row.names = FALSE
)

# Marker peaks per sample --
marker_peaks_s <- get_marker_peaks(proj, "Sample", peak_data, cut_off)

write.csv(
  marker_peaks_s$marker_peak_list,
  "marker_peaks_per_sample.csv",
  row.names = FALSE
)

write.csv(
  marker_peaks_s$total_peaks,
  "complete_peak_list_sample.csv",
  row.names = FALSE
)

# Marker peaks per treatment --
if (n_cond > 1) {

  for (i in seq_along(treatment)) {

    marker_peaks_t <- get_marker_peaks(proj, treatment[i], peak_data, cut_off)

    write.csv(
      marker_peaks_t$marker_peak_list,
      file = paste0("marker_peaks_per_condition-", i, ".csv"),
      row.names = FALSE
    )

    write.csv(
      marker_peaks_t$total_peaks,
      file = paste0("complete_peak_list_condition-", i, ".csv"),
      row.names = FALSE

    )
  }
}

# Get enriched motif data for clusters, write to disk ----
enriched_motifs_c <- get_enriched_motifs(
  proj, marker_peaks_c$marker_peaks, cut_off
)

write.csv(enriched_motifs_c$enrich_df, "enrichedMotifs_cluster.csv")
write.csv(enriched_motifs_c$heatmap_em, "motif_per_cluster_hm.csv")

# Create motif SeuratObjects ----
# Create motif count matrix --
proj <- addBgdPeaks(proj, force = TRUE)

proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)

markers_motifs <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "MotifMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = "z"
)

# Get deviation matrix for significant motifs --
marker_motifs_list <- getMarkers(
  markers_motifs, cutOff = "FDR < 0.9 & MeanDiff >= 0"
)

motifs <- list()
for (i in seq_len(length(marker_motifs_list))) {
  if (length(marker_motifs_list[[i]]$name) > 1) {
    motifs <- c(motifs, marker_motifs_list[[i]]$name)
  }
}

if (length(motifs) > 1) {
  motifs <- unlist(motifs)
  motifs <- paste0("z:", motifs)
  motifs <- unique(motifs)

  proj <- addImputeWeights(proj, reducedDims = "Spectral")

  dev_score <- getDeviation_ArchR(
    ArchRProj = proj,
    name = motifs,
    imputeWeights = getImputeWeights(proj)
  )

  dev_score[is.na(dev_score)] <- 0
}

dev_score2 <- dev_score[!is.infinite(rowSums(dev_score)), ]
colnames(dev_score2) <- rownames(getCellColData(proj))

# Remove 0 deviations per All samples --
all_zero <- names(which(rowSums(dev_score2) == 0))
dev_score2 <- dev_score2[which(!rownames(dev_score2) %in% c(all_zero)), ]

# Convert to dgCmatrix --
dev_score3 <- Matrix(as.matrix(dev_score2), sparse = TRUE)

empty_feat_idx_m <- which(rowSums(dev_score3) == 0)
empty_feat_m <- rownames(dev_score3)[empty_feat_idx_m]

# Create motif seurat objects --
seurat_objs_m <- c()
for (run in runs) {

  obj <- build_atlas_seurat_object(
    run_id = run[1],
    matrix = dev_score3,
    metadata = metadata,
    spatial_path = run[5]
  )

  saveRDS(obj, file = paste0(run[1], "_SeuratObjMotif.rds"))
  seurat_objs_m <- c(seurat_objs_m, obj)
}

print("Available seurat_objMotifs:")
seurat_objs_m

# Peak calling and motifs for Sample ----
if (n_samples > 1) {

  proj <- get_annotated_peaks(proj, "Sample", genome_size, genome)

  saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

  # Repeat getMarkerPeaks for new peak-set, enriched motifs per sample --
  sample_marker_peaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Sample",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    k = 100,
    testMethod = "wilcoxon"
  )

  enriched_motifs_s <- get_enriched_motifs(
    proj, sample_marker_peaks, cut_off
  )

  write.csv(enriched_motifs_s$enrich_df, "enrichedMotifs_sample.csv")
  write.csv(enriched_motifs_s$heatmap_em, "motif_per_sample_hm.csv")
}

# Peak calling and motif enrichment per treatment ----
if (n_cond > 1) {

  for (i in seq_along(treatment)) {

    proj <- get_annotated_peaks(proj, treatment[i], genome_size, genome)

    saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

    # Get marker peaks and Enriched motifs per treatment --
    treatment_marker_peaks <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = treatment[i],
      bias = c("TSSEnrichment", "log10(nFrags)"),
      k = 100,
      testMethod = "wilcoxon"
    )
    enriched_motifs_t <- get_enriched_motifs(
      proj, treatment_marker_peaks, cut_off
    )

    write.csv(
      enriched_motifs_t$enrich_df,
      paste0("enrichedMotifs_condition_", i, ".csv")
    )
    write.csv(
      enriched_motifs_t$heatmap_em,
      paste0("motif_per_condition_", i, "_hm.csv")
    )
  }
}

# Volcano plots for motifs ----
if (n_cond > 1) {
  for (j in seq_along(treatment)) {

    # Get motif markers for all clusters together --
    marker_motifs_df <- get_marker_df( # from archr.R
      proj = proj,
      group_by = treatment[j],
      matrix = "MotifMatrix",
      seq_names = "z",
      max_cells = 5000,
      test_method = "wilcoxon",
      diff_metric = "MeanDiff"
    )

    # Create a data from of marker genes for clusters for which no condition is
    # >90% of all cells --
    req_clusters <- get_required_clusters(proj, treatment[j]) # from archr.R
    marker_motifs_by_cluster_df <- get_marker_df_clusters(
      proj = proj,
      clusters = req_clusters,
      group_by =  treatment[j],
      seq_names = "z",
      matrix = "MotifMatrix",
      test_method = "ttest",
      diff_metric = "MeanDiff"
    )

    # Merge and cleanup data --
    conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
    for (cond in conditions) {

      volcano_table <- get_volcano_table( # from archr.R
        marker_motifs_df,
        marker_motifs_by_cluster_df,
        cond,
        "motif",
        empty_feat_m,
        fc_col = "MeanDiff"
      )
      write.table(
        volcano_table,
        paste0("volcanoMarkers_motifs_", j, "_", cond, ".csv"),
        sep = ",",
        quote = FALSE,
        row.names = FALSE
      )
      print(
        paste0("writing volcanoMarkers_motifs_", j, "_", cond, ".csv is done!")
      )

      features_m <- unique(volcano_table$cluster)
      others <- paste(conditions[conditions != cond], collapse = "|")
      volcano_plots_m <- list()
      for (i in seq_along(features_m)) {
        volcano_plots_m[[i]] <- scvolcano(
          volcano_table,  cond, others, features_m[[i]], fc_col = "MeanDiff"
        )
      }

      pdf(paste0("volcano_plots_motifs_", j, "_", cond, ".pdf"))
      for (plot in volcano_plots_m) {
        print(plot)
      }
      dev.off()
    }
  }
} else {
  print("There are not enough conditions to be compared with!")
}

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

all_m <- rename_cells(seurat_objs_m)

# Convert Seurat to h5ad and save ----
for (obj in all_m) {
  seurat_to_h5ad(obj, FALSE, paste0(unique(obj$Sample), "_m"))  # from utils.R
}
