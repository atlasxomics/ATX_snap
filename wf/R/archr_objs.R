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
args <- commandArgs(trailingOnly = TRUE)
print(args)

project_name <- args[1]
genome <- args[2]
metadata_path <- args[3]
tile_size <- 5000
min_tss <- 0  # Use filtering from SnapATAC2
min_frags <- 0  # Use filtering from SnapATAC2
lsi_iterations <- 2
lsi_resolution <- 0.7
lsi_varfeatures <- 25000
num_threads <- 50

runs <- strsplit(args[4:length(args)], ",")
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

proj <- create_archrproject( # from archr.R
  inputs, genome, min_tss, min_frags, tile_size, out_dir
)

# Add Conditions to CellColData ----
for (run in runs) {
  proj$Condition[proj$Sample == run[1]] <- run[3]
}
saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Copy over filtering from SnapATAC2
tryCatch({
  obs <- read.csv(metadata_path)
}, error = function(e) {
  message("Could not read obs.csv: ", e$message)
})

tryCatch({
  cells <- obs[["X"]]
}, error = function(e) {
  cells <- proj$cellName
  message(
    "Expected column 'cluster' not found in obs.csv; skipping filtering...",
    e$message
  )
})

tryCatch({
  clusters <- obs[["cluster"]]
}, error = function(e) {
  clusters <- rep("C0", length(cells))
  message(
    "No cluster information found in obs; assinging cluster 0",
    e$message
  )
})

proj <- proj[proj$cellNames %in% cells]

# Add clusters from SnapATAC2
proj <- addCellColData(
  ArchRProj = proj,
  data = clusters,
  cells = cells,
  name = "Clusters",
  force = TRUE
)

# Parse conditions into 'treatments', add as columns to CellColData ----
conds <- strsplit(proj$Condition, split = "\\s|-")

for (i in seq_along(conds[[1]])) {
  proj <- ArchR::addCellColData(
    proj,
    data = extract_nth_ele(conds, i), # from utils.R
    name = paste0("condition_", i),
    cells = proj$cellNames,
    force = TRUE
  )
}
treatment <- names(ArchR::getCellColData(proj))[
  grep("condition_", names(ArchR::getCellColData(proj)))
]

print(paste("Treatments:", treatment))

# Set global values for project data ----
n_samples <- length(unique(proj$Sample))
n_cond <- length(unique(proj$Condition))
n_cells <- length(proj$cellNames)

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Dimensionality reduction and clustering ----
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = lsi_iterations,
  clusterParams = list(
    resolution = c(lsi_resolution),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = lsi_varfeatures,
  dimsToUse = 1:30,
  force = TRUE
)

proj <- ArchR::addImputeWeights(proj)

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Gene Expression Analysis ----------------------------------------------------

# Create SeuratObjects for gene matrix ----
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

# Extract gene matrix for SeuratObject --
gene_matrix <- ArchR::getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  asMatrix = TRUE
)

# Identify empty features for filtering volcano plots --
print("Identifying empty features...")
gsm_mat <- SummarizedExperiment::assay(gene_matrix, "GeneScoreMatrix")
empty_feat_idx <- which(rowSums(gsm_mat) == 0)
empty_feat <- SummarizedExperiment::rowData(gene_matrix)$name[empty_feat_idx]

rm(gsm_mat)
gc()

matrix <- ArchR::imputeMatrix(
  mat = assay(gene_matrix),
  imputeWeights = getImputeWeights(proj)
)

gene_row_names <- gene_matrix@elementMetadata$name
rownames(matrix) <- gene_row_names

rm(gene_matrix)
gc()

# Create and save SeuratObjects --
print("Creating SeuratObjects...")

seurat_objs <- c()
for (run in runs) {

  obj <- build_atlas_seurat_object( # from seurat.R
    run_id = run[1],
    matrix = matrix,
    metadata = metadata,
    spatial_path = run[5]
  )

  saveRDS(obj, file = paste0(run[1], "_SeuratObj.rds"))
  seurat_objs <- c(seurat_objs, obj)
}

print("Available SeuratObjects:")
seurat_objs

# # Identify marker genes ----
# # Marker genes per cluster, save marker gene rds, csv, heatmap.csv --
# n_clust <- length(unique(proj$Clusters))

# cluster_marker_genes <- get_marker_genes( # from archr.R
#   proj,
#   group_by = "Clusters",
#   markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
#   heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
#   rd_name = "IterativeLSI"
# )

# saveRDS(cluster_marker_genes$markers_gs, "markersGS_clusters.rds")
# write.csv(
#   cluster_marker_genes$marker_list,
#   "ranked_genes_per_cluster.csv",
#   row.names = FALSE
# )
# write.csv(cluster_marker_genes$heatmap_gs, "genes_per_cluster_hm.csv")

# # Marker genes per sample, save marker gene rds, csv, heatmap.csv --
# if (n_samples > 1) {

#   sample_marker_genes <- get_marker_genes(
#     proj,
#     group_by = "Sample",
#     markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
#     heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
#     rd_name = "IterativeLSI"
#   )

#   saveRDS(sample_marker_genes$markers_gs, "markersGS_sample.rds")
#   write.csv(
#     sample_marker_genes$marker_list,
#     "ranked_genes_per_sample.csv",
#     row.names = FALSE
#   )
#   write.csv(sample_marker_genes$heatmap_gs, "genes_per_sample_hm.csv")

# }

# # Marker genes per treatment, save marker gene rds, csv, heatmap.csv --
# if (n_cond > 1) {

#   for (i in seq_along(treatment)) {

#     treatment_marker_genes <- get_marker_genes(
#       proj,
#       group_by = treatment[i],
#       markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
#       heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
#       rd_name = "IterativeLSI"
#     )

#     saveRDS(
#       treatment_marker_genes$markers_gs,
#       paste0("markersGS_condition_", i, ".rds")
#     )
#     write.csv(
#       treatment_marker_genes$marker_list,
#       paste0("ranked_genes_per_condition", i, ".csv"),
#       row.names = FALSE
#     )
#     write.csv(
#       treatment_marker_genes$heatmap_gs,
#       paste0("genes_per_conditions", i, "_hm.csv")
#     )
#   }
# }

# # Volcano plots for genes ----
# if (n_cond > 1) {
#   for (j in seq_along(treatment)) {

#     # Get gene markers df for all clusters together --
#     marker_genes_df <- get_marker_df(
#       proj = proj,
#       group_by = treatment[j],
#       matrix = "GeneScoreMatrix",
#       seq_names = NULL,
#       max_cells = n_cells,  # Equals total cells in project
#       test_method = "ttest"
#     )

#     # Create a merged marker genes df for clusters for which no condition is
#     # >90% of all cells --
#     req_clusters <- get_required_clusters(proj, treatment[j])
#     marker_genes_by_cluster_df <- get_marker_df_clusters(
#       proj, req_clusters, treatment[j]
#     )

#     # Per condition, merge dfs and cleanup data --
#     conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
#     for (cond in conditions) {

#       volcano_table <- get_volcano_table( # from archr.R
#         marker_genes_df, marker_genes_by_cluster_df, cond, "gene", empty_feat
#       )

#       write.table(
#         volcano_table,
#         paste0(
#           "volcanoMarkers_genes_", j, "_", cond, ".csv"
#         ),
#         sep = ",",
#         quote = FALSE,
#         row.names = FALSE
#       )
#       print(paste0("volcanoMarkers_genes_", j, "_", cond, ".csv is done!"))

#       features <- unique(volcano_table$cluster)
#       others <- paste(conditions[conditions != cond], collapse = "|")
#       volcano_plots <- list()
#       for (i in seq_along(features)) {
#         volcano_plots[[i]] <- scvolcano(
#           volcano_table, cond, others, features[[i]]
#         )
#       }

#       pdf(paste0("volcano_plots_", cond, ".pdf"))
#       for (plot in volcano_plots) {
#         print(plot)
#       }
#       dev.off()
#     }
#   }
# } else {
#   print("There are not enough conditions to be compared with!")
# }

# saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Peak and motif analysis -----------------------------------------------------

# Due to mcapply bug, too many threads in peak calling can lead to OOM error;
# decrease theads here, per specification from input parameters ----
addArchRThreads(threads = num_threads)

# Peak calling and motif enrichment for clusters ----
proj <- get_annotated_peaks(proj, "Clusters", genome_size, genome)

# saveArchRProject(ArchRProj = proj, archrproj_dir)

# # # Save run metrics in medians.csv ----
# # medians <- get_proj_medians(proj)
# # write.csv(medians, file = "medians.csv", row.names = FALSE)

# # # Get marker peaks for clusters, samples, treatments; save as csv ----

# # # Initialize base data frame, set significance cutoff --
# # peak_data <- data.frame(proj@peakSet@ranges, proj@peakSet@elementMetadata)
# # cut_off <- "Pval <= 0.05 & Log2FC >= 0.1"

# # # Marker peaks per clusters --
# # marker_peaks_c <- get_marker_peaks(proj, "Clusters", peak_data, cut_off)

# # write.csv(
# #   marker_peaks_c$marker_peak_list,
# #   "marker_peaks_per_cluster.csv",
# #   row.names = FALSE
# # )

# # write.csv(
# #   marker_peaks_c$total_peaks,
# #   "complete_peak_list_cluster.csv",
# #   row.names = FALSE
# # )

# # # Marker peaks per sample --
# # marker_peaks_s <- get_marker_peaks(proj, "Sample", peak_data, cut_off)

# # write.csv(
# #   marker_peaks_s$marker_peak_list,
# #   "marker_peaks_per_sample.csv",
# #   row.names = FALSE
# # )

# # write.csv(
# #   marker_peaks_s$total_peaks,
# #   "complete_peak_list_sample.csv",
# #   row.names = FALSE
# # )

# # # Marker peaks per treatment --
# # if (n_cond > 1) {

# #   for (i in seq_along(treatment)) {

# #     marker_peaks_t <- get_marker_peaks(proj, treatment[i], peak_data, cut_off)

# #     write.csv(
# #       marker_peaks_t$marker_peak_list,
# #       file = paste0("marker_peaks_per_condition-", i, ".csv"),
# #       row.names = FALSE
# #     )

# #     write.csv(
# #       marker_peaks_t$total_peaks,
# #       file = paste0("complete_peak_list_condition-", i, ".csv"),
# #       row.names = FALSE

# #     )
# #   }
# # }

# # # Get enriched motif data for clusters, write to disk ----
# # enriched_motifs_c <- get_enriched_motifs(
# #   proj, marker_peaks_c$marker_peaks, cut_off
# # )

# # write.csv(enriched_motifs_c$enrich_df, "enrichedMotifs_cluster.csv")
# # saveRDS(enriched_motifs_c$enrich_motifs, "enrichMotifs_clusters.rds")
# # write.csv(enriched_motifs_c$heatmap_em, "motif_per_cluster_hm.csv")


# Create motif SeuratObjects ----
# Create motif count matrix --
proj <- addBgdPeaks(proj, force = TRUE)

proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(ArchRProj = proj, outputDirectory = "adddevs_ArchRProject")
print("sleeping...")
Sys.sleep(7200)

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
  # motifs <- paste0("z:", motifs)
  motifs <- unique(motifs)

  # proj <- addImputeWeights(proj)

  dev_score <- getDeviation_ArchR(
    ArchRProj = proj,
    name = motifs,
    useAssay = "z",  # can be "deviations" or "z"ls
    log2Norm = TRUE
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

###
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
  saveRDS(enriched_motifs_s$enrich_motifs, "enrichMotifs_sample.rds")
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
    saveRDS(
      enriched_motifs_t$enrich_motifs,
      paste0("enrichMotifs_condition_", i, ".rds")
    )
    write.csv(
      enriched_motifs_t$heatmap_em,
      paste0("motif_per_condition_", i, "_hm.csv")
    )
  }
}

# Plot gene, peak, motif heatmaps by cluster, generate heatmaps .pdf ----
# Initiate heatmaps list and plot marker gene per cluster heatmap --
heatmaps <- list()

# Recompute gene heatmap for plotting (transpose, plotLog2FC different) --
heatmap_gs_plotting <- plotMarkerHeatmap(
  seMarker = cluster_marker_genes$markers_gs,
  cutOff = cut_off,
  transpose = TRUE
)

# Save for plotting with peaks and motifs --
gene_hm <- ComplexHeatmap::draw(
  heatmap_gs_plotting,
  heatmap_legend_side = "bot",
  annotation_legend_side = "bot",
  column_title = paste0("Marker genes (", cut_off, ")"),
  column_title_gp = gpar(fontsize = 12)
)
heatmaps[[1]] <- gene_hm

# Create peak heatmap by cluster --
heatmap_peaks <- plotMarkerHeatmap(
  seMarker = marker_peaks_c$marker_peaks,
  cutOff = cut_off,
  transpose = TRUE
)

peak_hm <- ComplexHeatmap::draw(
  heatmap_peaks,
  heatmap_legend_side = "bot",
  annotation_legend_side = "bot",
  column_title = paste0("Marker peaks (", cut_off, ")"),
  column_title_gp = gpar(fontsize = 12)
)

heatmaps[[2]] <- peak_hm

# Create motif heatmap by cluster --
heatmap_plot <- plotEnrichHeatmap(
  enriched_motifs_c$enrich_motifs,
  transpose = TRUE,
  n = 50,
  cutOff = 2
)

heatmap_motifs <- ComplexHeatmap::draw(
  heatmap_plot,
  heatmap_legend_side = "bot",
  column_title = paste0("Marker motifs (", cut_off, ")"),
  column_title_gp = gpar(fontsize = 12)
)

heatmaps[[3]] <- heatmap_motifs

# Save heatmaps to disk as pdf (per cluster: genes, peaks, motifs) --
print("Saving heatmap plots...")

pdf("heatmaps_all.pdf")
for (i in seq_along(heatmaps)) {
  print(heatmaps[[i]])
}
dev.off()

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
      test_method = "wilcoxon"
    )

    # Create a data from of marker genes for clusters for which no condition is
    # >90% of all cells --
    req_clusters <- get_required_clusters(proj, treatment[j]) # from archr.R
    marker_motifs_by_cluster_df <- get_marker_df_clusters(
      proj, req_clusters, treatment[j]
    )

    # Merge and cleanup data --
    conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
    for (cond in conditions) {

      volcano_table <- get_volcano_table( # from archr.R
        marker_motifs_df,
        marker_motifs_by_cluster_df,
        cond,
        "motif",
        empty_feat_m
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
          volcano_table,  cond, others, features_m[[i]]
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

all <- rename_cells(seurat_objs)
all_m <- rename_cells(seurat_objs_m)

# Convert Seurat to h5ad and save ----
for (obj in all) {
  seurat_to_h5ad(obj, FALSE, paste0(unique(obj$Sample), "_g"))  # from utils.R
}

# Convert Seurat to h5ad and save ----
for (obj in all_m) {
  seurat_to_h5ad(obj, FALSE, paste0(unique(obj$Sample), "_m"))  # from utils.R
}
