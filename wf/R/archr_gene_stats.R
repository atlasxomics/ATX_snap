library("ArchR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Rnorvegicus.UCSC.rn6")
library("ComplexHeatmap")
library("dplyr")
library("GenomicRanges")
library("grid")
library("gridExtra")
library("harmony")
library("plyr")
library("Seurat")
library("SummarizedExperiment")
library("tidyverse")

source("/root/wf/R/archr.R")
source("/root/wf/R/seurat.R")
source("/root/wf/R/utils.R")

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
print(args)

project_name <- args[1]
archrproj_path <- args[2]
output_dir <- args[3]
marker_threads <- max(1, as.integer(args[4]))
marker_max_cells <- 500

run_args <- if (length(args) < 5) character(0) else args[5:length(args)]
runs <- if (length(run_args) == 0) list() else strsplit(run_args, ",")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

sample_name_map <- c()
for (run in runs) {
  run_id <- run[1]
  sample_name <- if (length(run) >= 2) trimws(run[2]) else ""
  if (!is.na(sample_name) && nzchar(sample_name) &&
      tolower(sample_name) != "none") {
    sample_name_map[run_id] <- sample_name
  }
}

remap_sample_ids <- function(values, sample_name_map) {
  remapped <- values
  for (i in seq_along(values)) {
    value <- as.character(values[[i]])
    if (value %in% names(sample_name_map)) {
      remapped[[i]] <- sample_name_map[[value]]
    }
  }
  return(remapped)
}

get_gene_feature_names <- function(gene_matrix) {
  feature_names <- gene_matrix@elementMetadata$name
  if (is.null(feature_names)) {
    feature_names <- SummarizedExperiment::rowData(gene_matrix)$name
  }
  feature_names
}

ArchR::addArchRThreads(threads = marker_threads)

message("Loading ArchRProject from: ", archrproj_path)
proj <- ArchR::loadArchRProject(
  path = archrproj_path,
  force = TRUE,
  showLogo = FALSE
)

if ("Condition" %in% names(ArchR::getCellColData(proj))) {
  conds <- strsplit(proj$Condition, split = "\\s|-")
  for (i in seq_along(conds[[1]])) {
    proj <- ArchR::addCellColData(
      proj,
      data = extract_nth_ele(conds, i),
      name = paste0("condition_", i),
      cells = proj$cellNames,
      force = TRUE
    )
  }
}

treatment <- names(ArchR::getCellColData(proj))[
  grep("condition_", names(ArchR::getCellColData(proj)))
]

if (length(sample_name_map) == 0 && "sample_name" %in% names(ArchR::getCellColData(proj))) {
  sample_df <- as.data.frame(ArchR::getCellColData(proj, select = c("Sample", "sample_name")))
  sample_df <- sample_df[!is.na(sample_df$sample_name), , drop = FALSE]
  sample_df <- sample_df[nzchar(sample_df$sample_name), , drop = FALSE]
  sample_df <- sample_df[tolower(sample_df$sample_name) != "none", , drop = FALSE]
  if (nrow(sample_df) > 0) {
    sample_df <- unique(sample_df)
    sample_name_map <- stats::setNames(sample_df$sample_name, sample_df$Sample)
  }
}

print(paste("Treatments:", paste(treatment, collapse = ", ")))
print(paste(
  "Using maxCells =",
  marker_max_cells,
  "and threads =",
  marker_threads,
  "for ArchR gene marker tests."
))

n_samples <- length(unique(proj$Sample))
n_cond <- length(unique(proj$Condition))

# Identify empty features for filtering volcano plots.
print("Identifying empty gene features...")
gene_matrix <- ArchR::getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  asMatrix = TRUE
)
gene_row_names <- get_gene_feature_names(gene_matrix)
empty_feat_idx <- which(Matrix::rowSums(
  SummarizedExperiment::assay(gene_matrix, "GeneScoreMatrix")
) == 0)
empty_feat <- gene_row_names[empty_feat_idx]
print(paste("Found", length(empty_feat), "empty features"))
rm(gene_matrix, gene_row_names, empty_feat_idx)
gc(verbose = FALSE, full = TRUE)

# Marker genes per cluster, save marker gene csv, heatmap csv, and heatmap pdf.
cluster_marker_genes <- get_marker_genes(
  proj,
  group_by = "Clusters",
  markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
  heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
  max_cells = marker_max_cells,
  threads = marker_threads
)

write.csv(
  cluster_marker_genes$marker_list,
  file.path(tables_dir, "ranked_genes_per_cluster.csv"),
  row.names = FALSE
)
write.csv(
  cluster_marker_genes$heatmap_gs,
  file.path(tables_dir, "genes_per_cluster_hm.csv")
)

cut_off <- "Pval <= 0.05 & Log2FC >= 0.1"

if (
  is.null(cluster_marker_genes$markers_gs) ||
    nrow(cluster_marker_genes$heatmap_gs) == 0 ||
    ncol(cluster_marker_genes$heatmap_gs) == 0
) {
  message("Skipping gene heatmap plot because no cluster marker genes were available.")
  empty_pdf(file.path(figures_dir, "heatmap_genes.pdf"), "No cluster marker genes available")
} else {
  heatmap_gs_plotting <- plotMarkerHeatmap(
    seMarker = cluster_marker_genes$markers_gs,
    cutOff = cut_off,
    transpose = TRUE
  )

  gene_hm <- ComplexHeatmap::draw(
    heatmap_gs_plotting,
    heatmap_legend_side = "bot",
    annotation_legend_side = "bot",
    column_title = paste0("Marker genes (", cut_off, ")"),
    column_title_gp = gpar(fontsize = 12)
  )

  print("Saving gene heatmap...")
  pdf(file.path(figures_dir, "heatmap_genes.pdf"))
  print(gene_hm)
  dev.off()
}
rm(list = intersect(c("cluster_marker_genes", "heatmap_gs_plotting", "gene_hm"), ls()))
gc(verbose = FALSE, full = TRUE)

# Marker genes per sample.
if (n_samples > 1) {
  sample_marker_genes <- get_marker_genes(
    proj,
    group_by = "Sample",
    markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
    heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
    max_cells = marker_max_cells,
    threads = marker_threads
  )

  write.csv(
    sample_marker_genes$marker_list,
    file.path(tables_dir, "ranked_genes_per_sample.csv"),
    row.names = FALSE
  )
  write.csv(
    sample_marker_genes$heatmap_gs,
    file.path(tables_dir, "genes_per_sample_hm.csv")
  )

  if (length(sample_name_map) > 0) {
    sample_marker_list_named <- sample_marker_genes$marker_list
    marker_groups <- names(sample_marker_list_named)
    if (!is.data.frame(sample_marker_list_named) && !is.null(marker_groups)) {
      names(sample_marker_list_named) <- remap_sample_ids(
        marker_groups, sample_name_map
      )
    }
    write.csv(
      sample_marker_list_named,
      file.path(tables_dir, "ranked_genes_per_sample_name.csv"),
      row.names = FALSE
    )

    sample_heatmap_named <- sample_marker_genes$heatmap_gs
    if (!is.null(colnames(sample_heatmap_named))) {
      colnames(sample_heatmap_named) <- remap_sample_ids(
        colnames(sample_heatmap_named), sample_name_map
      )
    }
    if (!is.null(rownames(sample_heatmap_named))) {
      rownames(sample_heatmap_named) <- remap_sample_ids(
        rownames(sample_heatmap_named), sample_name_map
      )
    }
    write.csv(
      sample_heatmap_named,
      file.path(tables_dir, "genes_per_sample_name_hm.csv")
    )
  }

  rm(list = intersect(
    c("sample_marker_genes", "sample_marker_list_named", "sample_heatmap_named"),
    ls()
  ))
  gc(verbose = FALSE, full = TRUE)
}

# Marker genes per treatment.
if (n_cond > 1) {
  for (i in seq_along(treatment)) {
    treatment_marker_genes <- get_marker_genes(
      proj,
      group_by = treatment[i],
      markers_cutoff = "FDR <= 1 & Log2FC >= -Inf",
      heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
      max_cells = marker_max_cells,
      threads = marker_threads
    )

    write.csv(
      treatment_marker_genes$marker_list,
      file.path(tables_dir, paste0("ranked_genes_per_condition_", i, ".csv")),
      row.names = FALSE
    )
    write.csv(
      treatment_marker_genes$heatmap_gs,
      file.path(tables_dir, paste0("genes_per_condition_", i, "_hm.csv"))
    )
    rm(treatment_marker_genes)
    gc(verbose = FALSE, full = TRUE)
  }
}

# Volcano plots for genes.
if (n_cond > 1) {
  for (j in seq_along(treatment)) {
    marker_genes_df <- get_marker_df(
      proj = proj,
      group_by = treatment[j],
      matrix = "GeneScoreMatrix",
      seq_names = NULL,
      max_cells = marker_max_cells,
      test_method = "ttest",
      diff_metric = "Log2FC",
      threads = marker_threads
    )

    req_clusters <- get_required_clusters(proj, treatment[j])
    marker_genes_by_cluster_df <- get_marker_df_clusters(
      proj = proj,
      clusters = req_clusters,
      group_by = treatment[j],
      seq_names = "z",
      matrix = "GeneScoreMatrix",
      max_cells = marker_max_cells,
      test_method = "ttest",
      diff_metric = "Log2FC",
      threads = marker_threads
    )

    conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
    for (cond in conditions) {
      volcano_table <- get_volcano_table(
        marker_genes_df,
        marker_genes_by_cluster_df,
        cond,
        "gene",
        empty_feat,
        fc_col = "Log2FC"
      )

      volcano_csv <- file.path(
        tables_dir,
        paste0("volcanoMarkers_genes_", j, "_", cond, ".csv")
      )
      write.table(
        volcano_table,
        volcano_csv,
        sep = ",",
        quote = FALSE,
        row.names = FALSE
      )
      print(paste0(volcano_csv, " is done!"))

      features <- unique(volcano_table$cluster)
      others <- paste(conditions[conditions != cond], collapse = "|")
      volcano_plots <- list()
      for (i in seq_along(features)) {
        volcano_plots[[i]] <- scvolcano(
          volcano_table, cond, others, features[[i]], fc_col = "Log2FC"
        )
      }

      pdf(file.path(figures_dir, paste0("volcano_plots_", cond, ".pdf")))
      for (plot in volcano_plots) {
        print(plot)
      }
      dev.off()
      rm(volcano_table, volcano_plots)
      gc(verbose = FALSE, full = TRUE)
    }
    rm(marker_genes_df, marker_genes_by_cluster_df)
    gc(verbose = FALSE, full = TRUE)
  }
} else {
  print("There are not enough conditions to be compared with!")
}

print("Gene differential statistics complete.")
