#' Rountines for manipulation of and plotting with ArchR

library("ArchR")
library("ggplot2")
library("ggplot2")
library("gridExtra")
library("harmony")
library("patchwork")
library("SummarizedExperiment")
library("S4Vectors")

source("/root/wf/R/utils.R")


safe_get_marker_features <- function(..., context = "marker features") {
  tryCatch(
    ArchR::getMarkerFeatures(...),
    error = function(e) {
      message(
        "Skipping ", context,
        " because ArchR could not construct background matches: ",
        e$message
      )
      return(NULL)
    }
  )
}


get_group_levels <- function(proj, group_by) {
  groups <- sort(unique(proj@cellColData[group_by][, 1]))
  groups[!is.na(groups)]
}


empty_result_df <- function() {
  data.frame(stringsAsFactors = FALSE)
}


empty_marker_feature_table <- function(diff_metric = "Log2FC") {
  marker_df <- data.frame(
    metric = numeric(0),
    p_val = numeric(0),
    p_val_adj = numeric(0),
    gene = character(0),
    cluster = character(0),
    stringsAsFactors = FALSE
  )
  colnames(marker_df)[1] <- diff_metric
  marker_df
}


empty_pdf <- function(filename, label) {
  pdf(filename)
  plot.new()
  text(0.5, 0.5, label)
  dev.off()
}


add_motif_annotations <- function(proj, genome) {
  #' Wrapper for ArchR::addMotifAnnotations()

  motif_set <- list("mm10" = "cisbp", "hg38" = "cisbp", "rnor6" = "encode")
  species <- list(
    "mm10" = NULL, "hg38" = NULL, "rnor6" = ArchR::getGenome(ArchRProj = proj)
  )

  proj <- ArchR::addMotifAnnotations(
    ArchRProj = proj,
    motifSet = motif_set[[genome]],
    name = "Motif",
    force = TRUE,
    species = species[[genome]]
  )
  return(proj)
}

add_lsi <- function(proj, iterations, var_features) {
  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = iterations,
    clusterParams = list(
      resolution = c(resolution), sampleCells = 10000, n.start = 10
    ),
    varFeatures = var_features,
    dimsToUse = 1:30,
    force = TRUE
  )
  return(proj)
}


create_archrproject <- function(
  inputs, genome, min_tss, min_frags, tile_size, out_dir
) {
  #' Create ArrowFiles and ArchRProject; handles mm10, hg38, rnor6.  Inputs are
  #' a named vector mapping run_id to local fragment files path.

  if (genome %in% c("mm10", "hg38")) {

    ArchR::addArchRGenome(genome)
    geneAnnotation <- ArchR::getGeneAnnotation()
    genomeAnnotation <- ArchR::getGenomeAnnotation()

  } else if (genome == "rnor6") {

    # load in geneAnnotation and genomeAnotation from custom_ArchR repo
    load(
    "/root/custom_ArchR_genomes_and_annotations/rn6/rn6_liftoff_mm10NcbiRefSeq_ArchR_annotations.rda"
    )

  } else {

    stop(
      "Genome not one for 'mm10', 'hg38', 'rnor6'; please supply correct
      genome."
    )
  }

  arrow_files <- ArchR::createArrowFiles(
    inputFiles = inputs,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
    sampleNames = names(inputs),
    minTSS = min_tss,
    minFrags = min_frags,
    maxFrags = 1e+07,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    offsetPlus = 0,
    offsetMinus = 0,
    TileMatParams = list(tileSize = tile_size)
  )

  proj <- ArchR::ArchRProject(
    ArrowFiles = arrow_files,
    outputDirectory = out_dir,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation
  )

  return(proj)
}

get_annotated_peaks <- function(proj, group_by, genome_size, genome) {

  proj <- ArchR::addGroupCoverages(
    ArchRProj = proj,
    groupBy = group_by,
    maxCells = 1500,
    force = TRUE
  )

  proj <- ArchR::addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = group_by,
    pathToMacs2 = ArchR::findMacs2(),
    genomeSize = genome_size,
    maxPeaks = 300000,
    force = TRUE
  )
  proj <- ArchR::addPeakMatrix(proj, force = TRUE)

  # Add motif annotations -----
  proj <- add_motif_annotations(proj, genome)

  return(proj)
}

get_enriched_motifs <- function(proj, marker_peaks, cutoff) {

  empty_enrich_df <- data.frame(
    group_name = character(0),
    stringsAsFactors = FALSE
  )
  empty_heatmap_df <- data.frame(
    feature = character(0),
    value = numeric(0),
    stringsAsFactors = FALSE
  )

  if (is.null(marker_peaks)) {
    message(
      "Skipping motif enrichment for cutoff ",
      cutoff,
      " because no marker peaks were available."
    )
    return(list(
      enrich_df = empty_enrich_df,
      enrich_motifs = NULL,
      heatmap_em = empty_heatmap_df,
      has_enrichment = FALSE
    ))
  }

  enrich_motifs <- tryCatch(
    ArchR::peakAnnoEnrichment(
      seMarker = marker_peaks,
      ArchRProj = proj,
      peakAnnotation = "Motif",
      cutOff = cutoff
    ),
    error = function(e) {
      message(
        "peakAnnoEnrichment failed for cutoff ",
        cutoff,
        "; returning empty motif enrichment tables. ",
        e$message
      )
      return(NULL)
    }
  )

  if (is.null(enrich_motifs)) {
    return(list(
      enrich_df = empty_enrich_df,
      enrich_motifs = NULL,
      heatmap_em = empty_heatmap_df,
      has_enrichment = FALSE
    ))
  }

  enrich_df <- tryCatch(
    data.frame(enrich_motifs@assays@data),
    error = function(e) {
      message(
        "Could not extract motif enrichment assay data; returning empty table. ",
        e$message
      )
      empty_enrich_df
    }
  )

  motif_lst <- unique(rownames(enrich_motifs))
  if (length(motif_lst) > 0) {
    split_string <- strsplit(motif_lst, split = "\\(")

    req_motifs <- gsub("_", "-", extract_nth_ele(split_string)) # from utils.R
    req_motifs <- gsub(" ", "", req_motifs)

    rownames(enrich_motifs) <- req_motifs
  }

  heatmap_em <- tryCatch(
    ArchR::plotEnrichHeatmap(
      enrich_motifs, n = 50, transpose = FALSE, returnMatrix = TRUE, cutOff = 2
    ),
    error = function(e) {
      if (grepl("No enrichments found", e$message, fixed = TRUE)) {
        message(
          "No motif enrichments found for cutoff ",
          cutoff,
          "; returning empty motif heatmap table."
        )
        return(empty_heatmap_df)
      }
      message(
        "plotEnrichHeatmap failed; returning empty motif heatmap table. ",
        e$message
      )
      return(empty_heatmap_df)
    }
  )

  has_enrichment <- !is.null(heatmap_em) &&
    nrow(heatmap_em) > 0 &&
    ncol(heatmap_em) > 0

  if (has_enrichment) {
    motif_lst <- unique(rownames(heatmap_em))
    split_string <- strsplit(motif_lst, split = "\\(")

    req_motifs <- gsub("_", "-", extract_nth_ele(split_string))  # from utils.R
    req_motifs <- gsub(" ", "", req_motifs)

    rownames(heatmap_em) <- req_motifs
    heatmap_em <- as.data.frame(heatmap_em)
  }

  return(list(
    enrich_df = enrich_df,
    enrich_motifs = enrich_motifs,
    heatmap_em = heatmap_em,
    has_enrichment = has_enrichment
  ))
}

get_marker_df <- function(
  proj, group_by, matrix, seq_names, max_cells, test_method, diff_metric
) {
  #' Return data frame of ArhcR marker features with:
  #' "avg_log2FC", "p_val", "p_val_adj", "gene", "cluster"

  conditions <- get_group_levels(proj, group_by)
  markers_df <- setNames(
    lapply(conditions, function(x) empty_marker_feature_table(diff_metric)),
    conditions
  )

  # Get markers -----
  markers <- safe_get_marker_features(
    ArchRProj = proj,
    useMatrix = matrix,
    groupBy = group_by,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useSeqnames = seq_names,
    maxCells = max_cells,
    normBy = "none",
    testMethod = test_method,
    context = paste0("marker dataframe for ", group_by)
  )
  if (is.null(markers)) {
    return(markers_df)
  }

  # Create data frame with for MeanDiff, Pval, FDR -----
  markers_df1 <- tryCatch(
    SummarizedExperiment::assay(markers, diff_metric),
    error = function(e) NULL
  )
  markers_df2 <- tryCatch(
    SummarizedExperiment::assay(markers, "Pval"),
    error = function(e) NULL
  )
  markers_df3 <- tryCatch(
    SummarizedExperiment::assay(markers, "FDR"),
    error = function(e) NULL
  )
  if (is.null(markers_df1) || is.null(markers_df2) || is.null(markers_df3)) {
    message("Unable to extract marker assays for ", group_by, "; returning empty tables.")
    return(markers_df)
  }

  for (conds in conditions) {

    markers_df[[conds]] <- S4Vectors::DataFrame(
      markers_df1[[conds]],
      markers_df2[[conds]],
      markers_df3[[conds]]
    )
    markers_df[[conds]] <- as.data.frame(markers_df[[conds]])
    markers_df[[conds]]$genes <- SummarizedExperiment::rowData(markers)$name
    markers_df[[conds]]$cluster <- rep(
      "All", length(rownames(markers_df[[conds]]))
    )
    colnames(markers_df[[conds]]) <- c(
      diff_metric, "p_val", "p_val_adj", "gene", "cluster"
    )
  }
  return(markers_df)
}

get_marker_df_clusters <- function(
  proj, clusters, group_by, matrix, seq_names, test_method, diff_metric
) {
  conditions <- get_group_levels(proj, group_by)
  markers_df_by_cluster <- setNames(
    lapply(conditions, function(x) empty_marker_feature_table(diff_metric)),
    conditions
  )

  for (cluster in clusters) {

    proj_filter <- BiocGenerics::which(proj$Clusters == cluster)
    cells_subset <- proj$cellNames[proj_filter]
    proj_subset <- proj[cells_subset, ]
    n_cells <- length(proj_subset$cellNames)

    markers_by_cluster <- safe_get_marker_features(
      ArchRProj = proj_subset,
      useMatrix = matrix,
      groupBy = group_by,
      bias = c("TSSEnrichment", "log10(nFrags)"),
      maxCells = n_cells,
      useSeqnames = seq_names,
      normBy = "none",
      testMethod = test_method,
      context = paste0("cluster-level marker dataframe for ", group_by, " in ", cluster)
    )
    if (is.null(markers_by_cluster)) {
      next
    }

    markerlist_df1 <- tryCatch(
      SummarizedExperiment::assay(markers_by_cluster, diff_metric),
      error = function(e) NULL
    )
    markerlist_df2 <- tryCatch(
      SummarizedExperiment::assay(markers_by_cluster, "Pval"),
      error = function(e) NULL
    )
    markerlist_df3 <- tryCatch(
      SummarizedExperiment::assay(markers_by_cluster, "FDR"),
      error = function(e) NULL
    )
    if (is.null(markerlist_df1) || is.null(markerlist_df2) || is.null(markerlist_df3)) {
      message(
        "Unable to extract cluster-level marker assays for ",
        group_by,
        " in ",
        cluster,
        "; skipping this cluster."
      )
      next
    }

    cluster_conditions <- intersect(colnames(markerlist_df1), conditions)
    for (conds in cluster_conditions) {
      marker_df <- S4Vectors::DataFrame(
        markerlist_df1[, conds],
        markerlist_df2[, conds],
        markerlist_df3[, conds]
      )
      marker_df <- as.data.frame(marker_df)
      marker_df$genes <- SummarizedExperiment::rowData(markers_by_cluster)$name
      marker_df$cluster <- rep(cluster, dim(marker_df)[1])
      colnames(marker_df) <- c(
        diff_metric, "p_val", "p_val_adj", "gene", "cluster"
      )
      markers_df_by_cluster[[conds]] <- rbind(
        markers_df_by_cluster[[conds]],
        marker_df
      )
    }
  }

  return(markers_df_by_cluster)
}

get_marker_genes <- function(
  proj, group_by, markers_cutoff, heatmap_cutoff
) {

  markers_gs <- safe_get_marker_features(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = group_by,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "ttest",
    context = paste0("marker genes for ", group_by)
  )
  if (is.null(markers_gs)) {
    return(
      list(
        markers_gs = NULL,
        marker_list = empty_result_df(),
        heatmap_gs = empty_result_df()
      )
    )
  }

  marker_list <- tryCatch(
    ArchR::getMarkers(markers_gs, cutOff = markers_cutoff),
    error = function(e) {
      message("Unable to extract marker genes for ", group_by, ": ", e$message)
      empty_result_df()
    }
  )
  marker_pvals <- tryCatch(
    SummarizedExperiment::assay(markers_gs, "Pval"),
    error = function(e) NULL
  )

  if (!is.null(marker_pvals) && !is.data.frame(marker_list)) {
    marker_names <- SummarizedExperiment::rowData(markers_gs)$name
    marker_groups <- names(marker_list)
    if (is.null(marker_groups)) {
      marker_groups <- seq_along(marker_list)
    }

    for (group in marker_groups) {
      marker_group <- marker_list[[group]]
      if (nrow(marker_group) == 0) {
        marker_list[[group]]$Pval <- numeric(0)
        next
      }

      pval_col <- if (is.character(group)) {
        match(group, colnames(marker_pvals))
      } else {
        as.integer(group)
      }

      if (is.na(pval_col)) {
        next
      }

      marker_list[[group]]$Pval <- marker_pvals[
        match(marker_group$name, marker_names), pval_col
      ]
    }
  }

  heatmap_gs <- tryCatch(
    ArchR::plotMarkerHeatmap(
      seMarker = markers_gs,
      cutOff = heatmap_cutoff,
      plotLog2FC = TRUE,
      returnMatrix = TRUE
    ),
    error = function(e) {
      message(
        "Unable to generate marker heatmap matrix for ",
        group_by,
        ": ",
        e$message
      )
      empty_result_df()
    }
  )
  return(
    list(
      markers_gs = markers_gs,
      marker_list = marker_list,
      heatmap_gs = heatmap_gs
    )
  )
}

get_marker_peaks <- function(proj, group_by, peak_data, cut_off) {

  marker_peaks <- safe_get_marker_features(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = group_by,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    k = 100,
    testMethod = "wilcoxon",
    context = paste0("marker peaks for ", group_by)
  )
  if (is.null(marker_peaks)) {
    return(list(
      marker_peaks = NULL,
      marker_peak_list = empty_result_df(),
      total_peaks = empty_result_df()
    ))
  }

  marker_peak_list <- tryCatch(
    ArchR::getMarkers(marker_peaks, cutOff = cut_off),
    error = function(e) {
      message("Unable to extract marker peaks for ", group_by, ": ", e$message)
      empty_result_df()
    }
  )

  # Merge all peaks with significant marker peaks, write to csv -----
  total_peaks <- if (
    is.data.frame(marker_peak_list) &&
      (!all(c("start", "end") %in% colnames(marker_peak_list)))
  ) {
    empty_result_df()
  } else {
    merge(peak_data, marker_peak_list, by = c("start", "end"))
  }

  return(list(
    marker_peaks = marker_peaks,
    marker_peak_list = marker_peak_list,
    total_peaks = total_peaks
  ))
}

get_proj_medians <- function(proj) {
  #' Get dataframe with median TSS, nFrags, FRIP value per run.

  medians <- ArchR::getCellColData(ArchRProj = proj)
  tss <- aggregate(
    medians@listData$TSSEnrichment,
    by = list(medians@listData$Sample),
    FUN = median
  )
  nfrags <- aggregate(
    medians@listData$nFrags,
    by = list(medians@listData$Sample),
    FUN = median
  )
  conditions <- aggregate(
    medians@listData$Condition,
    by = list(medians@listData$Sample),
    FUN = max
  )
  frip <- aggregate(
    medians@listData$FRIP,
    by = list(medians@listData$Sample),
    FUN = median
  )
  frip$x <- round(frip$x, 4)
  list_dfs <- list(tss, nfrags, frip, conditions)
  merged_df <- Reduce(
    function(x, y) merge(x, y, by = "Group.1", all = TRUE), list_dfs
  )
  names(merged_df) <- c(
    "run_id", "median_TSS", "median_fragments", "median_FRIP", "condition"
  )

  return(merged_df)
}

get_required_clusters <- function(proj, group_by) {

  req_df <- as.data.frame(ArchR::getCellColData(proj))
  df1 <- table(req_df$Clusters, req_df[, group_by])
  distr <- as.data.frame.matrix(round(prop.table(as.matrix(df1), 1), 2))

  lst <- list()
  for (i in 1:nrow(distr)) {
    row <- distr[i, ]
    if (
      sum(unname(unlist(row)) >= 0.90) == 1
    ) {
      lst[[i]] <- rownames(row)
    }
  }
  not_req_list <- unlist(lst)

  req_clusters <- unique(proj$Clusters)
  req_clusters <- req_clusters[order(as.numeric(gsub("C", "", req_clusters)))]
  req_clusters <- req_clusters[which(!req_clusters %in% not_req_list)]

  return(req_clusters)
}

get_topn_hm_feats <- function(heatmap, n_groups, n_feats) {
  #' Return a vector with top n features from heatmap for each grouping
  #' (cluster, samples, conditions); vector will have length n_groups * nfeats.

  lst <- list()
  for (i in seq_along(1: n_groups)) {
    lst[[i]] <- heatmap[, c(1, i + 1)]
    lst[[i]] <- lst[[i]][
      order(lst[[i]][, 2], decreasing = TRUE),
    ][1:n_feats, 1]
  }

  # Flatten to vector, remove duplicated features, NAs -----
  req_feats <- unlist(lst)
  req_feats <- req_feats[!duplicated(req_feats)]
  req_feats <- na.omit(req_feats)

  return(req_feats)
}

get_volcano_table <- function(
  markers_df,
  markers_by_cluster_df,
  condition,
  feature,
  empty_feat,
  fc_col = "avg_log2FC"
) {

  # Merge df with all clusters with df for each cluster -----
  marker_df <- markers_df[[condition]]
  if (is.null(marker_df)) {
    marker_df <- empty_marker_feature_table(fc_col)
  }
  marker_df_cluster <- markers_by_cluster_df[[condition]]
  if (is.null(marker_df_cluster)) {
    marker_df_cluster <- empty_marker_feature_table(fc_col)
  }
  merged_df <- rbind(
    marker_df, marker_df_cluster
  )

  # Remove empty features -----
  if (!is.null(empty_feat)) {
    merged_df <- merged_df[which(!merged_df$gene %in% empty_feat), ]
  }

  # Remove NA values -----
  merged_df <- na.omit(merged_df)

  if (nrow(merged_df) == 0) {
    return(merged_df)
  }

  # Check if fc_col exists
  if (!fc_col %in% colnames(merged_df)) {
    stop(paste("Column", fc_col, "not found in merged_df. Available columns:",
               paste(colnames(merged_df), collapse = ", ")))
  }

  # Convert fold change column to numeric if needed -----
  if (!is.numeric(merged_df[[fc_col]])) {
    converted_values <- as.numeric(as.character(merged_df[[fc_col]]))

    # Safety check - only assign if conversion was successful
    if (length(converted_values) > 0 && !all(is.na(converted_values))) {
      merged_df[[fc_col]] <- converted_values
    } else {
      warning("Could not convert fc_col to numeric - keeping original values")
    }
  }
  # Remove FDR equal to 0 -----
  merged_df <- merged_df[which(!merged_df$p_val_adj == 0), ]

  # Get string of other conditions
  others <- paste(
    names(markers_df)[condition != names(markers_df)], collapse = "|"
  )

  return(merged_df)
}

plot_umap <- function(proj, name) {
  plot <- ArchR::plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = name,
    embedding = "UMAP"
  )
  return(plot)
}

save_umap <- function(proj, color_by) {
  #' Save a plots of UMAP embeddings as a single .pdf, color by each feature in
  #' character vector 'color_by'.

  umap_plots <- c()
  for (feature in color_by) {
    umap <- plot_umap(proj, feature)
    umap_plots[[feature]] <- umap
  }

  pdf("umap_plots.pdf")
  gridExtra::grid.arrange(grobs = umap_plots, ncol = 2)
  dev.off()
}

performPairwiseMarkerComparisons <- function(
    ArchRProj,
    groupBy = "Condition",
    useMatrix = "GeneScoreMatrix",
    testMethod = "ttest",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    maxCells = 500,
    scaleTo = 10^4,
    threads = getArchRThreads(),
    outputDir = "PairwiseMarkers",
    saveResults = TRUE,
    removeSameSignReciprocals = TRUE,  # New parameter
    verbose = TRUE
) {
  
  if(saveResults && !dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }
  
  allGroups <- unique(getCellColData(ArchRProj, groupBy, drop = TRUE))
  allGroups <- allGroups[!is.na(allGroups)]
  
  if(verbose) {
    cat("Found", length(allGroups), "groups:", paste(allGroups, collapse = ", "), "\n")
  }
  
  # Generate all directional pairwise comparisons (A vs B and B vs A)
  directionalCombos <- expand.grid(group1 = allGroups, group2 = allGroups, stringsAsFactors = FALSE)
  directionalCombos <- directionalCombos[directionalCombos$group1 != directionalCombos$group2, ]
  
  if(verbose) {
    cat("Performing", nrow(directionalCombos), "directional pairwise comparisons\n")
  }
  
  pairwiseResults <- list()
  rawResults <- list()  # Store raw results for filtering
  
  for(i in seq_len(nrow(directionalCombos))) {
    
    group1 <- directionalCombos$group1[i]
    group2 <- directionalCombos$group2[i]
    comparisonName <- paste0(group1, "_vs_", group2)
    
    if(verbose) {
      cat("Processing comparison", i, "of", nrow(directionalCombos), ":", comparisonName, "\n")
    }
    
    tryCatch({
      markerFeatures <- getMarkerFeatures(
        ArchRProj = ArchRProj,
        groupBy = groupBy,
        useGroups = group1,
        bgdGroups = group2,
        useMatrix = useMatrix,
        bias = bias,
        testMethod = testMethod,
        maxCells = maxCells,
        scaleTo = scaleTo,
        threads = threads,
        verbose = FALSE
      )
      
      markers_df <- processMarkerFeatures(markerFeatures, group1, ArchRProj)
      
      pairwiseResults[[comparisonName]] <- markerFeatures
      rawResults[[comparisonName]] <- markers_df
      
    }, error = function(e) {
      if(verbose) {
        cat("Error in comparison", comparisonName, ":", e$message, "\n")
      }
      pairwiseResults[[comparisonName]] <- NULL
      rawResults[[comparisonName]] <- NULL
    })
  }
  
  # Filter same-signed reciprocals if requested
  if(removeSameSignReciprocals && length(rawResults) > 0) {
    if(verbose) {
      cat("Filtering same-signed reciprocal features...\n")
    }
    
    filteredResults <- filterSameSignReciprocals(rawResults, verbose = verbose)
    
    # Save filtered results
    if(saveResults) {
      for(comparisonName in names(filteredResults)) {
        if(!is.null(filteredResults[[comparisonName]]) && !is.null(comparisonName)) {
          splitResult <- strsplit(comparisonName, "_vs_")[[1]]
          if(length(splitResult) >= 2) {
            group1 <- splitResult[1]
            group2 <- splitResult[2]
            filename <- paste0(group1, "-vs-", group2, "_filtered.csv")
            write.csv(
              filteredResults[[comparisonName]], 
              file.path(outputDir, filename), 
              row.names = FALSE
            )
          }
        }
      }
      
      # Save summary of filtering
      filteringSummary <- generateFilteringSummary(rawResults, filteredResults)
      write.csv(
        filteringSummary,
        file.path(outputDir, "filtering_summary.csv"),
        row.names = FALSE
      )
    }
  }
  
  # Save raw results regardless
  if(saveResults) {
    for(comparisonName in names(rawResults)) {
      if(!is.null(rawResults[[comparisonName]]) && !is.null(comparisonName)) {
        splitResult <- strsplit(comparisonName, "_vs_")[[1]]
        if(length(splitResult) >= 2) {
          group1 <- splitResult[1]
          group2 <- splitResult[2]
          filename <- paste0(group1, "-vs-", group2, "_raw.csv")
          write.csv(
            rawResults[[comparisonName]], 
            file.path(outputDir, filename), 
            row.names = FALSE
          )
        }
      }
    }
    
    saveRDS(pairwiseResults, file.path(outputDir, "all_directional_pairwise_comparisons.rds"))
  }
  
  if(verbose) {
    cat("Completed all directional pairwise comparisons!\n")
    cat("Results saved to:", outputDir, "\n")
  }

  # Return filtered results if filtering was applied, otherwise raw results
  if(removeSameSignReciprocals && exists("filteredResults")) {
    return(list(raw = pairwiseResults, filtered = filteredResults))
  } else {
    return(pairwiseResults)
  }
}

# Function to filter same-signed reciprocal features
filterSameSignReciprocals <- function(rawResults, verbose = TRUE) {
  
  filteredResults <- list()
  comparisonNames <- names(rawResults)
  
  # Track features removed for summary
  removalStats <- data.frame(
    comparison = character(),
    total_features = integer(),
    removed_features = integer(),
    remaining_features = integer(),
    stringsAsFactors = FALSE
  )

  for (comparisonName in comparisonNames) {

    if (is.null(rawResults[[comparisonName]])) next

    # Parse comparison name
    groups <- strsplit(comparisonName, "_vs_")[[1]]
    group1 <- groups[1]
    group2 <- groups[2]

    # Find reciprocal comparison
    reciprocalName <- paste0(group2, "_vs_", group1)
    
    if(reciprocalName %in% comparisonNames && !is.null(rawResults[[reciprocalName]])) {
      
      df1 <- rawResults[[comparisonName]]
      df2 <- rawResults[[reciprocalName]]
      
      # Find common genes
      commonGenes <- intersect(df1$gene, df2$gene)
      
      if(length(commonGenes) > 0) {
        
        # Get log2FC for common genes
        df1_common <- df1[df1$gene %in% commonGenes, ]
        df2_common <- df2[df2$gene %in% commonGenes, ]
        
        # Match by gene name
        df1_common <- df1_common[order(df1_common$gene), ]
        df2_common <- df2_common[order(df2_common$gene), ]
        
        # Find same-signed features (both positive or both negative)
        sameSigned <- sign(df1_common$avg_log2FC) == sign(df2_common$avg_log2FC)
        sameSigned <- sameSigned & !is.na(sameSigned)
        
        # Remove same-signed features
        genesToRemove <- df1_common$gene[sameSigned]
        df1_filtered <- df1[!df1$gene %in% genesToRemove, ]
        
        filteredResults[[comparisonName]] <- df1_filtered
        
        if(verbose) {
          cat("Comparison", comparisonName, ": removed", sum(sameSigned), 
              "same-signed features out of", nrow(df1), "total features\n")
        }
        
        # Store removal stats
        removalStats <- rbind(removalStats, data.frame(
          comparison = comparisonName,
          total_features = nrow(df1),
          removed_features = sum(sameSigned),
          remaining_features = nrow(df1_filtered),
          stringsAsFactors = FALSE
        ))
        
      } else {
        # No common genes, keep all features
        filteredResults[[comparisonName]] <- rawResults[[comparisonName]]
        
        removalStats <- rbind(removalStats, data.frame(
          comparison = comparisonName,
          total_features = nrow(rawResults[[comparisonName]]),
          removed_features = 0,
          remaining_features = nrow(rawResults[[comparisonName]]),
          stringsAsFactors = FALSE
        ))
      }
      
    } else {
      # No reciprocal comparison found, keep all features
      filteredResults[[comparisonName]] <- rawResults[[comparisonName]]
      
      removalStats <- rbind(removalStats, data.frame(
        comparison = comparisonName,
        total_features = nrow(rawResults[[comparisonName]]),
        removed_features = 0,
        remaining_features = nrow(rawResults[[comparisonName]]),
        stringsAsFactors = FALSE
      ))
    }
  }

  # Store removal stats as attribute
  attr(filteredResults, "removal_stats") <- removalStats

  return(filteredResults)
}

# Function to generate filtering summary
generateFilteringSummary <- function(rawResults, filteredResults) {
  
  removalStats <- attr(filteredResults, "removal_stats")
  
  if(is.null(removalStats)) {
    return(data.frame(
      comparison = names(rawResults),
      total_features = sapply(rawResults, nrow),
      removed_features = 0,
      remaining_features = sapply(rawResults, nrow),
      removal_percentage = 0
    ))
  }
  
  removalStats$removal_percentage <- round(
    (removalStats$removed_features / removalStats$total_features) * 100, 2
  )
  
  return(removalStats)
}

# Analysis code for processing marker features (unchanged)
processMarkerFeatures <- function(markers, groups, proj) {
  df <- S4Vectors::DataFrame(
    markers_df1 <- SummarizedExperiment::assay(markers, "MeanDiff"),
    markers_df2 <- SummarizedExperiment::assay(markers, "Pval"),
    markers_df3 <- SummarizedExperiment::assay(markers, "FDR"),
    markers_df4 <- SummarizedExperiment::assay(markers, "Log2FC"),
    markers_df5 <- SummarizedExperiment::assay(markers, "Mean")
  )
  df$genes <- SummarizedExperiment::rowData(markers)$name
  colnames(df) <- c(
    "mean_diff", "p_val", "p_val_adj", "avg_log2FC", "mean", "gene"
  )
  return(df)
}
