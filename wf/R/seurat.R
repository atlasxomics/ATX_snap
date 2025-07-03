#' Rountines for manipulation of and plotting with Seurat

library("BPCells")
library("ggplot2")
library("ggrepel")
library("gridExtra")
library("harmony")
library("Matrix")
library("patchwork")
library("purrr")
library("Seurat")
library("SeuratDisk")

source("/root/wf/R/utils.R")

build_atlas_seurat_object <- function(
  run_id,
  matrix,
  metadata,
  spatial_path
) {
  #' Prepare and combine gene matrix, metadata, and image for SeuratObject
  #' for runs within a project.


  # First garbage collect to free up memory before starting
  gc(verbose = FALSE, full = TRUE)

  # Read the image first, as it's likely smaller than the matrix manipulation
  image <- Seurat::Read10X_Image(
    image.dir = spatial_path,
    filter.matrix = TRUE
  )

  # Filter metadata
  metadata <- metadata[metadata$Sample == run_id, ]

  # Find indices for subsetting instead of creating a temporary vector with grep
  col_indices <- grep(pattern = run_id, colnames(matrix))

  # Create sparse matrix directly from the subset to minimize memory usage
  message("Subsetting matrix and converting to sparse format...")

  # Option 1: If matrix is already in memory and very large
  matrix_sparse <- as(matrix[, col_indices], "dgCMatrix")

  # Immediately remove references to free memory
  rm(col_indices)
  # Don't remove matrix here as it's an input argument

  # Force garbage collection
  gc(verbose = FALSE, full = TRUE)

  # Set column names
  colnames(matrix_sparse) <- rownames(metadata)

  # Create assay object and immediately remove the source matrix
  message("Creating Seurat assay object...")
  matrix_assay <- Seurat::CreateAssayObject(counts = matrix_sparse)
  rm(matrix_sparse)
  gc(verbose = FALSE, full = TRUE)

  # Create Seurat object
  message("Creating Seurat object...")
  object <- Seurat::CreateSeuratObject(
    counts = matrix_assay,
    assay = "Spatial",
    meta.data = as.data.frame(metadata)
  )

  # Clean up
  rm(matrix_assay)
  gc(verbose = FALSE, full = TRUE)

  # Add image
  message("Adding spatial image...")
  image <- image[Seurat::Cells(x = object)]
  Seurat::DefaultAssay(object = image) <- "Spatial"
  object[["slice1"]] <- image

  # Final cleanup
  rm(image)
  gc(verbose = FALSE, full = TRUE)

  return(object)
}

combine_objs <- function(
  seurat_lst, samples, spatial, project_name
) {

  # ----------- Filter SeuratObjs  -----------
  # Remove cell with NA counts from each SeuratObj
  to_remove <- lapply(seurat_lst, function(x) {
    names(which(colSums(is.na(x@assays[[1]]@counts)) > 0))
  })
  filtered <- mapply(
    function(x, y) x[, !colnames(x) %in% y], seurat_lst, to_remove
  )

  # ----------- Prepare coordinates -----------
  # Make combinations of spatial coordinates, return as matrix; if n=1, simply
  # return df as matix.
  spatial_all <- lapply(seq_along(spatial), function(i) {
    tmp <- dplyr::bind_rows(spatial[-i])
    tmp$Spatial_1 <- 0
    tmp$Spatial_2 <- 0
    as.matrix(rbind(spatial[[i]], tmp))
  })

  # ----------- Prepare Metadata -----------
  # remove assay-specific suffixes in obj meta.data (nCount_Spatial -> nCount)
  filtered <- lapply(filtered, function(x) {
    colnames(x@meta.data) <- gsub(
      paste0("_", Seurat::Assays(x)),
      "",
      colnames(x@meta.data)
    )
    x
  })

  # get the list of metadata from seurat objects
  list_of_metadata <- lapply(filtered, function(x) {
    x@meta.data
  })

  # combine meta data
  meta.data <- do.call("rbind", list_of_metadata)
  write.csv(meta.data, "req_meta_data.csv", row.names = TRUE)

  # ----------- Create Combined SeuratObj -----------

  ncells  <- as.double(sum(sapply(filtered, ncol)))
  nfeats <- as.double(nrow(seurat_lst[[1]]))
  matrix_size <- ncells * nfeats

  if (matrix_size < 2^31 - 1) {

    print("Feature matrix less than 2^31 -1...")

    # Convert list of SeuratObjs to list of feat counts as data frames.
    filtered_dfs <- lapply(filtered, function(x) {
      df <- as.data.frame(x@assays[[1]]@counts)
      colnames(df) <- Seurat::Cells(x)
      df$region <- rownames(df)
      df
    })

    # Merge counts dfs, drop region
    combined_mat <- purrr::reduce(filtered_dfs, full_join, by = "region")
    rownames(combined_mat) <- combined_mat$region
    combined_mat$region <- NULL

    combined <- Seurat::CreateSeuratObject(
      counts = combined_mat,
      assay = "scATAC",
      meta.data = meta.data
    )

  } else {

    print("Feature matrix greater than 2^31 -1; using BPCells...")

    # Use BPCells to handle large matrices
    # Convert SeuratObjs to BPCells files
    counts_mat <- c()
    for (i in seq_along(filtered)) {

      path <- paste0(project_name, "/", samples[[i]], "_BP")

      BPCells::write_matrix_dir(
        mat = filtered[[i]]@assays[["Spatial"]]@counts,
        dir = path,
        overwrite = TRUE
      )
      counts_mat[[i]] <- BPCells::open_matrix_dir(dir = path)
    }

    # Create combiend SeuratObject
    combined <- Seurat::CreateSeuratObject(
      counts = counts_mat,
      assay = "scATAC",
      meta.data = meta.data
    )

    combined <- SeuratObject::JoinLayers(combined)
  }

  # Normalize and calculate variable features
  combined <- Seurat::NormalizeData(
    combined, normalization.method = "LogNormalize", scale.factor = 10000
  )

  combined <- Seurat::FindVariableFeatures(
    combined, selection.method = "vst", nfeatures = 2000, slot = "counts"
  )

  # Make clusters factors
  combined@meta.data$Clusters <- factor(
    combined@meta.data$Clusters,
    levels = c(paste0("C", seq_along(unique(combined@meta.data$Clusters))))
  )

  # Add spatial coordinates as embeddings
  # Remove cells that remain in coordinates but not in combined (had NA counts)
  include <- intersect(Seurat::Cells(combined), rownames(spatial_all[[1]]))
  combined <- subset(combined, cells = include)

  # Add embeddings to combined
  for (i in seq_along(samples)) {
    spatial_all[[i]] <- spatial_all[[i]][colnames(combined), ]
    combined[[samples[i]]] <- Seurat::CreateDimReducObject(
      embeddings = spatial_all[[i]],
      key = samples[i],
      assay = Seurat::DefaultAssay(combined)
    )
  }

  return(combined)
}

plot_feature <- function(seurat_obj, feature, name) {
  # Wrapper of Seurat's SpatialFeaturePlot with specific aesthetics

  Seurat::SpatialFeaturePlot(
    object = seurat_obj,
    features = feature,
    alpha = c(0.2, 1),
    pt.size.factor = 1,
    crop = FALSE
  ) +
    ggtitle(name) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 15, hjust = 0.5),
      text = element_text(size = 10)
    )
}

plot_spatial <- function(seurat_object, name) {
  #' Wrapper of Seurat's SpatialDimPlot with specific aesthetics

  clusters <- sort(unique(seurat_object$Clusters))
  colors <- ArchR::ArchRPalettes$stallion2[seq_len(length(clusters))]
  names(colors) <- clusters
  Seurat::SpatialDimPlot(
    seurat_object,
    group.by = "Clusters",
    pt.size.factor = 1,
    image.alpha = 0,
    crop = FALSE,
    cols = colors,
    stroke = 0
  ) +
    ggtitle(name) +
    theme(
      plot.title = element_text(size = 15),
      text = element_text(size = 10),
      legend.position = "bottom"
    )
}

plot_geneset <- function(seurat_obj, marker_genes, name, title) {
  #' Return a Seurat SpatialFeaturePlot of the average expression score for a
  #' set of genes.

  seurat_obj <- Seurat::AddModuleScore(
    object = seurat_obj,
    features = marker_genes,
    name = name
  )
  plot <- Seurat::SpatialFeaturePlot(
    object = seurat_obj,
    pt = 1,
    features = paste0(name, 1),
    crop = FALSE
  ) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))

  return(plot)
}

rename_cells <- function(seurat_list) {
  #' For each SeuratObj is a array of SeuratObjs, rename the cells from barcode
  #' to run_id#barcode-1; ensures cell names are unique when combining
  #' SeuratObjs

  all <-  list()
  for (i in seq_along(seurat_list)) {
    all[[i]] <- seurat_list[[i]]
    all[[i]] <- Seurat::RenameCells(
      all[[i]],
      new.names = paste0(
        unique(all[[i]]@meta.data$Sample),
        "#",
        colnames(all[[i]]),
        "-1"
      )
    )
  }
  return(all)
}

# Plot clusters ontop of spatial coordinates -----
save_spatial_cluster_plots <- function(seurat_objs) {
  #' Save plot of spatially arranged tixels colored by cluster identiy to .pdf.

  spatial_cluster_plots <- list()
  for (i in seq_along(seurat_objs)) {

    run_id <- unique(seurat_objs[[i]]$Sample)
    plot <- plot_spatial(seurat_objs[[i]], run_id)
    spatial_cluster_plots[[i]] <- plot
  }

  # Devide plots four per page -----
  spatial_lists <- split(
    spatial_cluster_plots, ceiling(seq_along(spatial_cluster_plots) / 4)
  )

  pdf("spatial_plots.pdf")
  for (i in seq_along(spatial_lists)) {
    print(Seurat::CombinePlots(spatial_lists[[i]], legend = "bottom"))
  }
  dev.off()
}

save_qc_plots <- function(seurat_objs, metrics) {
  #' Create QC plots with plot_feature

  all_qc_plots <- list()
  for (i in seq_along(metrics)) {

    spatial_qc_plots <- list()

    for (j in seq_along(seurat_objs)) {

      run_id <- unique(seurat_objs[[j]]$Sample)
      plot <- plot_feature(seurat_objs[[j]], metrics[i], run_id)
      spatial_qc_plots[[j]] <- plot
    }
    all_qc_plots[[i]] <- spatial_qc_plots
  }

  pdf("qc_plots.pdf")
  for (i in seq_along(metrics)) {
    lists <- split(all_qc_plots[[i]], ceiling(seq_along(all_qc_plots[[i]]) / 6))
    for (list in lists) {
      gridExtra::grid.arrange(grobs = list, ncol = 2)
    }
  }
  dev.off()
}

sctheme <- function(base_size = 24, xy_val = TRUE, x_ang = 0, x_jus_h = 0.5) {
  oup_theme <- theme(
    text = element_text(size = base_size, family = "Helvetica"),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black", size = base_size / 20),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = base_size),
    axis.text.x = element_text(angle = x_ang, hjust = x_jus_h),
    legend.position = "bottom",
    legend.key = element_rect(colour = NA, fill = NA),
  )
  if (!xy_val) {
    oup_theme <- oup_theme + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank()
    )
  }
  return(oup_theme)
}

scvolcano <- function(
  inpMarkers, condition1, condition2, feature = "All", fc_col = "Log2FC"
) {

  # Subset by cluster
  ggData <- inpMarkers[inpMarkers$cluster == feature, ]
  minfdr <- 0.09
  minfdr1 <- 10^-(1 / 6 * (-log10(min(ggData$p_val_adj))))

  # Add Significance column
  ggData$Significance <- ifelse(
    ggData$p_val_adj < minfdr,
    ifelse(
      ggData[[fc_col]] > 0.0,
      condition1,
      condition2
    ),
    "Not siginficant"
  )

  ggData$Significance <- factor(
    ggData$Significance,
    levels = c(condition1, condition2, "Not siginficant")
  )

  # Avoid log10(0)
  ggData[ggData$p_val_adj < 1e-300, "p_val_adj"] <- 1e-300
  ggData$log10fdr <- -log10(ggData$p_val_adj)

  # Make volcano plot
  ggOut <- ggplot(ggData, aes(x = .data[[fc_col]], y = log10fdr)) +
    geom_point() +
    sctheme() +
    ylab("-log10(FDR)") +
    geom_point(aes(color = Significance)) +
    scale_color_manual(values = c("#F8766D", "#619CFF", "gray")) +
    geom_text_repel(
      data = subset(ggData, p_val_adj < minfdr1), aes(label = gene)
    ) +
    ggtitle(paste("Markers:", feature)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 18)
    )

  return(ggOut)
}

seurat_to_h5ad <- function(seuratobj, spatial, prefix) {

  if (spatial) {
    # Reformat spatial info for squidpy compatibility
    seuratobj <- convert_spatial(seuratobj)  # From utils.R
  }

  SaveH5Seurat(seuratobj, filename = "temp.h5Seurat", overwrite = TRUE)
  Convert("temp.h5Seurat", dest = "h5ad", overwrite = TRUE)
  file.rename("temp.h5ad", paste0(prefix, "_converted.h5ad"))

  message("Seurat object converted to h5ad file.")
}