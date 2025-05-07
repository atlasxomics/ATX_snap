#' Generic rountines

library("dplyr")
library("Seurat")


add_def_treatments <- function(def, treatment, feature, combined) {
  #' Edit ShinyApp scXdef.rds for Conditions.

  for (i in seq_along(treatment)) {
    def[[paste0("meta", (2 + i))]] <- treatment[i]

    def[[treatment[i]]] <- read.csv(
      paste0("req_", feature, "2_", i, ".csv")
    )$x

    def[[paste0(treatment[i], "_1")]] <- sort(
      unique(combined@meta.data[[treatment[i]]])
    )[1]

    def[[paste0(treatment[i], "_2")]] <- sort(
      unique(combined@meta.data[[treatment[i]]])
    )[2]
  }

  return(def)
}

copy_ui_server_files <- function(folder, ui_file, server_file) {

  file.copy(file.path(folder, ui_file), "/root/ui.R")
  file.copy(file.path(folder, server_file), "/root/server.R")
}

convert_to_v3 <- function(seurat_obj_v5) {

  if (!inherits(seurat_obj_v5, "Seurat")) {
    stop("Input is not a valid Seurat object.")
  }

  for (assay_name in Assays(seurat_obj_v5)) {
    assay_obj <- seurat_obj_v5[[assay_name]]
    if (inherits(assay_obj, "Assay5")) {

      counts <- assay_obj@layers[["counts"]]
      colnames(counts) <- colnames(seurat_obj_v5)
      rownames(counts) <- rownames(seurat_obj_v5)

      data <- assay_obj@layers[["data"]]
      colnames(data) <- colnames(seurat_obj_v5)
      rownames(data) <- rownames(seurat_obj_v5)

      new_assay <- CreateAssayObject(counts = counts)
      if (!is.null(data)) new_assay@data <- data

      # Add the converted assay back to the Seurat object
      seurat_obj_v5[[assay_name]] <- new_assay
    }
  }

  seurat_obj_v5@tools <- list()
  seurat_obj_v5@commands <- list()
  class(seurat_obj_v5) <- "Seurat"

  return(seurat_obj_v5)
}

convert_spatial <- function(seurat_obj) {

  # Try to extract spatial coordinates from the 'images' slot
  if (length(seurat_obj@images) > 0) {
    message("Attempting spatial coords from image slot...")

    # Check if the image is Visium format
    if (inherits(seurat_obj@images[[1]], "VisiumV1")) {
      # Extract coordinates
      coords <- data.frame(seurat_obj@images[[1]]@coordinates[, c(5, 4)])
      coords[[2]] <- -coords[[2]]
      colnames(coords) <- c("spatial_1", "spatial_2")

      # Create and add the spatial reduction
      spatial_embedding <- Seurat::CreateDimReducObject(
        embeddings = as.matrix(coords),
        key = "spatial",
        assay = Seurat::DefaultAssay(seurat_obj)
      )
      seurat_obj@reductions[["spatial"]] <- spatial_embedding

      message("Successfully created spatial coordinates from image slot.")
      return(seurat_obj)
    } else {
      message("Image slots must be class VisiumV1; trying alternative method.")
      # Continue to next method (fall through)
    }
  }

  # Try to extract spatial coordinates from reductions if Sample metadata exists
  if ("Sample" %in% colnames(seurat_obj@meta.data)) {

    message("Attempting spatial coords from reductions...")

    samples <- unique(seurat_obj@meta.data[["Sample"]])
    reduction_names <- names(seurat_obj@reductions)

    # Check if all samples have corresponding reductions
    if (all(samples %in% reduction_names)) {
      # Extract and combine spatial coordinates for each sample
      reductions <- list()
      for (sample in samples) {
        mat <- seurat_obj@reductions[[sample]]@cell.embeddings
        sample_cells <- grepl(paste0("^", sample), rownames(mat))
        if (sum(sample_cells) > 0) {
          sample_mat <- mat[sample_cells, ]
          reductions[[sample]] <- sample_mat
        } else {
          message(paste0("No matching cells found for sample: ", sample))
        }
      }

      # Only proceed if we have data to work with
      if (length(reductions) > 0) {
        spatial <- do.call(rbind, reductions)
        spatial_embedding <- Seurat::CreateDimReducObject(
          embeddings = spatial,
          key = "spatial",
          assay = Seurat::DefaultAssay(seurat_obj)
        )
        seurat_obj@reductions[["spatial"]] <- spatial_embedding
        message("Successfully reformatted spatial SeuratObject coordinates.")
        return(seurat_obj)
      }
    } else {
      message("Mismatch between sample names and
              embeddings; skipping conversion.")
      # Continue to final return
    }
  }

  # If we reach here, we couldn't find spatial coordinates
  message("Unable to find spatial coordinates in SeuratObject;
          aborting spatial conversion")
  return(seurat_obj)
}

determine_files_to_copy <- function(n_samples, n_cond, proj) {

  if (n_samples == 1) {
    return(list(ui_file = "ui_3.R", server_file = "server_3.R"))

  } else if (n_cond <= 1) {
    return(list(ui_file = "ui_2.R", server_file = "server_2.R"))

  } else if (n_cond > 1 && (
    length(unique(proj$condition_1)) > 2 || length(unique(proj$condition_2)) > 2
  )) {
    return(list(ui_file = "ui_4.R", server_file = "server_4.R"))

  } else {
    return(list(ui_file = "ui.R", server_file = "server.R"))
  }
}

edit_shiny_conf <- function(prefix, treatment) {
  #' Set default (top) values in ShinyApp/scXconf.rds files.

  conf <- readRDS(paste0("/root/shinyApp/", prefix, "conf.rds"))

  fav <- which(conf$ID %in% c("Clusters", "Sample", "SampleName", treatment))
  rest <- which(!conf$ID %in% c("Clusters", "Sample", "SampleName", treatment))
  conf <- conf[c(fav, rest), ]

  return(conf)
}

edit_shiny_def <- function(
  prefix, samples, combined, req_feats1, req_feats3, xlim, ylim
) {
  #' Edit default values in ShinyApp/scXdef.rds files.

  def <- readRDS(paste0("/root/shinyApp/", prefix, "def.rds"))

  # Add coordinate limits --
  def$limits <- list()
  for (i in seq_along(samples)) {
    def[["limits"]][[samples[i]]] <- c(
      min(xlim[[i]]), max(xlim[[i]]), min(ylim[[i]]), max(ylim[[i]])
    )
  }

  # Set defaults for meta --
  def$meta1 <- "Clusters"
  def$meta2 <- "Sample"
  def$meta3 <- "SampleName"

  # Add default feats for Clusters, Sample --
  def$Clusters <- req_feats1
  def$Sample <- req_feats3

  # Set reductions values --
  def$dimred[3] <- paste0(names(combined[3]@reductions)[1], "_1")
  def$dimred[4] <- paste0(names(combined[3]@reductions)[1], "_2")
  def$dimred[5] <- paste0(names(combined[3]@reductions)[2], "_1")
  def$dimred[6] <- paste0(names(combined[3]@reductions)[2], "_2")

  return(def)
}

find_func <- function(tempdir, pattern) {
  list.files(
    path = tempdir,     #' Replace tempdir with the directory you want, followed
    pattern = pattern,  #' by 0 or more characters, then ".csv", and then
    full.names = TRUE,  #' nothing else ($) include the directory in the result.
    recursive = TRUE
  )
}

find_sample_names <- function(seurat_lst) {
  #' Extract list of sample names from list of SeuratObjs.
  sapply(
    seq_along(seurat_lst), function(i) {
      unique(seurat_lst[[i]]@meta.data$Sample)
    }
  )
}

extract_nth_ele <- function(lst, n = 1) {
  #' Extract nth element for each item in an array; return as a vector.
  sapply(lst, `[`, n)
}

metrics_to_csv <- function(metric_list) {
  #' Write inpute metrics to csv
  names(metric_list) <- c(
    "project_name",
    "genome",
    "tile_size",
    "minimum_tss",
    "minimum_fragments",
    "lsi_iterations",
    "lsi_resolution",
    "lsi_varFeatures",
    "clustering_resolution",
    "umap_minimum_distance",
    "number_threads",
    "min_cells_cluster",
    "max_clusters"
  )
  write.csv(metric_list, file = "exe_metadata.csv", row.names = FALSE)
}