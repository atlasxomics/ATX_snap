

# Function to process gene matrix in chunks
processArchRProjectInChunks <- function(
  proj,
  chunk_size = 5000,  # Adjust this value based on your system's memory
  output_dir = getwd(),
  random_seed = 42,   # For reproducibility in random sampling
  stratify = TRUE     # Stratify by sample/run to prevent bias
) {
  # Extract metadata first (smaller operation)
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
  
  # Get all cell names to chunk
  all_cells <- rownames(proj@cellColData)
  n_cells <- length(all_cells)
  
  # Set seed for reproducibility
  set.seed(random_seed)
  
  # Create chunks in a way that ensures even representation
  if (stratify && "Sample" %in% colnames(proj@cellColData)) {
    print("Stratifying chunks by sample to prevent bias...")
    
    # Get sample information - extract as simple character vector
    # Fix: Convert S4 object to a simple vector
    sample_info <- as.character(proj@cellColData$Sample)
    names(sample_info) <- all_cells
    
    unique_samples <- unique(sample_info)
    n_samples <- length(unique_samples)
    
    print(paste0("Found ", n_samples, " unique samples/runs for stratified chunking"))
    
    # Determine cells per sample
    cells_per_sample <- list()
    for (sample in unique_samples) {
      # Get cell names (indexes) for this sample
      sample_cell_names <- names(sample_info)[sample_info == sample]
      cells_per_sample[[sample]] <- sample_cell_names
      print(paste0("Sample ", sample, " has ", length(cells_per_sample[[sample]]), " cells"))
    }
    
    # Random shuffle cells within each sample to prevent any bias
    for (sample in unique_samples) {
      cells_per_sample[[sample]] <- sample(cells_per_sample[[sample]], size = length(cells_per_sample[[sample]]))
    }
    
    # Calculate number of chunks needed
    n_chunks <- ceiling(n_cells / chunk_size)
    print(paste0("Processing ", n_cells, " cells in ", n_chunks, " chunks of up to ", chunk_size, " cells each"))
    
    # Create empty list to hold chunks
    chunks <- vector("list", n_chunks)
    for (i in 1:n_chunks) {
      chunks[[i]] <- character(0)  # Initialize with empty character vector
    }
    
    # Distribute cells from each sample evenly across chunks
    for (sample in unique_samples) {
      sample_cells <- cells_per_sample[[sample]]
      n_sample_cells <- length(sample_cells)
      
      # Calculate how many cells from this sample should go in each chunk
      cells_per_chunk <- ceiling(n_sample_cells / n_chunks)
      
      for (i in 1:n_chunks) {
        start_idx <- (i - 1) * cells_per_chunk + 1
        end_idx <- min(i * cells_per_chunk, n_sample_cells)
        
        if (start_idx <= end_idx) {
          # Add cells from this sample to this chunk
          chunk_cells <- sample_cells[start_idx:end_idx]
          chunks[[i]] <- c(chunks[[i]], chunk_cells)
        }
      }
    }
  } else {
    # Simple random chunking as an alternative
    print("Using random chunking to prevent bias...")
    
    # Shuffle all cells randomly
    shuffled_cells <- sample(all_cells, size = n_cells)
    
    # Calculate number of chunks
    n_chunks <- ceiling(n_cells / chunk_size)
    print(paste0("Processing ", n_cells, " cells in ", n_chunks, " chunks of up to ", chunk_size, " cells each"))
    
    # Create chunks
    chunks <- list()
    for (i in 1:n_chunks) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, n_cells)
      
      chunks[[i]] <- shuffled_cells[start_idx:end_idx]
    }
  }
  
  # Verify chunk sizes
  for (i in 1:length(chunks)) {
    print(paste0("Chunk ", i, " contains ", length(chunks[[i]]), " cells"))
  }
  
  # Initialize an empty list to store gene matrices
  chunks_list <- list()
  gene_row_names <- NULL
  
  # Get imputation weights for the full project
  # This ensures consistent imputation across chunks
  print("Calculating imputation weights from full project for consistent imputation...")
  full_impute_weights <- getImputeWeights(proj)
  
  # Loop through chunks
  for (i in 1:n_chunks) {
    print(paste0("Processing chunk ", i, " of ", n_chunks))
    
    # Get cells for current chunk
    chunk_cells <- chunks[[i]]
    
    # Create a temporary sub-project with these cells
    chunk_proj <- subsetArchRProject(
      ArchRProj = proj,
      cells = chunk_cells,
      outputDirectory = paste0(output_dir, "/temp_chunk_", i),
      force = TRUE
    )
    
    # Extract gene matrix for this chunk
    chunk_matrix <- getMatrixFromProject(
      ArchRProj = chunk_proj,
      useMatrix = "GeneScoreMatrix"
    )
    
    # Save the row names from the first chunk (gene names should be consistent across chunks)
    if (i == 1) {
      gene_row_names <- chunk_matrix@elementMetadata$name
    }
    
    # Impute the chunk matrix using the global imputation weights for consistency
    print("Imputing chunk with global imputation weights for consistency across chunks...")
    chunk_imputed <- imputeMatrix(
      mat = assay(chunk_matrix),
      imputeWeights = full_impute_weights
    )
    rownames(chunk_imputed) <- gene_row_names
    
    # Store the chunk
    chunks_list[[i]] <- chunk_imputed
    
    # Clean up to free memory
    rm(chunk_proj, chunk_matrix, chunk_imputed)
    gc(verbose = FALSE, full = TRUE)
  }
  
  # Return the chunks and metadata
  return(list(chunks_list = chunks_list, metadata = metadata, gene_row_names = gene_row_names))
}


# Modified processArchRProjectInChunks function with added debugging
processArchRProjectInChunks <- function(
  proj,
  chunk_size = 5000,  # Adjust this value based on your system's memory
  output_dir = getwd(),
  random_seed = 42,   # For reproducibility in random sampling
  stratify = TRUE     # Stratify by sample/run to prevent bias
) {
  # Extract metadata first (smaller operation)
  print("Extracting metadata...")
  metadata <- getCellColData(ArchRProj = proj)
  
  # Set metadata rownames to barcodes
  print("Formatting metadata rownames...")
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
  
  # Get all cell names to chunk
  all_cells <- rownames(proj@cellColData)
  n_cells <- length(all_cells)
  
  print(paste0("Total cells to process: ", n_cells))
  
  # Set seed for reproducibility
  set.seed(random_seed)
  
  # Create chunks in a way that ensures even representation
  if (stratify && "Sample" %in% colnames(proj@cellColData)) {
    print("Stratifying chunks by sample to prevent bias...")
    
    # Get sample information - extract as simple character vector
    # Fix: Convert S4 object to a simple vector
    sample_info <- as.character(proj@cellColData$Sample)
    names(sample_info) <- all_cells
    
    unique_samples <- unique(sample_info)
    n_samples <- length(unique_samples)
    
    print(paste0("Found ", n_samples, " unique samples/runs for stratified chunking"))
    
    # Determine cells per sample
    cells_per_sample <- list()
    for (sample in unique_samples) {
      # Get cell names (indexes) for this sample
      sample_cell_names <- names(sample_info)[sample_info == sample]
      cells_per_sample[[sample]] <- sample_cell_names
      print(paste0("Sample ", sample, " has ", length(cells_per_sample[[sample]]), " cells"))
    }
    
    # Random shuffle cells within each sample to prevent any bias
    for (sample in unique_samples) {
      cells_per_sample[[sample]] <- sample(cells_per_sample[[sample]], size = length(cells_per_sample[[sample]]))
    }
    
    # Calculate number of chunks needed
    n_chunks <- ceiling(n_cells / chunk_size)
    print(paste0("Processing ", n_cells, " cells in ", n_chunks, " chunks of up to ", chunk_size, " cells each"))
    
    # Create empty list to hold chunks
    chunks <- vector("list", n_chunks)
    for (i in 1:n_chunks) {
      chunks[[i]] <- character(0)  # Initialize with empty character vector
    }
    
    # Distribute cells from each sample evenly across chunks
    for (sample in unique_samples) {
      sample_cells <- cells_per_sample[[sample]]
      n_sample_cells <- length(sample_cells)
      
      # Calculate how many cells from this sample should go in each chunk
      cells_per_chunk <- ceiling(n_sample_cells / n_chunks)
      
      for (i in 1:n_chunks) {
        start_idx <- (i - 1) * cells_per_chunk + 1
        end_idx <- min(i * cells_per_chunk, n_sample_cells)
        
        if (start_idx <= end_idx) {
          # Add cells from this sample to this chunk
          chunk_cells <- sample_cells[start_idx:end_idx]
          chunks[[i]] <- c(chunks[[i]], chunk_cells)
        }
      }
    }
  } else {
    # Simple random chunking as an alternative
    print("Using random chunking to prevent bias...")
    
    # Shuffle all cells randomly
    shuffled_cells <- sample(all_cells, size = n_cells)
    
    # Calculate number of chunks
    n_chunks <- ceiling(n_cells / chunk_size)
    print(paste0("Processing ", n_cells, " cells in ", n_chunks, " chunks of up to ", chunk_size, " cells each"))
    
    # Create chunks
    chunks <- list()
    for (i in 1:n_chunks) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, n_cells)
      
      chunks[[i]] <- shuffled_cells[start_idx:end_idx]
    }
  }
  
  # Verify chunk sizes
  for (i in 1:length(chunks)) {
    print(paste0("Chunk ", i, " contains ", length(chunks[[i]]), " cells"))
  }
  
  # Initialize an empty list to store gene matrices
  print("Initializing chunks_list to store gene matrices...")
  chunks_list <- list()
  gene_row_names <- NULL
  
  # Get imputation weights for the full project
  # This ensures consistent imputation across chunks
  print("Calculating imputation weights from full project for consistent imputation...")
  tryCatch({
    full_impute_weights <- getImputeWeights(proj)
    print("Successfully calculated imputation weights")
  }, error = function(e) {
    print(paste0("Error calculating imputation weights: ", e$message))
    return(NULL)
  })
  
  # Loop through chunks
  print("Beginning processing of individual chunks...")
  for (i in 1:n_chunks) {
    print(paste0("Processing chunk ", i, " of ", n_chunks))
    
    # Get cells for current chunk
    chunk_cells <- chunks[[i]]
    
    # Create a temporary sub-project with these cells
    print(paste0("Creating subset project for chunk ", i, "..."))
    tryCatch({
      chunk_proj <- subsetArchRProject(
        ArchRProj = proj,
        cells = chunk_cells,
        outputDirectory = paste0(output_dir, "/temp_chunk_", i),
        force = TRUE
      )
      print(paste0("Successfully created subset project for chunk ", i))
    }, error = function(e) {
      print(paste0("Error creating subset project for chunk ", i, ": ", e$message))
      return(NULL)
    })
    
    # Extract gene matrix for this chunk
    print(paste0("Extracting gene matrix for chunk ", i, "..."))
    tryCatch({
      chunk_matrix <- getMatrixFromProject(
        ArchRProj = chunk_proj,
        useMatrix = "GeneScoreMatrix"
      )
      print(paste0("Successfully extracted gene matrix for chunk ", i))
    }, error = function(e) {
      print(paste0("Error extracting gene matrix for chunk ", i, ": ", e$message))
      return(NULL)
    })
    
    # Save the row names from the first chunk (gene names should be consistent across chunks)
    if (i == 1) {
      gene_row_names <- chunk_matrix@elementMetadata$name
      print(paste0("Saved ", length(gene_row_names), " gene row names from first chunk"))
    }
    
    # Impute the chunk matrix using the global imputation weights for consistency
    print(paste0("Imputing chunk ", i, " with global imputation weights..."))
    tryCatch({
      chunk_imputed <- imputeMatrix(
        mat = assay(chunk_matrix),
        imputeWeights = full_impute_weights
      )
      rownames(chunk_imputed) <- gene_row_names
      print(paste0("Successfully imputed chunk ", i))
    }, error = function(e) {
      print(paste0("Error imputing chunk ", i, ": ", e$message))
      return(NULL)
    })
    
    # Store the chunk
    print(paste0("Storing chunk ", i, " in chunks_list"))
    chunks_list[[i]] <- chunk_imputed
    
    # Clean up to free memory
    rm(chunk_proj, chunk_matrix, chunk_imputed)
    gc(verbose = FALSE, full = TRUE)
  }
  
  print("Final chunks_list length: ", length(chunks_list))
  print("Gene row names length: ", length(gene_row_names))
  print("Metadata dimensions: ", dim(metadata)[1], "x", dim(metadata)[2])
  
  # Return the chunks and metadata
  print("Returning results...")
  results <- list(
    chunks_list = chunks_list, 
    metadata = metadata, 
    gene_row_names = gene_row_names
  )
  
  return(results)
}