calculateFC <- function(seurat_obj, 
                        cell_type_column = NULL, 
                        group_column, 
                        ident.1, 
                        ident.2, 
                        logfc.threshold = 0, 
                        min.pct = 0) {
  
  # Function to estimate time
  estimate_time <- function(expr) {
    start_time <- Sys.time()
    expr
    end_time <- Sys.time()
    return(end_time - start_time)
  }

  # Print a message indicating the start of the analysis
  if (!is.null(cell_type_column)) {
    message("Performing differential expression analysis by cell type.")
  } else {
    message("Performing global differential expression analysis across all cells.")
  }
  
  results_list <- list()

  if (!is.null(cell_type_column)) {
    # If cell_type_column is provided, split the data by cell type and analyze each separately
    cell_types <- unique(seurat_obj@meta.data[[cell_type_column]])
    
    for (cell_type in cell_types) {
      message(paste("Analyzing cell type:", cell_type))
      
      # Subset Seurat object for the specific cell type
      cells_of_type <- seurat_obj@meta.data[[cell_type_column]] == cell_type
      seurat_subset <- subset(seurat_obj, cells = which(cells_of_type))
      
      # Get expression data and metadata
      data <- as.matrix(seurat_subset@assays$RNA@data)
      metadata <- seurat_subset@meta.data
      
      # Group indices
      group1_idx <- which(metadata[[group_column]] == ident.1)
      group2_idx <- which(metadata[[group_column]] == ident.2)
      
      # Calculate mean expression for each gene in both groups
      group1_means <- rowMeans(data[, group1_idx])
      group2_means <- rowMeans(data[, group2_idx])
      
      # Calculate log fold change
      logFC <- log2(group1_means + 1) - log2(group2_means + 1)
      
      # Calculate p-values using Wilcoxon test
      p_values <- apply(data, 1, function(gene_expr) {
        wilcox.test(gene_expr[group1_idx], gene_expr[group2_idx])$p.value
      })
      
      # Adjust p-values for multiple testing
      p_adj <- p.adjust(p_values, method = "bonferroni")
      
      # Filter results based on logfc.threshold and min.pct
      pct1 <- rowMeans(data[, group1_idx] > 0)
      pct2 <- rowMeans(data[, group2_idx] > 0)
      
      results <- data.frame(
        gene = rownames(data),
        logFC = logFC,
        p_val = p_values,
        p_val_adj = p_adj,
        pct1 = pct1,
        pct2 = pct2
      )
      
      results <- subset(results, abs(logFC) > logfc.threshold & (pct1 > min.pct | pct2 > min.pct))
      
      results_list[[cell_type]] <- results
    }
    
  } else {
    # Global analysis without separating by cell type
    
    # Get expression data and metadata
    data <- as.matrix(seurat_obj@assays$RNA@data)
    metadata <- seurat_obj@meta.data
    
    # Group indices
    group1_idx <- which(metadata[[group_column]] == ident.1)
    group2_idx <- which(metadata[[group_column]] == ident.2)
    
    # Calculate mean expression for each gene in both groups
    group1_means <- rowMeans(data[, group1_idx])
    group2_means <- rowMeans(data[, group2_idx])
    
    # Calculate log fold change
    logFC <- log2(group1_means + 1) - log2(group2_means + 1)
    
    # Calculate p-values using Wilcoxon test
    p_values <- apply(data, 1, function(gene_expr) {
      wilcox.test(gene_expr[group1_idx], gene_expr[group2_idx])$p.value
    })
    
    # Adjust p-values for multiple testing
    p_adj <- p.adjust(p_values, method = "bonferroni")
    
    # Filter results based on logfc.threshold and min.pct
    pct1 <- rowMeans(data[, group1_idx] > 0)
    pct2 <- rowMeans(data[, group2_idx] > 0)
    
    results <- data.frame(
      gene = rownames(data),
      logFC = logFC,
      p_val = p_values,
      p_val_adj = p_adj,
      pct1 = pct1,
      pct2 = pct2
    )
    
    results <- subset(results, abs(logFC) > logfc.threshold & (pct1 > min.pct | pct2 > min.pct))
    
    results_list[["global"]] <- results
  }

  return(results_list)
}
