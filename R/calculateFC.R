calculateFC <- function(seurat_obj, cell_type_column = NULL, group_column, ident.1, ident.2, logfc.threshold = 0, min.pct = 0, test.use = "wilcox") {
  
  # Check if cell_type_column is provided
  if (!is.null(cell_type_column)) {
    # Split the Seurat object by the cell type
    cell_types <- unique(seurat_obj@meta.data[[cell_type_column]])
    results_list <- list()
    
    # Perform differential expression for each cell type
    for (cell_type in cell_types) {
      cat("Processing cell type:", cell_type, "\n")
      
      # Subset the Seurat object to include only the current cell type
      seurat_subset <- subset(seurat_obj, subset = seurat_obj@meta.data[[cell_type_column]] == cell_type)
      
      # Perform differential expression analysis
      de_results <- FindMarkers(
        object = seurat_subset,
        ident.1 = ident.1,
        ident.2 = ident.2,
        group.by = group_column,
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        test.use = test.use
      )
      
      # Store the results
      results_list[[cell_type]] <- de_results
    }
    
    # Return the list of results for each cell type
    return(results_list)
    
  } else {
    # No cell type specified, perform global differential expression analysis
    cat("Performing global differential expression analysis across all cells.\n")
    
    # Perform differential expression analysis
    de_results <- FindMarkers(
      object = seurat_obj,
      ident.1 = ident.1,
      ident.2 = ident.2,
      group.by = group_column,
      logfc.threshold = logfc.threshold,
      min.pct = min.pct,
      test.use = test.use
    )
    
    # Return the global differential expression results
    return(de_results)
  }
}
