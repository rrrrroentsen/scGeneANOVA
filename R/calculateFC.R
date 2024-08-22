#' calculateFC
#'
#' This function performs differential gene expression analysis between two specified groups within a Seurat object. It calculates log fold changes, p-values, and adjusted p-values for the genes, with an option to subset by cell type. It is a key step in the scGeneANOVA pipeline.
#'
#' @param seurat_obj A Seurat object containing the scRNA-seq data.
#' @param cell_type_column The column name in the metadata that contains cell type information. If NULL, analysis is performed on all cells.
#' @param group_column The column name in the metadata that contains group information (e.g., treatment groups).
#' @param ident.1 The identity class to compare as group 1.
#' @param ident.2 The identity class to compare as group 2.
#' @param cell_type A specific cell type to subset the data for comparison. If NULL, all cells are used.
#' @param features A vector of features (e.g., gene names) to include in the analysis. If NULL, all features are included.
#' @param slot The data slot to use ("data", "scale.data", or "counts"). Default is "data".
#' @param pseudocount.use A numeric value to add to the data to avoid division by zero. Default is 1.
#' @param base The base of the logarithm used for fold change calculations. Default is 2.
#'
#' @return A dataframe containing the log fold changes, percentage of cells expressing the gene in each group, p-values, adjusted p-values, group comparison, cell type, and gene names.
#'
#' @examples
#' \dontrun{
#' seurat_obj <- readRDS("path/to/your/seurat_obj.rds")
#' results <- calculateFC(seurat_obj, group_column = "Patient_ident", ident.1 = "INR", ident.2 = "IR")
#' }
#' @export
calculateFC <- function(seurat_obj, 
                        cell_type_column, 
                        group_column, 
                        ident.1, 
                        ident.2, 
                        cell_type = NULL, 
                        features = NULL, 
                        slot = "data", 
                        pseudocount.use = 1, 
                        base = 2) {

    # Start time for execution
    start_time <- Sys.time()

    # Subset by cell type if specified
    if (!is.null(cell_type)) {
        seurat_obj <- subset(seurat_obj, seurat_obj[[cell_type_column]] == cell_type)
        cell_type_name <- cell_type
    } else {
        cell_type_name <- "All"
    }

    # Set identities to the grouping column
    Idents(seurat_obj) <- group_column

    # Select cells from the two groups
    cells.1 <- WhichCells(seurat_obj, ident = ident.1)
    cells.2 <- WhichCells(seurat_obj, ident = ident.2)

    # Get expression matrix
    data <- GetAssayData(object = seurat_obj, slot = slot)

    # Define mean function based on slot and pseudocount
    log1pdata.mean.fxn <- function(x) {
        return(log(x = (rowSums(x = expm1(x = x)) + pseudocount.use) / NCOL(x), base = base))
    }

    scaledata.mean.fxn <- rowMeans
    counts.mean.fxn <- function(x) {
        return(log(x = (rowSums(x = x) + pseudocount.use) / NCOL(x), base = base))
    }

    # Select mean function based on the slot
    mean.fxn <- switch(
        EXPR = slot,
        'data' = log1pdata.mean.fxn,
        'scale.data' = scaledata.mean.fxn,
        'counts' = counts.mean.fxn,
        log1pdata.mean.fxn
    )

    # Calculate fold change
    fc.results <- FoldChange(
        object = data,
        cells.1 = cells.1,
        cells.2 = cells.2,
        features = features,
        mean.fxn = mean.fxn,
        pseudocount.use = pseudocount.use,
        base = base,
        fc.name = "avg_logFC"
    )

    # Add p-values and adjusted p-values
    test.results <- FindMarkers(seurat_obj, ident.1 = ident.1, ident.2 = ident.2, features = features, logfc.threshold = 0, min.pct = 0, only.pos = FALSE, slot = slot)
    
    # Merge fold change with test results
    combined.results <- merge(fc.results, test.results[, c("p_val", "p_val_adj")], by = "row.names")
    colnames(combined.results)[1] <- "Gene" # Rename the first column to "Gene"
    
    # Add columns for group identities and cell type
    combined.results$Group <- paste(ident.1, "vs", ident.2)
    combined.results$Cell_Type <- cell_type_name

    # Calculate time taken
    end_time <- Sys.time()
    time_taken <- end_time - start_time
    
    # Print the execution time
    print(paste("Execution Time:", time_taken))
    
    return(combined.results)
}
