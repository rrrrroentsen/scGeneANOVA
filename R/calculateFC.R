#' calculateFC
#'
#' This function calculates the fold change (FC) in gene expression between two groups of cells.
#'
#' @param seurat_obj A Seurat object containing the scRNA-seq data.
#' @param cell_type_column The column name in the metadata that contains cell type information.
#' @param group_column The column name in the metadata that contains group information.
#' @param ident.1 The first group to compare.
#' @param ident.2 The second group to compare.
#' @param cell_type A specific cell type to analyze. If NULL, all cell types are analyzed together.
#' @param features A vector of gene names to analyze. If NULL, all genes are analyzed.
#' @param slot The data slot to use for expression values (e.g., "data" or "scale.data").
#' @param pseudocount.use A small constant to add to expression values to avoid log transformation issues.
#' @param base The logarithm base used for fold change calculation.
#'
#' @return A dataframe containing fold change (FC) results for the specified genes and groups.
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

    # Subset the Seurat object by cell type if specified
    if (!is.null(cell_type)) {
        seurat_obj <- subset(seurat_obj, seurat_obj[[cell_type_column]] == cell_type)
        cell_type_name <- cell_type
    } else {
        cell_type_name <- "All_Cells"
    }

    # Set identities to the grouping column
    Idents(seurat_obj) <- group_column

    # Select cells from the two groups
    cells.1 <- WhichCells(seurat_obj, ident = ident.1)
    cells.2 <- WhichCells(seurat_obj, ident = ident.2)

    # Get the expression matrix from the specified slot
    data <- GetAssayData(object = seurat_obj, slot = slot)

    # Define the default mean function based on slot and pseudocount
    default.mean.fxn <- function(x) {
        return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
    }

    # Check normalization method if slot is "data"
    norm.method <- NULL
    if (slot == "data") {
        norm.command <- paste0("NormalizeData.", DefaultAssay(object = seurat_obj))
        if (norm.command %in% Command(object = seurat_obj)) {
            norm.method <- Command(
                object = seurat_obj,
                command = norm.command,
                value = "normalization.method"
            )
        }
    }

    # Select mean function based on the slot and normalization method
    mean.fxn <- switch(
        EXPR = slot,
        'data' = switch(
            EXPR = norm.method %||% '',
            'LogNormalize' = function(x) {
                return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
            },
            default.mean.fxn
        ),
        'scale.data' = rowMeans,
        default.mean.fxn
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

    # Add group identities and cell type to the results
    if (!is.null(fc.results) && nrow(fc.results) > 0) {
        fc.results$Group <- paste(ident.1, "vs", ident.2)
        fc.results$Cell_Type <- cell_type_name
        
        fc.results$Comparison <- gsub(" vs ", "-", fc.results$Group)
        fc.results$Comparison <- as.character(fc.results$Comparison)
    }

    return(fc.results)
}
