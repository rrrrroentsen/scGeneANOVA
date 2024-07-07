#' scGeneANOVA
#'
#' This function performs gene expression analysis and ANOVA on scRNA-seq data.
#'
#' @param seurat_obj A Seurat object containing the scRNA-seq data.
#' @param gene_list A list of genes to analyze. If NULL, all genes are used.
#' @param cell_type_column The column name in the metadata that contains cell type information. If NULL, all cells are treated as one group.
#' @param group_column The column name in the metadata that contains group information.
#' @param sample_column The column name in the metadata that contains sample information.
#' @param chunk_size The number of genes to process in each chunk. Default is 100.
#'
#' @return A dataframe containing ANOVA and Tukey's test results.
#' @examples
#' \dontrun{
#' seurat_obj <- readRDS("path/to/your/seurat_obj.rds")
#' results <- scGeneANOVA(seurat_obj, group_column = "Patient_ident", sample_column = "orig.ident")
#' }
#' @export
scGeneANOVA <- function(seurat_obj, gene_list = NULL, cell_type_column = NULL, group_column, sample_column, chunk_size = 100) {
  # Function implementation here
}



scGeneANOVA <- function(seurat_obj, gene_list = NULL, cell_type_column = NULL, group_column, sample_column, chunk_size = 100) {
  # If gene_list is NULL, use all genes
  if (is.null(gene_list)) {
    gene_list <- rownames(seurat_obj@assays$RNA@data)
  }

  # If cell_type_column is NULL, analyze the entire dataset without cell type differentiation
  if (is.null(cell_type_column)) {
    seurat_obj@meta.data$Default_Cell_Type <- "All_Cells"
    cell_type_column <- "Default_Cell_Type"
  }

  # Function to calculate average expression and expression proportion
  analyze_gene_expression <- function(seurat_obj, genes_chunk, cell_type_column, sample_column) {
    expression_data <- FetchData(seurat_obj, vars = c(genes_chunk, cell_type_column, sample_column))
    avg_exp_cols <- paste0("avg.exp_", genes_chunk)
    expression_data %>%
      group_by(!!sym(sample_column), !!sym(cell_type_column)) %>%
      summarize(across(all_of(genes_chunk), ~ mean(expm1(.), na.rm = TRUE), .names = "avg.exp_{.col}"), .groups = 'drop') %>%
      pivot_longer(cols = all_of(avg_exp_cols), names_to = "Gene", values_to = "avg.exp") %>%
      mutate(Gene = sub("avg.exp_", "", Gene))
  }

  # Initialize result dataframe and progress bar
  total_chunks <- ceiling(length(gene_list) / chunk_size)
  pb <- progress_bar$new(
    format = "  Analyzing genes [:bar] :percent eta: :eta",
    total = total_chunks, clear = FALSE, width = 60
  )

  final_result <- bind_rows(lapply(seq(1, length(gene_list), by = chunk_size), function(start_idx) {
    pb$tick()
    end_idx <- min(start_idx + chunk_size - 1, length(gene_list))
    genes_chunk <- gene_list[start_idx:end_idx]
    analyze_gene_expression(seurat_obj, genes_chunk, cell_type_column, sample_column)
  }))

  # Ensure Patient column is correct size
  patient_column <- seurat_obj@meta.data[[group_column]][match(final_result[[sample_column]], seurat_obj@meta.data[[sample_column]])]
  final_result <- cbind(final_result, Patient = patient_column)

  # Function to perform ANOVA and Tukey's test for single gene and cell type
  anova_single_gene_celltype <- function(data, gene, cell_type, cell_type_column) {
    anova_data <- filter(data, Gene == gene, !!sym(cell_type_column) == cell_type)

    # Check if there are at least 2 groups with more than one data point
    if (n_distinct(anova_data$Patient) < 2) return(NULL)

    tryCatch({
      aov_result <- aov(avg.exp ~ Patient, data = anova_data)
      tukey_result <- TukeyHSD(aov_result)
      list(anova = summary(aov_result), tukey = tukey_result)
    }, error = function(e) NULL)
  }

  # Perform ANOVA and Tukey's test for all genes and cell types
  anova_all_genes <- function(data, gene_list, cell_types, cell_type_column) {
    pb <- progress_bar$new(
      format = "  Performing ANOVA and Tukey's test [:bar] :percent eta: :eta",
      total = length(gene_list) * length(cell_types), clear = FALSE, width = 60
    )
    results <- lapply(gene_list, function(gene) {
      lapply(cell_types, function(cell_type) {
        pb$tick()
        result <- anova_single_gene_celltype(data, gene, cell_type, cell_type_column)
        if (is.null(result)) {
          data.frame(Comparison = NA, diff = NA, lwr = NA, upr = NA, p.adj = NA, Gene = gene, Cell_Type = cell_type, P_Value = NA)
        } else {
          tukey_df <- as.data.frame(result$tukey$Patient)
          tukey_df$Comparison <- rownames(tukey_df)
          tukey_df$Gene <- gene
          tukey_df$Cell_Type <- cell_type
          tukey_df$P_Value <- result$anova[[1]][["Pr(>F)"]][1]
          tukey_df
        }
      })
    })
    bind_rows(results)
  }

  # Get unique cell types
  cell_types <- unique(final_result[[cell_type_column]])

  # Perform ANOVA and Tukey's test for all genes and cell types
  anova_tukey_results_df <- anova_all_genes(final_result, gene_list, cell_types, cell_type_column)

  return(anova_tukey_results_df)
}

