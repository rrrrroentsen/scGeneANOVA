#' scGeneANOVA
#'
#' This function performs gene expression analysis and ANOVA on scRNA-seq data,
#' with additional functionality to calculate fold change (FC) between specified groups.
#'
#' @param seurat_obj A Seurat object containing the scRNA-seq data.
#' @param gene_list A list of genes to analyze. If NULL, all genes are used.
#' @param cell_type_column The column name in the metadata that contains cell type information. If NULL, all cells are treated as one group.
#' @param group_column The column name in the metadata that contains group information.
#' @param sample_column The column name in the metadata that contains sample information.
#' @param chunk_size The number of genes to process in each chunk. Default is 100.
#'
#' @return A dataframe containing ANOVA, Tukey's test results, and fold change (FC) calculations.
#' @examples
#' \dontrun{
#' seurat_obj <- readRDS("path/to/your/seurat_obj.rds")
#' results <- scGeneANOVA(seurat_obj, group_column = "Patient_ident", sample_column = "orig.ident")
#' }
#' @export
scGeneANOVA <- function(seurat_obj, gene_list = NULL, cell_type_column = NULL, group_column, sample_column, chunk_size = 100) {
  
  # If gene_list is NULL, use all genes in the Seurat object
  if (is.null(gene_list)) {
    gene_list <- rownames(seurat_obj@assays$RNA@data)
  }

  # If cell_type_column is NULL, treat all cells as a single group
  if (is.null(cell_type_column)) {
    seurat_obj@meta.data$Default_Cell_Type <- "All_Cells"
    cell_type_column <- "Default_Cell_Type"
  }

  # Function to calculate average gene expression for each sample and cell type
  analyze_gene_expression <- function(seurat_obj, genes_chunk, cell_type_column, sample_column) {
    # Fetch data for specified genes, cell types, and samples
    expression_data <- FetchData(seurat_obj, vars = c(genes_chunk, cell_type_column, sample_column))
    avg_exp_cols <- paste0("avg.exp_", genes_chunk)
    
    # Calculate mean expression for each gene across samples and cell types
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

  # Process genes in chunks and calculate average expression
  final_result <- bind_rows(lapply(seq(1, length(gene_list), by = chunk_size), function(start_idx) {
    pb$tick()
    end_idx <- min(start_idx + chunk_size - 1, length(gene_list))
    genes_chunk <- gene_list[start_idx:end_idx]
    analyze_gene_expression(seurat_obj, genes_chunk, cell_type_column, sample_column)
  }))

  # Add patient information to the result dataframe
  patient_column <- seurat_obj@meta.data[[group_column]][match(final_result[[sample_column]], seurat_obj@meta.data[[sample_column]])]
  final_result <- cbind(final_result, Patient = patient_column)

  # Function to perform ANOVA and Tukey's test for a single gene and cell type
  anova_single_gene_celltype <- function(data, gene, cell_type, cell_type_column) {
    anova_data <- filter(data, Gene == gene, !!sym(cell_type_column) == cell_type)

    # Perform ANOVA and Tukey's test if there are at least 2 groups with data
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

  # Get unique cell types for analysis
  cell_types <- unique(final_result[[cell_type_column]])

  # Perform ANOVA and Tukey's test for all genes and cell types
  anova_tukey_results_df <- anova_all_genes(final_result, gene_list, cell_types, cell_type_column)

  # Get unique groups for pairwise comparisons
  unique_groups <- unique(seurat_obj@meta.data[[group_column]])
  group_combinations <- combn(unique_groups, 2, simplify = FALSE)

  # Initialize list to hold FC results
  fc_results_list <- list()

  # Calculate FC for each pairwise group combination
  total_combinations <- length(group_combinations)
  pb_fc <- progress_bar$new(
    format = "  Calculating FC [:bar] :percent eta: :eta",
    total = total_combinations, clear = FALSE, width = 60
  )

  for (comb in group_combinations) {
    ident1 <- comb[1]
    ident2 <- comb[2]
    fc_result <- calculateFC(seurat_obj, 
                             cell_type_column = cell_type_column, 
                             group_column = group_column, 
                             ident.1 = ident1, 
                             ident.2 = ident2, 
                             cell_type = NULL, 
                             features = gene_list, 
                             slot = "data", 
                             pseudocount.use = 1, 
                             base = 2)
    pb_fc$tick()

    # Check if FC result is not empty
    if (!is.null(fc_result) && nrow(fc_result) > 0) {
      fc_result$Gene <- rownames(fc_result)

      # Standardize group comparison naming convention
      fc_result$Comparison <- gsub(" vs ", "-", fc_result$Group)

      # Correct reversed comparisons if necessary
      reversed_comparisons <- sapply(fc_result$Comparison, function(x) {
        parts <- strsplit(x, "-")[[1]]
        if (paste0(parts, collapse = "-") %in% anova_tukey_results_df$Comparison) {
          return(FALSE)
        } else if (paste0(rev(parts), collapse = "-") %in% anova_tukey_results_df$Comparison) {
          return(TRUE)
        } else {
          return(NA)
        }
      })

      # Adjust reversed comparisons
      fc_result$Comparison[reversed_comparisons] <- sapply(fc_result$Comparison[reversed_comparisons], function(x) {
        parts <- strsplit(x, "-")[[1]]
        paste0(rev(parts), collapse = "-")
      })
      fc_result$avg_logFC[reversed_comparisons] <- -fc_result$avg_logFC[reversed_comparisons]
      temp_pct <- fc_result$pct.1[reversed_comparisons]
      fc_result$pct.1[reversed_comparisons] <- fc_result$pct.2[reversed_comparisons]
      fc_result$pct.2[reversed_comparisons] <- temp_pct

      fc_result$Comparison <- as.character(fc_result$Comparison)

      fc_results_list[[paste(ident1, ident2, sep = ":")]] <- fc_result
    }
  }

  # Combine FC results with ANOVA results
  all_fc_results <- bind_rows(fc_results_list)

  # Merge FC results with ANOVA and Tukey's test results
  if (nrow(all_fc_results) > 0) {
    anova_tukey_results_df$Comparison <- as.character(anova_tukey_results_df$Comparison)
    all_fc_results$Comparison <- as.character(all_fc_results$Comparison)

    final_merged_output <- merge(anova_tukey_results_df, all_fc_results, by = c("Gene", "Cell_Type", "Comparison"), all = TRUE)

    # Remove the Group column from the final output
    final_merged_output <- final_merged_output %>% select(-Group)
  } else {
    final_merged_output <- anova_tukey_results_df
  }

  return(final_merged_output)
}
