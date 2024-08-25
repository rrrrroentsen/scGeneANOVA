run_scGeneANOVA <- function(seurat_obj, gene_list = NULL, cell_type_column = NULL, group_column, sample_column, chunk_size = 100) {
    
    # If gene_list is NULL, use all genes in the Seurat object
    if (is.null(gene_list)) {
        gene_list <- rownames(seurat_obj@assays$RNA@data)
    }
    
    # If cell_type_column is NULL, treat all cells as a single group
    if (is.null(cell_type_column)) {
        seurat_obj@meta.data$Default_Cell_Type <- "All_Cells"
        cell_type_column <- "Default_Cell_Type"
    }
    
    # Define a function to calculate average gene expression
    analyze_gene_expression <- function(seurat_obj, genes_chunk, cell_type_column, sample_column) {
        # Fetch data for the specified genes, cell types, and samples
        expression_data <- Seurat::FetchData(seurat_obj, vars = c(genes_chunk, cell_type_column, sample_column))
        avg_exp_cols <- paste0("avg.exp_", genes_chunk)
        
        # Calculate the average expression for each gene across samples and cell types
        expression_data %>%
            dplyr::group_by(!!rlang::sym(sample_column), !!rlang::sym(cell_type_column)) %>%
            dplyr::summarize(dplyr::across(dplyr::all_of(genes_chunk), ~ mean(expm1(.), na.rm = TRUE), .names = "avg.exp_{.col}"), .groups = 'drop') %>%
            tidyr::pivot_longer(cols = dplyr::all_of(avg_exp_cols), names_to = "Gene", values_to = "avg.exp") %>%
            dplyr::mutate(Gene = sub("avg.exp_", "", Gene))
    }
    
    # Initialize the results dataframe and progress bar
    total_chunks <- ceiling(length(gene_list) / chunk_size)
    pb <- progress::progress_bar$new(
        format = "  Analyzing genes [:bar] :percent eta: :eta",
        total = total_chunks, clear = FALSE, width = 60
    )
    
    # Process genes in chunks and calculate average expression
    final_result <- dplyr::bind_rows(lapply(seq(1, length(gene_list), by = chunk_size), function(start_idx) {
        pb$tick()
        end_idx <- min(start_idx + chunk_size - 1, length(gene_list))
        genes_chunk <- gene_list[start_idx:end_idx]
        analyze_gene_expression(seurat_obj, genes_chunk, cell_type_column, sample_column)
    }))
    
    # Add patient information to the results dataframe
    patient_column <- seurat_obj@meta.data[[group_column]][match(final_result[[sample_column]], seurat_obj@meta.data[[sample_column]])]
    final_result <- cbind(final_result, Patient = patient_column)
    
    # Define a function for ANOVA and Tukey's test for a single gene and cell type
    anova_single_gene_celltype <- function(data, gene, cell_type, cell_type_column) {
        anova_data <- dplyr::filter(data, Gene == gene, !!rlang::sym(cell_type_column) == cell_type)
        
        # Perform ANOVA and Tukey's test if there are at least 2 groups with data
        if (dplyr::n_distinct(anova_data$Patient) < 2) return(NULL)
        
        tryCatch({
            aov_result <- stats::aov(avg.exp ~ Patient, data = anova_data)
            tukey_result <- stats::TukeyHSD(aov_result)
            list(anova = summary(aov_result), tukey = tukey_result)
        }, error = function(e) NULL)
    }
    
    # Perform ANOVA and Tukey's test for all genes and cell types
    anova_all_genes <- function(data, gene_list, cell_types, cell_type_column) {
        pb <- progress::progress_bar$new(
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
        dplyr::bind_rows(results)
    }
    
    # Get unique cell types for analysis
    cell_types <- unique(final_result[[cell_type_column]])
    
    # Perform ANOVA and Tukey's test for all genes and cell types
    anova_tukey_results_df <- anova_all_genes(final_result, gene_list, cell_types, cell_type_column)
    
    # Get unique groups for pairwise comparisons
    unique_groups <- unique(seurat_obj@meta.data[[group_column]])
    group_combinations <- utils::combn(unique_groups, 2, simplify = FALSE)
    
    # Initialize list to hold FC results
    fc_results_list <- list()
    
    # Calculate FC for each pairwise group combination
    total_combinations <- length(group_combinations)
    pb_fc <- progress::progress_bar$new(
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
    all_fc_results <- dplyr::bind_rows(fc_results_list)
    
    # Merge FC results with ANOVA and Tukey's test results
    if (nrow(all_fc_results) > 0) {
        anova_tukey_results_df$Comparison <- as.character(anova_tukey_results_df$Comparison)
        all_fc_results$Comparison <- as.character(all_fc_results$Comparison)
        
        final_merged_output <- merge(anova_tukey_results_df, all_fc_results, by = c("Gene", "Cell_Type", "Comparison"), all = TRUE)
        
        # Remove the Group column from the final output
        final_merged_output <- dplyr::select(final_merged_output, -Group)
    } else {
        final_merged_output <- anova_tukey_results_df
    }
    
    return(final_merged_output)
}

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
        seurat_obj <- Seurat::subset(seurat_obj, subset = seurat_obj[[cell_type_column]] == cell_type)
        cell_type_name <- cell_type
    } else {
        cell_type_name <- "All_Cells"
    }
    
    # Set identities to the grouping column
    Seurat::Idents(seurat_obj) <- group_column
    
    # Select cells from the two groups
    cells.1 <- Seurat::WhichCells(seurat_obj, ident = ident.1)
    cells.2 <- Seurat::WhichCells(seurat_obj, ident = ident.2)
    
    # Get the expression matrix from the specified slot
    data <- Seurat::GetAssayData(object = seurat_obj, slot = slot)
    
    # If features is NULL, use all genes
    if (is.null(features)) {
        features <- rownames(data)
    }
    
    # Select the data for the specified cells and features
    data.1 <- data[features, cells.1, drop = FALSE]
    data.2 <- data[features, cells.2, drop = FALSE]
    
    # Define the mean function with apply to replace rowMeans
    scGeneANOVA_fc <- function(x, pseudocount.use, base) {
        if (is.vector(x)) {
            x <- matrix(x, nrow = 1)
        }
        return(log(apply(x, 1, mean, na.rm = TRUE) + pseudocount.use, base = base))
    }
    
    # Calculate average expression for both groups
    avg.exp.1 <- scGeneANOVA_fc(data.1, pseudocount.use, base)
    avg.exp.2 <- scGeneANOVA_fc(data.2, pseudocount.use, base)
    
    # Calculate fold change
    fold_changes <- avg.exp.1 - avg.exp.2
    
    # Create the result dataframe
    fc.results <- data.frame(
        Gene = features,
        avg_logFC = fold_changes,
        pct.1 = apply(data.1 > 0, 1, mean),
        pct.2 = apply(data.2 > 0, 1, mean)
    )
    
    # Add group identities and cell type to the results
    if (nrow(fc.results) > 0) {
        fc.results$Group <- paste(ident.1, "vs", ident.2)
        fc.results$Cell_Type <- cell_type_name
        fc.results$Comparison <- gsub(" vs ", "-", fc.results$Group)
    }
    
    return(fc.results)
}
