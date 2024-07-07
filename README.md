# scGeneANOVA
## Introduction
`scGeneANOVA` is an R package for performing gene expression analysis and ANOVA on single-cell RNA-seq data using Seurat.

## Installation
To install `scGeneANOVA`, you can use the `devtools` package to install directly from GitHub:

devtools::install_github("rrrrroentsen/scGeneANOVA")

## Usage

library(Seurat)

library(tidyr)

library(dplyr)

library(progress)

library(scGeneANOVA)

seurat_obj <- readRDS("path/to/your/seurat_obj.rds")

gene_list <- c("CD40", "CD80", "CD83")

results <- scGeneANOVA(seurat_obj, gene_list = gene_list, cell_type_column = "Cell_Type", group_column = "Patient_ident", sample_column = "orig.ident")

print(results)

write.csv(results, "anova_tukey_results.csv", row.names = FALSE)

## Function Parameters
•	seurat_obj: A Seurat object containing single-cell RNA-seq data.

•	gene_list: A list of genes to analyze. If NULL, all genes are used.

•	cell_type_column: The column name in the metadata that contains cell type information. If NULL, all cells are treated as one group.

•	group_column: The column name in the metadata that contains group information.

•	sample_column: The column name in the metadata that contains sample information.

•	chunk_size: The number of genes to process in each chunk. Default is 100.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
•	Built using Seurat and ggplot2.

•	Inspired by the need for comprehensive single-cell RNA-seq data analysis.

## Contributing
Contributions are welcome! Please fork the repository and create a pull request for any enhancements or bug fixes.
