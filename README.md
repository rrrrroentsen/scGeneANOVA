# scGeneANOVA
## Introduction
`scGeneANOVA` is an R package for performing gene expression analysis and ANOVA on single-cell RNA-seq data using Seurat.

## Installation
To install `scGeneANOVA`, you can use the `devtools` package to install directly from GitHub:
# Install devtools if you haven't
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
# Install scGeneANOVA from GitHub
devtools::install_github("rrrrroentsen/scGeneANOVA")

## Usage
# Load necessary libraries
library(Seurat)
library(scGeneANOVA)

# Example usage with a Seurat object
seurat_obj <- readRDS("path/to/your/seurat_obj.rds")
results <- scGeneANOVA(seurat_obj, gene_list = NULL, cell_type_column = "Cell_Type", group_column = "Patient_ident", sample_column = "orig.ident")
# View ANOVA and Tukey's test results
print(results)

# Function Parameters
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
