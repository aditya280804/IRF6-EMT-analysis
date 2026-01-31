library(depmap)
library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)

# Load data
ccle_data_initial <- depmap::depmap_TPM()
ccle_data <- ccle_data_initial[, c(3, 5, 1)]
ccle_data <- ccle_data %>%
  pivot_wider(names_from = depmap_id, values_from = rna_expression)
ccle_data <- as.matrix(ccle_data)
rownames(ccle_data) <- ccle_data[, 1]
ccle_data <- ccle_data[, -1]

# Genes of interest
irf_genes <- paste0("IRF", 1:9)
marker_genes <- c("ZEB1", "ZEB2", "ELF3", "KLF4", "OVOL1", "OVOL2", 
                  "GRHL2", "SNAI1", "SNAI2", "CD274", "ESR1", "CDH1")

genes_to_extract <- unique(c(irf_genes, marker_genes))

# Filter only available genes
genes_available <- intersect(genes_to_extract, rownames(ccle_data))
filtered_data <- ccle_data[genes_available, ]

# Calculate Spearman correlation
cor_matrix <- matrix(nrow = length(irf_genes), ncol = length(marker_genes))
rownames(cor_matrix) <- irf_genes
colnames(cor_matrix) <- marker_genes

for (i in seq_along(irf_genes)) {
  for (j in seq_along(marker_genes)) {
    if (irf_genes[i] %in% rownames(filtered_data) & marker_genes[j] %in% rownames(filtered_data)) {
      cor_matrix[i, j] <- cor(as.numeric(filtered_data[irf_genes[i], ]),
                              as.numeric(filtered_data[marker_genes[j], ]),
                              method = "spearman")
    } else {
      cor_matrix[i, j] <- NA
    }
  }
}

# Plot heatmap
pheatmap(cor_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = TRUE, 
         main = "Correlation between IRF1-9 and EMT Markers",
         angle_col = 0)
