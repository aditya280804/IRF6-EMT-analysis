# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(corto)
library(depmap)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(GSVA)

# Get expression data from DepMap
ccle_data <- depmap::depmap_TPM()
data_final <- dplyr::select(ccle_data, -entrez_id, -gene, -cell_line)
ssgsea_data <- data_final %>%
  pivot_wider(names_from = depmap_id, values_from = rna_expression)
ssgsea_data <- as.data.frame(ssgsea_data)
rownames(ssgsea_data) <- as.character(ssgsea_data$gene_name)
# nesmat <- ssgsea(as.matrix(ssgsea_data[, -1]), hallmark_EMT_genes, scale = TRUE, minsize = 4)
tpm <- as.matrix(ssgsea_data[-1])
param <- GSVA::ssgseaParam(tpm,hallmark_EMT_genes)
nesmat <- GSVA::gsva(param)
colnames(nesmat)<- colnames(ssgsea_data[,-1])
nesmat <- rbind(ssgsea_data[c("IRF6","CD274","ESR1"),-1],nesmat)
nesmat <- t(nesmat)
nesmat <- as.data.frame(nesmat)
gene <- "IRF6"
datasets <- c("CD274","ESR1","ER_EARLY","ER_LATE","PDL1")
plot_list <- list()

# Loop through the datasets
for(i in datasets) {
  dataset <- i
  
  # Calculate Spearman correlation
  correlation_value <- cor(nesmat[[gene]], nesmat[[dataset]], method = "spearman")
  
  # Perform Welch's t-test (two-tailed)
  ttest_result <- t.test(nesmat[[gene]], nesmat[[dataset]], var.equal = FALSE)
  
  # Extract p-value from the t-test result
  p_value <- ttest_result$p.value
  
  # Create the plot
  p <- ggplot(nesmat, aes_string(x = gene, y = dataset)) +
    geom_point() +  # Add points to the plot
    labs(title = paste(gene, "vs", dataset,"\np =", round(p_value, 10))) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Add annotation with correlation value and p-value
  p <- p + annotate("text", x = max(nesmat[[gene]], na.rm = TRUE), y = max(nesmat[[dataset]], na.rm = TRUE),
                    label = paste("R =", round(correlation_value, 2)), hjust = 1.2, vjust = 1.2, size = 5)
  
  # Store the plot in the list
  plot_list[[dataset]] <- p
}

# Combine all plots using patchwork
combined_plots <- wrap_plots(plot_list)

# Print the combined plot
print(combined_plots)
