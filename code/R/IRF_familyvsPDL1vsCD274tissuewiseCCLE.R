library(readr)
library(tidyr)
library(dplyr)
library(corto)
library(ggplot2)
library(patchwork)
library(ggrepel)

# Load the data
merged_data <- read_csv("/home/aditya/Work/data/epimesscores,tissue/merged_data.csv")

# Get the unique lineages
lineages <- unique(merged_data$lineage)
lineages <- setdiff(lineages, c("epidermoid_carcinoma", "unknown", "adrenal_cortex"))
# Define genes and datasets
genes <- c("IRF6")
datasets <- c("PDL1", "CD274")

# Initialize list to store plots for all lineages
all_plots <- list()

# Loop through each lineage
for (i in lineages) {
  # Subset data for the current lineage
  merged_data_lineage <- subset(merged_data, lineage == i)
  merged_data_lineage <- merged_data_lineage[,c(1,3,5)]
  merged_data_lineage <- pivot_wider(merged_data_lineage, names_from = depmap_id, values_from = rna_expression)
  merged_data_lineage <- as.data.frame(merged_data_lineage)
  rownames(merged_data_lineage) <- merged_data_lineage[,1]
  
  # Perform ssGSEA
  nesmat <- ssgsea(as.matrix(merged_data_lineage[, -1]), hallmark_EMT_genes, scale = TRUE, minsize = 4)
  colnames(nesmat) <- colnames(merged_data_lineage[,-1])
  nesmat <- rbind(merged_data_lineage[c("IRF6","CD274","ESR1","IRF1","IRF2","IRF3","IRF4","IRF5","IRF7","IRF8","IRF9"),-1], nesmat)
  nesmat <- t(nesmat)
  nesmat <- as.data.frame(nesmat)
  
  # Initialize list to store plots for the current lineage
  plot_list <- list()
  
  # Loop through genes and datasets, create plots, and store them
  for (dataset in datasets) {
    for (gene in genes) {
      correlation_value <- cor(nesmat[[gene]], nesmat[[dataset]], method = "spearman")
      ttest_result <- t.test(nesmat[[gene]], nesmat[[dataset]], var.equal = FALSE)
      p_value <- ttest_result$p.value
      
      p <- ggplot(nesmat, aes_string(x = gene, y = dataset)) +
        geom_point() +
        labs(title = paste(gene, "vs", dataset, "\nLineage:", i, "\np-value =", round(p_value, 10))) +
        theme(plot.title = element_text(hjust = 0.5)) +
        annotate("text", x = max(nesmat[[gene]], na.rm = TRUE), y = max(nesmat[[dataset]], na.rm = TRUE),
                 label = paste("Spearman correlation =", round(correlation_value, 2)), hjust = 1.2, vjust = 1.2, size = 5)
      
      plot_list[[paste(i, gene, dataset, sep = "_")]] <- p
    }
  }
  
  # Combine all plots for the current lineage and store them
  all_plots[[i]] <- wrap_plots(plot_list)
  # Update message after processing each lineage
  cat("Finished processing lineage:", i, "\n")
}

# Combine all plots for all lineages and print
combined_plot <- wrap_plots(all_plots)
pdf(paste0("plots_IRF6.pdf"), width = 20, height = 16)
print(combined_plot)
dev.off()