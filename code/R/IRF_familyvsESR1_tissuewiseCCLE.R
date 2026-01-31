library(readr)
library(tidyr)
# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(corto)
library(depmap)
library(ggplot2)
library(patchwork)
library(ggrepel)
merged_data <- read_csv("/home/aditya/Work/data/epimesscores,tissue/merged_data.csv")
merged_data_breast <- subset(merged_data, lineage == "breast")
merged_data_breast <- merged_data_breast[,c(1,3,5)]
merged_data_breast <- pivot_wider(merged_data_breast,names_from = depmap_id,values_from = rna_expression)
merged_data_breast <- as.data.frame(merged_data_breast)
rownames(merged_data_breast) <- merged_data_breast[,1]


nesmat <- ssgsea(as.matrix(merged_data_breast[, -1]), hallmark_EMT_genes, scale = TRUE, minsize = 4)
colnames(nesmat)<- colnames(merged_data_breast[,-1])
nesmat <- rbind(merged_data_breast[c("IRF6","CD274","ESR1","IRF1","IRF2","IRF3","IRF4","IRF5","IRF7","IRF8","IRF9"),-1],nesmat)
nesmat <- t(nesmat)
nesmat <- as.data.frame(nesmat)
genes <- c("IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9")
datasets <- c("ER_EARLY", "ER_LATE")

# Initialize list to store plots
plot_list <- list()

# Loop through genes and datasets, create plots, and store them
for (dataset in datasets) {
  for (gene in genes) {
    correlation_value <- cor(nesmat[[gene]], nesmat[[dataset]], method = "spearman")
    ttest_result <- t.test(nesmat[[gene]], nesmat[[dataset]], var.equal = FALSE)
    p_value <- ttest_result$p.value
    
    p <- ggplot(nesmat, aes_string(x = gene, y = dataset)) +
      geom_point() +
      labs(title = paste(gene, "vs", dataset, "\np-value =", round(p_value, 10))) +
      theme(plot.title = element_text(hjust = 0.5)) +
      annotate("text", x = max(nesmat[[gene]], na.rm = TRUE), y = max(nesmat[[dataset]], na.rm = TRUE),
               label = paste("Spearman correlation =", round(correlation_value, 2)), hjust = 1.2, vjust = 1.2, size = 5)
    
    plot_list[[paste(gene, dataset, sep = "_")]] <- p
  }
}

# Combine all plots using patchwork and print
wrap_plots(plot_list)
