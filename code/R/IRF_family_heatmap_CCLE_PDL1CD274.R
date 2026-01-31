library(readr)
library(tidyr)
library(dplyr)
library(corto)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(pheatmap)
library(gridExtra)

# Load t# Load t# Load t# Load t# Load the data
merged_data <- read_csv("/home/aditya/Work/data/epimesscores,tissue/merged_data.csv")

# Get the unique lineages
lineages <- unique(merged_data$lineage)
lineages <- setdiff(lineages, c("epidermoid_carcinoma", "adrenal_cortex", "unknown"))
genes <- c("IRF6","IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9")
datasets <- c("PDL1", "CD274","hallmark_EMT","epithelial","mesenchymal")
plot_list <- list()
for (i in 1:length(genes)) {
  # Generate the dataframe name
  df_name <- paste0("df_IRF", i)
  corr_df_name <- paste0("signif_df_IRF", i)
  # Create an empty dataframe with specified column names and row names
  df <- data.frame(matrix(nrow = length(datasets), ncol = length(lineages)))
  colnames(df) <- lineages
  rownames(df) <- datasets
  
  # Assign the dataframe to a variable with the generated name
  assign(df_name, df)
  df_2 <- data.frame(matrix(nrow = length(datasets), ncol = length(lineages)))
  colnames(df_2) <- lineages
  rownames(df_2) <- datasets
  
  # Assign the dataframe to a variable with the generated name
  assign(corr_df_name, df_2)
}
for (i in lineages){
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
    for (gene in genes) {
      for (dataset in datasets) {
        correlation_value <- cor(nesmat[[gene]], nesmat[[dataset]], method = "spearman")
        df_name <- paste0("df_", gene)
        df <- get(df_name)
        df[dataset, i] <- correlation_value
        assign(df_name, df)
        significance_value <- t.test(nesmat[[gene]], nesmat[[dataset]],var.equal = FALSE)$p.value
        df_name <- paste0("signif_df_", gene)
        df <- get(df_name)
        df[dataset, i] <- significance_value
        assign(df_name, df)
      }
    }
    cat("Finished processing lineage:", i, "\n")
}

# Initialize a list to store the heatmap grobs
heatmap_list <- list()
for (i in  1:9) {
  df_name <- paste0("df_", genes[i])
  p_df_name <- paste0("signif_df_", genes[i])
  variable <- get(df_name)
  signif_variable <- get(p_df_name)
  signif_variable <- apply(signif_variable, c(1, 2), format_pvalue)
  
  p <- pheatmap(as.matrix(variable), 
                main = paste0(genes[i], " Correlations"), 
                display_numbers = signif_variable,  # Do not display correlation values in cells
                cluster_rows = TRUE,      # Cluster rows
                cluster_cols = FALSE,  
                fontsize = 10,# Cluster column
                color = colorRampPalette(c("navy", "white", "firebrick3"))(100))  # Color scheme
  heatmap_list[[i]] <- p$gtable
}

# Arrange the heatmaps using grid.arrange
do.call(grid.arrange, c(heatmap_list, ncol = 1))

