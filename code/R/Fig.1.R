# Load necessary libraries
library(readr)
library(GSVA)
library(dplyr)
library(tidyr)
library(corto)
library(depmap)
library(ggplot2)
library(ggrepel)

# Get expression data from DepMap
merged_data <- depmap::depmap_TPM()
# Create a named list with gene categories


# Fig.1A----------------------------------------------------------------


# Function to calculate scores using single sample gene set enrichment analysis (ssgsea)
calculate_scores <- function() {
  # Select relevant columns and pivot data
  data_final <- select(merged_data, -entrez_id, -gene, -cell_line)
  ssgsea_data <- data_final %>%
    pivot_wider(names_from = depmap_id, values_from = rna_expression)
  
  # Transform data for ssgsea
  ssgsea_data <- as.data.frame(ssgsea_data)
  rownames(ssgsea_data) <- as.character(ssgsea_data$gene_name)
  
  # Run ssgsea
  
  ssgsea_data <- as.matrix(ssgsea_data[, -1])
  param <- GSVA::ssgseaParam(ssgsea_data, hallmark_EMT_genes)
  nesmat <- GSVA::gsva(param)
  colnames(nesmat) <- colnames(ssgsea_data)
  result_matrix <- rbind(nesmat, ssgsea_data)
  
  # Clean up intermediate data
  rm(data_final, ssgsea_data)
  
  return(result_matrix)
}

# Calculate scores using the defined function
result_matrix <- calculate_scores()
result_matrix <- as.matrix(result_matrix)

# Initialize a data frame to store correlation results
correlation_matrix <- data.frame(gene = character(), mes_correlation = double(), epi_correlation = double(), stringsAsFactors = FALSE)
num_rows <- nrow(result_matrix)

# Calculate correlations for each gene with mesenchymal and epithelial gene sets
for (row_index in 3:num_rows) {
  mes_correlation <- cor(result_matrix[3, ], result_matrix[row_index,], method = "spearman")
  epi_correlation <- cor(result_matrix[2, ], result_matrix[row_index,], method = "spearman")
  gene_name <- rownames(result_matrix)[row_index]
  
  # Append results to the correlation matrix
  correlation_matrix[nrow(correlation_matrix) + 1, ] <- c(gene_name, mes_correlation, epi_correlation)
}

# Convert correlation values to numeric
correlation_matrix$mes_correlation <- as.numeric(correlation_matrix$mes_correlation)
correlation_matrix$epi_correlation <- as.numeric(correlation_matrix$epi_correlation)
correlation_matrix <- correlation_matrix %>%
  mutate(color = case_when(
    epi_correlation < 0 & mes_correlation > 0 ~ "blue",
    epi_correlation > 0 & mes_correlation < 0 ~ "yellow",
    epi_correlation < 0 & mes_correlation < 0 ~ "green",
    epi_correlation > 0 & mes_correlation > 0 ~ "purple",
    TRUE ~ "grey"
  ))
# Genes to label in the scatterplot
genes_to_label <- c("ELF3","KLF4","SNAI2","VIM","SNAI1","ZEB1","CDH1","GRHL2","OVOL2", "IRF1","IRF2","IRF3","IRF4","IRF7","IRF8","IRF9")


Fig1a <- ggplot(correlation_matrix, aes(x = epi_correlation, y = mes_correlation, label = gene)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_point(color = "grey") +
  geom_point(color = "orange", data = subset(correlation_matrix, (color == "yellow" & gene %in% genes_to_label))) +
  geom_point(color = "darkblue", data = subset(correlation_matrix, (color == "blue" & gene %in% genes_to_label))) +
  geom_point(color = "darkgreen", data = subset(correlation_matrix, (color == "green" & gene %in% genes_to_label))) +
  geom_point(color = "purple", data = subset(correlation_matrix, (color == "purple" & gene %in% genes_to_label))) +
  geom_point(color = "red", data = subset(correlation_matrix, gene == "IRF6")) +
  geom_text_repel(
    aes(color = gene),
    data = subset(correlation_matrix, gene %in% genes_to_label),
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'black',
    segment.size = 0.2,
    color = "black"
  ) +
  geom_text_repel(
    data = subset(correlation_matrix, gene == "IRF6"),
    size = 5,
    box.padding = 3,
    point.padding = 2,
    segment.color = 'red',
    segment.size = 1,
    color = "red",
    fontface = "bold",
    xlim = c(0.4, 0.55),
    ylim = c(-0.2, 0.55)
  ) +
  labs(x = "epithelial scores",
       y = "mesenchymal scores") +
  theme_minimal() +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 15, family = "Arial"),
    axis.text  = element_text(size = 14, family = "Arial"),
    plot.title = element_text(family = "Arial"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # <-- fixed
  ) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.2, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.2, linetype = "dashed")

Fig1a
# Fig.1B ------------------------------------------------------------------
library(pheatmap)
ccle_data <- merged_data[, c(3, 5, 1)]
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

# --- Compute correlation and p-value matrices ---
cor_matrix <- matrix(NA,
                     nrow = length(irf_genes),
                     ncol = length(marker_genes),
                     dimnames = list(irf_genes, marker_genes))

pval_matrix <- cor_matrix

for (i in seq_along(irf_genes)) {
  for (j in seq_along(marker_genes)) {
    if (irf_genes[i] %in% rownames(filtered_data) &
        marker_genes[j] %in% rownames(filtered_data)) {
      
      test <- cor.test(
        as.numeric(filtered_data[irf_genes[i], ]),
        as.numeric(filtered_data[marker_genes[j], ]),
        method = "spearman"
      )
      
      cor_matrix[i, j]  <- test$estimate
      pval_matrix[i, j] <- test$p.value
    }
  }
}
# --- FDR correction ---
pval_fdr <- matrix(
  p.adjust(as.vector(pval_matrix), method = "BH"),
  nrow = nrow(pval_matrix),
  dimnames = dimnames(pval_matrix)
)

## To use without FDR uncomment 
#pval_fdr <- pval_matrix

# --- Convert FDR p-values to stars ---
stars_fdr <- ifelse(pval_fdr < 0.001, "***",
                    ifelse(pval_fdr < 0.01, "**",
                           ifelse(pval_fdr < 0.05, "*", "")))

# --- Display matrix (correlation + FDR stars) ---
display_fdr <- matrix(
  paste0(stars_fdr),
  nrow = nrow(cor_matrix),
  dimnames = dimnames(cor_matrix)
)

# --- Plot ---
Fig1b <- pheatmap(
  cor_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = display_fdr,
  fontsize = 14,
  angle_col = "45",
  treeheight_row = 10,
  treeheight_col = 10,
  height = 3,
  width = 7,
  legend = TRUE,
  filename = "~/Work/IRF6_codes_figs/figures/Fig1b.png"
)
Fig1b


# Fig.1C -----------------------------------------------------------------

# Define the function ELF3_vs_median
ELF3_vs_median <- function(cancer_type, find_gene) {
  # Generate file path based on cancer type
  file_path <- paste0("/home/aditya/Work/data/TCGA_datasets/TCGA.", toupper(cancer_type), ".sampleMap_HiSeqV2.gz")
  
  # Read the data
  TCGA_data <- read_tsv(file_path)
  
  
  
  # Prepare TCGA data for analysis
  TCGA_data <- as.data.frame(TCGA_data)
  rownames(TCGA_data) <- TCGA_data[, 1]
  TCGA_data <- TCGA_data[, -1]
  
  # Perform ssGSEA
  
  param <- GSVA::ssgseaParam(as.matrix(TCGA_data), hallmark_EMT_genes)
  nesmat <- GSVA::gsva(param)
  #nesmat <- ssgsea(as.matrix(TCGA_data), hallmark_EMT_genes, scale = TRUE, minsize = 4)
  colnames(nesmat) <- colnames(TCGA_data)
  nesmat <- rbind(nesmat,TCGA_data[find_gene, ])
  nesmat <- as.matrix(nesmat)
  epi_cor <- cor(nesmat[2, ], nesmat[7, ], method = "spearman")
  mes_cor <- cor(nesmat[3, ], nesmat[7, ], method = "spearman")
  # Create a result data frame
  result_df <- data.frame(
    type = cancer_type,
    epi_correlation = epi_cor,
    mes_correlation = mes_cor,
    stringsAsFactors = FALSE
  )
  return(result_df)
}

cancer_types <- c("BRCA", "ACC", "BLCA", "CESC", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUSC", "MESO", "OV", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PAAD", "LUAD", "LUNG")
genes <- c("IRF6")
for (find_gene in genes) {
  result_df <- data.frame(type = character(), epi_correlation = double(), mes_correlation = double(), stringsAsFactors = FALSE)
  for (cancer_type in cancer_types) {
    # Get the correlation results for the current cancer type and gene
    result <- ELF3_vs_median(cancer_type, find_gene = "IRF6")
    
    # Append the results to the result_df
    result_df <- rbind(result_df, result)
  }
}
# Define the list of cancer types to label
label_list <- c("BRCA", "LUAD", "OV", "COAD", "LUSC","PAAD","STAD")

# Create a new column to indicate if the point should be labeled
result_df$label <- ifelse(result_df$type %in% label_list, result_df$type, "")

# Create a new column to specify the color
result_df$color <- ifelse(result_df$type %in% label_list, result_df$type, "unlabeled")

# Create a scatterplot with correlation values
Fig1c <- ggplot(result_df, aes(x = epi_correlation, y = mes_correlation)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3)+
  geom_point(aes(color = color)) +
  geom_text_repel(
    aes(label = ifelse(label != "", label, NA), color = color),
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'black',
    segment.size = 0.2,
    max.overlaps = Inf
  ) +
  labs(x = "epithelial scores",
       y = "mesenchymal scores") +
  theme_minimal() +
  scale_x_continuous(limits = c(-0.6, 1)) +
  scale_y_continuous(limits = c(-1, 0.6))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 15, family = "Arial"),
    axis.text  = element_text(size = 14, family = "Arial"),
    plot.title = element_text(family = "Arial"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Remove the legend
  ) +
  scale_color_manual(values = c(setNames(rainbow(length(unique(result_df$type[result_df$type %in% label_list]))), unique(result_df$type[result_df$type %in% label_list])), "unlabeled" = "grey"))


# Fig.1D ------------------------------------------------------------------

library(affy)
library(hgu133plus2.db)
source("~/Work/microarrayexpr.R")
source("~/Work/RNAseqexpr.R")
source("~/Work/CountsToTPM/CountsToTPM/count_to_TPM.R")
library(AnnotationDbi)
# Load all datasets
GSE36081 <- get_expression_matrix_microarray("GSE36081")
GSE118407 <- get_expr_htp_human("GSE118407")
GSE61220 <- get_expr_htp_human("GSE61220") 
GSE59922 <- get_expression_matrix_microarray("GSE59922")
GSE40690 <- get_expression_matrix_microarray("GSE40690")
cel_files_path <- "/home/aditya/Downloads/GSE58252_RAW/"
data <- ReadAffy(celfile.path = cel_files_path)
normalized_data <- rma(data)
GSE58252 <- exprs(normalized_data)
probe_ids <- rownames(GSE58252)
gene_symbols <- mapIds(hgu133plus2.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
rownames(GSE58252) <- gene_symbols

# Put all expression matrices into a named list
expr_list <- list(
  GSE36081 = GSE36081,
  GSE118407 = GSE118407,
  GSE61220 = GSE61220[,c(1,2,3,10,11,12)],
  GSE59922 = GSE59922,
  GSE40690 = GSE40690[,c(1:9)],
  GSE58252 = GSE58252
)

# Function to extract IRF6 expression from a dataset
get_IRF6_expression <- function(expr_mat) {
  if (!"IRF6" %in% rownames(expr_mat)) stop("IRF6 not found in dataset")
  t(expr_mat["IRF6", , drop = FALSE])  # return as matrix with samples as columns
}

# Extract IRF6 expression for all datasets
IRF6_expr <- lapply(expr_list, get_IRF6_expression)

library(ggplot2)
library(dplyr)
library(patchwork)

# JCO palette (extend as needed)
jco_palette <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF")

# Generalized function for multiple groups
plot_expression_bar_multi <- function(expr_matrix, group_names, group_sizes, title = NULL) {
  if(length(group_names) != length(group_sizes)) stop("group_names and group_sizes must be same length")
  
  df <- data.frame(
    expr = as.numeric(expr_matrix),
    sample = colnames(expr_matrix),
    stringsAsFactors = FALSE
  )
  
  # Assign groups
  df$group <- rep(group_names, times = group_sizes)
  
  # Summary for bar plot
  summary_df <- df %>%
    group_by(group) %>%
    summarise(
      mean_expr = mean(expr),
      sd_expr = sd(expr),
      .groups = "drop"
    )
  summary_df <- summary_df[match(c("control", setdiff(summary_df$group, "control")), summary_df$group), ]
  # Pairwise t-tests vs first group (control)
  sig_texts <- sapply(2:length(group_names), function(i) {
    ttest <- t.test(expr ~ group, data = df %>% filter(group %in% c(group_names[1], group_names[i])))
    p <- ttest$p.value
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  })
  
  # Compute positions for stars
  sig_df <- summary_df[-1, ]   # skip control
  sig_df$sig <- sig_texts
  sig_df$y_pos <- sig_df$mean_expr + sig_df$sd_expr + 0.2  # place above the error bar
  
  # Then use geom_text with proper x
  
  
  
  # Plot
  ggplot(summary_df, aes(x = group, y = mean_expr, fill = group)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_errorbar(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr), width = 0.2) +
    scale_fill_manual(values = jco_palette[1:length(group_names)]) +
    # Add significance above each treatment (skip control)
    geom_text(data = sig_df, aes(x = group, y = y_pos, label = sig), size = 6) +
    labs(title = title, y = "IRF6 Expression", x = "") +
    theme_minimal()+
    theme(legend.position = "none", axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))
}

# Plot 1–5
plot1 <- plot_expression_bar_multi(IRF6_expr$GSE36081,
                                   group_names = c("control", "GRHL2"),
                                   group_sizes = c(3, 3),
                                   title = "GSE36081")

plot2 <- plot_expression_bar_multi(IRF6_expr$GSE118407,
                                   group_names = c("control", "shGRHL2"),
                                   group_sizes = c(4, 4),
                                   title = "GSE118407")

plot3 <- plot_expression_bar_multi(IRF6_expr$GSE61220,
                                   group_names = c("control", "TGFB"),
                                   group_sizes = c(3, 3),
                                   title = "GSE61220")

plot4 <- plot_expression_bar_multi(IRF6_expr$GSE59922,
                                   group_names = c("control", "TGFB"),
                                   group_sizes = c(3, 3),
                                   title = "GSE59922")

# Plot 5 has 2 treatments
plot5 <- plot_expression_bar_multi(IRF6_expr$GSE40690,
                                   group_names = c("control", "SNAI1", "SNAI2"),
                                   group_sizes = c(3, 3, 3),
                                   title = "GSE40690")
plot6 <- plot_expression_bar_multi(IRF6_expr$GSE58252,
                                   group_names = c("SNAI1","control"),
                                   group_sizes = c(3, 3),
                                   title = "GSE58252")

# Plotting ----------------------------------------------------------------


library(patchwork)

# Use your layout
layout <- c(
  area(1, 1, 3, 6),  # Fig1a
  area(4, 1, 6, 6),  # Fig1b_wrapped
  area(7, 1, 9, 6),  # Fig1c
  area(1, 7, 9, 16)   # Fig1d_combined
)
#plot(layout)
Fig1b_wrapped <- wrap_elements(Fig1b$gtable)
Fig1d_combined <- wrap_plots(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 3) |> 
  wrap_elements()
# Combine all plots using the custom layout
Fig1 <- Fig1a + Fig1b_wrapped + Fig1c + Fig1d_combined +
  plot_layout(design = layout) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 30)  # <- panel letters
    )
  )

Fig1

