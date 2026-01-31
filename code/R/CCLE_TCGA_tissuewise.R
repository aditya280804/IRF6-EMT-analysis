library(readr)
library(dplyr)
library(tidyr)
library(corto)
library(GSVA)
library(ggplot2)
library(ggrepel)

source("~/Work/EMT_signatures_storage.R")

# -----------------------------
# Parameters
# -----------------------------
gene_name <- "IRF6"

merged_data <- read_csv(
  "/home/aditya/Work/data/epimesscores,tissue/merged_data.csv"
)

# Keep only required columns
merged_data <- merged_data[, c(1, 3, 5, 7)]

# -----------------------------
# Initialize results table
# -----------------------------
correlation_table <- data.frame(
  tissue = character(),
  epi_rho = double(),
  epi_p = double(),
  mes_rho = double(),
  mes_p = double(),
  stringsAsFactors = FALSE
)

# -----------------------------
# Loop over tissues
# -----------------------------
tissue_list <- unique(merged_data$lineage)
tissue_list <- tissue_list[c(-28, -29)]   # same exclusions as before

for (tissue in tissue_list) {
  message("Processing tissue: ", tissue)
  
  temp <- merged_data %>%
    filter(lineage == tissue) %>%
    select(-lineage)
  
  temp <- pivot_wider(
    temp,
    names_from = depmap_id,
    values_from = rna_expression
  )
  
  temp <- as.data.frame(temp)
  rownames(temp) <- temp[, 1]
  temp <- temp[, -1]
  
  # ssGSEA
  param <- GSVA::ssgseaParam(as.matrix(temp), hallmark_EMT_genes)
  nesmat <- GSVA::gsva(param)
  #nesmat <- ssgsea(temp, EM_genesets_cell_lines)
  colnames(nesmat) <- colnames(temp)
  
  # Append IRF6 expression
  nesmat <- rbind(nesmat, temp[gene_name, ])
  nesmat <- as.matrix(nesmat)
  
  # -----------------------------
  # Correlations + p-values
  # -----------------------------
  epi_test <- cor.test(
    nesmat["epithelial", ],
    nesmat[gene_name, ],
    method = "spearman"
  )
  
  mes_test <- cor.test(
    nesmat["mesenchymal", ],
    nesmat[gene_name, ],
    method = "spearman"
  )
  
  correlation_table[nrow(correlation_table) + 1, ] <- list(
    tissue  = tissue,
    epi_rho = epi_test$estimate,
    epi_p   = epi_test$p.value,
    mes_rho = mes_test$estimate,
    mes_p   = mes_test$p.value
  )
}

# -----------------------------
# Final formatting
# -----------------------------
result_df <- correlation_table %>%
  mutate(
    epi_rho = round(epi_rho, 3),
    mes_rho = round(mes_rho, 3),
    epi_p   = signif(epi_p, 3),
    mes_p   = signif(mes_p, 3)
  )

# -----------------------------
# Save correlation table
# -----------------------------
write_tsv(
  result_df,
  "/home/aditya/Work/IRF6_codes_figs/IRF6_epi_mes_correlations.tsv"
)

# -----------------------------
# Plotting
# -----------------------------
label_list <- c("ovary", "colorectal", "lung", "breast", "blood", "gastric")

plot_df <- result_df %>%
  mutate(
    label = ifelse(tissue %in% label_list, tissue, NA),
    color = ifelse(tissue %in% label_list, tissue, "unlabeled")
  )

ggplot(plot_df, aes(x = epi_rho, y = mes_rho)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(aes(color = color), size = 2) +
  geom_text_repel(
    aes(label = label, color = color),
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "black",
    segment.size = 0.25,
    max.overlaps = Inf
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "epithelial scores",
       y = "mesenchymal scores",
        title = "Correlation of IRF6 with Epithelial and Mesenchymal scores in CCLE
") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11),
    legend.position = "none"
  ) +
  scale_color_manual(
    values = c(
      setNames(
        rainbow(length(label_list)),
        label_list
      ),
      "unlabeled" = "grey70"
    )
  )

# TCGA --------------------------------------------------------------------

# ===============================
# Libraries
# ===============================
library(readr)
library(dplyr)
library(tidyr)
library(GSVA)
library(ggplot2)
library(ggrepel)

# ===============================
# Inputs
# ===============================
gene_name <- "IRF6"

tcga_dir <- "/home/aditya/Work/data/TCGA_datasets"

cancer_types <- c(
  "BRCA","ACC","BLCA","CESC","COAD","DLBC","ESCA","GBM","HNSC",
  "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUSC","MESO","OV",
  "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA",
  "THYM","UCEC","UCS","UVM","PAAD","LUAD","LUNG"
)

label_list <- c("BRCA","LUAD","OV","COAD","LUSC","PAAD","STAD")

# hallmark_EMT_genes must already be loaded
# source("path_to_EMT_gene_sets.R")

# ===============================
# Function: one cancer type
# ===============================
compute_correlations <- function(cancer_type, gene) {
  
  file_path <- paste0(
    tcga_dir,
    "/TCGA.",
    toupper(cancer_type),
    ".sampleMap_HiSeqV2.gz"
  )
  
  expr <- read_tsv(file_path, show_col_types = FALSE)
  expr <- as.data.frame(expr)
  rownames(expr) <- expr[,1]
  expr <- expr[,-1]
  
  # ssGSEA
  param <- GSVA::ssgseaParam(as.matrix(expr), hallmark_EMT_genes)
  nes <- GSVA::gsva(param)
  colnames(nes) <- colnames(expr)
  
  # add gene expression
  nes <- rbind(nes, expr[gene, ])
  nes <- as.matrix(nes)
  
  # correlations
  epi_test <- cor.test(
    nes["epithelial", ],
    nes[gene, ],
    method = "spearman"
  )
  
  mes_test <- cor.test(
    nes["mesenchymal", ],
    nes[gene, ],
    method = "spearman"
  )
  
  data.frame(
    cancer_type = cancer_type,
    epi_rho = epi_test$estimate,
    epi_p = epi_test$p.value,
    mes_rho = mes_test$estimate,
    mes_p = mes_test$p.value,
    stringsAsFactors = FALSE
  )
}

# ===============================
# Compute all cancers
# ===============================
results <- bind_rows(
  lapply(cancer_types, compute_correlations, gene = gene_name)
)

# ===============================
# Format table for submission
# ===============================
results_table <- results %>%
  mutate(
    epi_rho = round(epi_rho, 3),
    mes_rho = round(mes_rho, 3),
    epi_p = signif(epi_p, 3),
    mes_p = signif(mes_p, 3)
  )

write.csv(
  results_table,
  "IRF6_TCGA_EMT_correlations.csv",
  row.names = FALSE
)

# ===============================
# Prepare plotting dataframe
# ===============================
plot_df <- results_table %>%
  mutate(
    label = ifelse(cancer_type %in% label_list, cancer_type, ""),
    color = ifelse(cancer_type %in% label_list, cancer_type, "unlabeled")
  )

# ===============================
# Plot
# ===============================
Fig1c <- ggplot(plot_df, aes(x = epi_rho, y = mes_rho)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(aes(color = color), size = 2) +
  geom_text_repel(
    aes(label = label, color = color),
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "black",
    segment.size = 0.25,
    max.overlaps = Inf
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    x = "Epithelial score (Spearman ρ)",
    y = "Mesenchymal score (Spearman ρ)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "none"
  ) +
  scale_color_manual(
    values = c(
      setNames(
        rainbow(length(label_list)),
        label_list
      ),
      "unlabeled" = "grey70"
    )
  )

ggsave(
  "Fig1c_IRF6_TCGA_EMT_correlations.pdf",
  Fig1c,
  width = 5,
  height = 4,
  dpi = 300
)

