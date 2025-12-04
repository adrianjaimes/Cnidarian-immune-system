#######################################Correlations_SingleCell (Figure 2C,D,E,F)#######################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)

seurat_obj <- readRDS("seurat_clust.RDS")

# 1. Load the metacell mapping
cell_to_mc <- read.table("CellID_to_MC_table.txt", header = FALSE, sep = "\t", col.names = c("CellID", "Metacell"))

# 2. Add Metacell info to Seurat object metadata
seurat_obj@meta.data$CellID <- rownames(seurat_obj@meta.data)  # ensure CellID is correctly set
seurat_obj@meta.data <- left_join(seurat_obj@meta.data, cell_to_mc, by = "CellID")


# 1. Define gene IDs
rlrb_gene <- "Nvec-vc1.1-XM-048731786.1"
cardib_gene <- "Nvec-vc1.1-XM-048731785.1"

# 2. Extract normalized expression data (genes x cells)
expr_mat <- GetAssayData(seurat_obj, layer = "data")

# 3. Extract and format expression of genes of interest (cells x genes)
expr_df <- as.data.frame(t(as.matrix(expr_mat[c(rlrb_gene, cardib_gene), ])))
colnames(expr_df) <- c("RLRb", "CARDIB")
expr_df$CellID <- rownames(expr_df)

# 4. Get metadata (must include Metacell and Condition)
meta_df <- seurat_obj@meta.data %>%
  select(CellID, Metacell, Condition)


# 5. Merge expression and metadata
merged_df <- inner_join(meta_df, expr_df, by = "CellID")

# 6. Average per metacell × condition
meta_summary <- merged_df %>%
  group_by(Metacell, Condition) %>%
  summarise(
    RLRb = mean(RLRb),
    CARDIB = mean(CARDIB),
    .groups = "drop"
  )

condition_colors <- c("Ctrl" = "#0072B2", "iHCl" = "#D55E00", "tPIC" = "#009E73")


plot_correlation <- function(df, title, condition_colors) {
  cor_test <- cor.test(df$RLRb, df$CARDIB, method = "pearson")
  R_val <- round(cor_test$estimate, 2)
  p_val <- signif(cor_test$p.value, 2)
  label_text <- paste0("R = ", R_val, "\np < ", p_val)
  
  ggplot(df, aes(x = RLRb, y = CARDIB, color = Condition)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", fill = "gray70", alpha = 0.3) +
    scale_color_manual(values = condition_colors, drop = FALSE) +
    annotate("text",
             x = max(df$RLRb, na.rm = TRUE) * 0.95,
             y = min(df$CARDIB, na.rm = TRUE) * 1.05,
             label = label_text,
             hjust = 1, vjust = 0,
             size = 5) +
    labs(
      title = title,
      x = "RLRb UMI fraction",
      y = "CARDIB UMI fraction",
      color = "Condition"
    ) +
    theme_minimal(base_size = 14)
}



p_all  <- plot_correlation(meta_summary, title = "All Conditions Combined", condition_colors = condition_colors)

p_ctrl <- plot_correlation(filter(meta_summary, Condition == "Ctrl"),  "Ctrl",  condition_colors)
p_ihcl <- plot_correlation(filter(meta_summary, Condition == "iHCl"),  "iHCl",  condition_colors)
p_tpic <- plot_correlation(filter(meta_summary, Condition == "tPIC"),  "tPIC",  condition_colors)



# View
p_all
p_ctrl
p_ihcl
p_tpic

pdf("RLRb_vs_CARDIB_colored_by_condition_consistent.pdf", width = 6, height = 5)
print(p_all)
print(p_ctrl)
print(p_ihcl)
print(p_tpic)
dev.off()

write.csv(meta_summary, "meta_summary.csv")