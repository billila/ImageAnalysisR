# ESTIMATE tumor purity with tidyestimate
# install.packages("tidyestimate")
library(tidyestimate)

expr_mat <- assay(se_ov, "counts")
rownames(expr_mat) <- gene_info$gene_name
gene <- gene_info$gene_name
expr_mat <- as.data.frame(expr_mat)
expr_mat <- expr_mat[, !duplicated(colnames(expr_mat))]


expr_df <- as.data.frame(expr_mat) %>%
  rownames_to_column("gene_symbol") %>%
  as_tibble()

expr_df_clean <- expr_df %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), mean), .groups = "drop")

expr_mat_clean <- as.matrix(expr_df_clean[,-1])
rownames(expr_mat_clean) <- expr_df_clean$gene_symbol


scores <- expr_mat_clean |> 
  tidyestimate::filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE) |> 
  tidyestimate::estimate_score(is_affymetrix = TRUE)



scores$patientID <- scores$sample
estimate_purity <- merge(
  ov_full,
  scores,
  by = "patientID",
  all.x = TRUE
)

cor(estimate_purity$ESTIMATE, estimate_purity$purity, 
    use = "pairwise.complete.obs")

plot(estimate_purity$ESTIMATE, estimate_purity$purity)

cor(estimate_purity$stromal.y, estimate_purity$strom_perc,
    use = "pairwise.complete.obs")

plot(estimate_purity$stromal.y, estimate_purity$strom_perc)

estimate_purity$stromal.y_scaled <- (estimate_purity$stromal.y - min(estimate_purity$stromal.y, na.rm = TRUE)) /
  (max(estimate_purity$stromal.y, na.rm = TRUE) - min(estimate_purity$stromal.y, na.rm = TRUE))

summary(estimate_purity$stromal.y_scaled)
range(estimate_purity$stromal.y_scaled, na.rm = TRUE)

cor(estimate_purity$stromal.y_scaled, estimate_purity$strom_perc,
    use = "pairwise.complete.obs")

cor.test(estimate_purity$stromal.y_scaled, estimate_purity$strom_perc,
         use = "pairwise.complete.obs")

plot(estimate_purity$stromal.y_scaled, estimate_purity$strom_perc)


library(ggplot2)
library(ggpubr) # per stat_cor

# Scatter plot
cor_plot <- ggplot(estimate_purity, aes(x = strom_perc, y = stromal.y_scaled)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  stat_cor(method = "pearson", label.x = 0.55, label.y = 0.9) + # mostra r e p-value
  theme_minimal() +
  labs(
    x = "Stromal cells (%) from Hover-Net segmentation",
    y = "Scaled ESTIMATE stromal score",
    title = "Correlation between ESTIMATE stromal score and H&E stromal fraction (Hover-Net)"
  )

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/estimate_stromal_corr.png", width=12, height=8, units = "in", res = 800)
cor_plot
dev.off()
