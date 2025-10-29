#######script to reprocduce plot and test between  ########
####### clustering - purity - proprtion- nuclei type ########



# clustering embeddings 
library(survival)
library(survminer)
ov_full <- ov_full[!is.na(ov_full$V1), ]
embedding_cols <- grep("^V", names(ov_full), value = TRUE)
set.seed(123)
km <- kmeans(ov_full[, embedding_cols], centers = 4)
ov_full$cluster <- factor(km$cluster)
table(ov_full$cluster)

ov_full <- ov_full %>%
  mutate(
    time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup),
    event = ifelse(!is.na(days_to_death), 1, 0)
  )


fit <- survfit(Surv(time, vital_status) ~ cluster, data = ov_full)
km_plot <- ggsurvplot(fit, data = ov_full, risk.table = TRUE, pval = TRUE, xlim = c(0,3000))
pairwise_survdiff(Surv(time, vital_status) ~ cluster, data = ov_full)

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/km_clusters.png", width=12, height=8, units = "in", res = 800)
km_plot
dev.off()

library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)


# umap clustering

library(uwot)
ov_full <- ov_full[!is.na(ov_full$V1), ]
umap_res <- umap(ov_full[, grepl("^V", colnames(ov_full))])
plot(umap_res, col = as.factor(ov_full$cluster), pch = 19)

library(uwot)
library(ggplot2)
library(cowplot)   

library(Rtsne)
library(uwot)
library(ggplot2)
library(cowplot)

#  embeddings columns 
embedding_cols <- grep("^V", colnames(ov_full), value = TRUE)
X <- ov_full[, embedding_cols]

# t-SNE con Rtsne
set.seed(123)
tsne_out <- Rtsne(X, perplexity = 30, check_duplicates = FALSE)
tsne_df <- as.data.frame(tsne_out$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df$cluster <- as.factor(ov_full$cluster)

p_tsne <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 3.0, alpha = 0.9) +
  labs(title = "t-SNE (k-means clusters)") +
  #scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# UMAP uwot
set.seed(123)
umap_out <- umap(X)
umap_df <- as.data.frame(umap_out)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$cluster <- as.factor(ov_full$cluster)

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 3.0, alpha = 0.9) +
  labs(title = "UMAP (k-means clusters)") +
  #scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Appaiati
umap_tsne <- plot_grid(p_tsne, p_umap, labels = c("A", "B"), ncol = 2, align = "h")

p_tsne_noleg <- p_tsne + theme(legend.position = "none")
p_umap_noleg <- p_umap + theme(legend.position = "none")


legend_shared <- get_legend(
  p_tsne + theme(legend.position = "bottom")  # prendi la legenda in basso
)

# no legend
plots <- plot_grid(p_tsne_noleg, p_umap_noleg,
                   labels = c("A", "B"), ncol = 2, align = "h")

# add legend 
umap_tsne <- plot_grid(plots, legend_shared, ncol = 1, rel_heights = c(1, 0.1))
png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/tsne_umap_clusters.png", width=12, height=8, units = "in", res = 800)
umap_tsne
dev.off()


ov_full_prop <- ov_full %>%
  mutate(across(c(inflammatory, benign_epithelial, necrotic, 
                  neoplastic, stromal, no_label),
                ~ .x / num_nuclei, .names = "prop_{.col}"))

# long format per plotting
ov_full_long <- ov_full_prop %>%
  select(cluster, starts_with("prop_")) %>%
  pivot_longer(cols = starts_with("prop_"),
               names_to = "cell_type", values_to = "proportion") %>%
  mutate(cell_type = gsub("prop_", "", cell_type))

#  Boxplot proprortion cluster
boxplot_cluster <- ggplot(ov_full_long, aes(x = cluster, y = proportion, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~cell_type, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of cell type proportions across clusters",
       x = "Cluster",
       y = "Proportion of nuclei")

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/prop_clusters_boxplot.png", width=12, height=8, units = "in", res = 800)
boxplot_cluster
dev.off()

summary_table <- ov_full_long %>%
  group_by(cluster, cell_type) %>%
  summarise(mean_proportion = mean(proportion, na.rm = TRUE)) %>%
  pivot_wider(names_from = cell_type, values_from = mean_proportion)


summary_table

# 2. Barplot mean prorportion cluster
cluster_summary <- ov_full_long %>%
  group_by(cluster, cell_type) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop")

ggplot(cluster_summary, aes(x = cluster, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Average cell type composition per cluster",
       x = "Cluster",
       y = "Average proportion of nuclei",
       fill = "Cell type")



###################################### 


ov_full_prop <- ov_full %>%
  mutate(across(c(inflammatory, benign_epithelial, necrotic, 
                  neoplastic, stromal, no_label),
                ~ .x / num_nuclei, .names = "prop_{.col}"))


ov_full_long <- ov_full_prop %>%
  select(cluster, starts_with("prop_")) %>%
  pivot_longer(cols = starts_with("prop_"),
               names_to = "cell_type", values_to = "proportion") %>%
  mutate(cell_type = gsub("prop_", "", cell_type))

#  (Kruskal-Wallis) 
stat_tests <- ov_full_long %>%
  group_by(cell_type) %>%
  kruskal_test(proportion ~ cluster) %>%
  ungroup()

print(stat_tests)



library(dplyr)
library(rstatix)

# Dunn's test post-hoc per ogni cell type
posthoc_tests <- ov_full_long %>%
  group_by(cell_type) %>%
  dunn_test(proportion ~ cluster, p.adjust.method = "BH") %>%  # BH = FDR correction
  ungroup()


print(posthoc_tests)

significant_pairs <- posthoc_tests %>%
  filter(p.adj < 0.05)

print(significant_pairs)

# Cluster 1: immune-rich, con molti nuclei no-label, mix di neoplastici.
# 
# Cluster 2: necrotic + neoplastic, immune-poor.
# 
# Cluster 3: immune-rich, meno necrotico, meno no-label.
# 
# Cluster 4: stromal-rich, pochi neoplastici.

# Violin plot  cluster e tipi cellulari
ggplot(ov_full_long, aes(x = cluster, y = proportion, fill = cell_type)) +
  geom_violin(position = position_dodge(width = 0.8), alpha = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Cell type proportions across clusters",
       x = "Cluster",
       y = "Proportion of nuclei",
       fill = "Cell type")



### purity ##########

purity <- c("purity_hovernet","tcga_purity_mean", "ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE" )
library(dplyr)
library(ggpubr)
library(rstatix)

# Test Kruskal-Wallis  purity
kruskal_purity <- kruskal_test(purity_hovernet ~ cluster, data = ov_full)
print(kruskal_purity)

# Post-hoc Dunn test  
dunn_purity <- dunn_test(purity_hovernet ~ cluster, data = ov_full, p.adjust.method = "BH")
print(dunn_purity)

# Boxplot pval
p_purity <- ggboxplot(ov_full, x = "cluster", y = "purity_hovernet", fill = "cluster",
                      palette = "Set2", add = "jitter") +
  stat_compare_means(method = "kruskal.test", label.y = 1.05) + # test globale
  stat_compare_means(method = "dunn.test", label = "p.signif", 
                     ref.group = "1", hide.ns = TRUE) + # confronti post-hoc vs C1
  labs(title = "Tumor purity across clusters",
       x = "Cluster",
       y = "Purity") +
  theme_bw()

p_purity




#######


library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)


purity_vars <- c("purity_hovernet","tcga_purity_mean","ESTIMATE",
                 "ABSOLUTE","LUMP","IHC","CPE")

# Reshape long
ov_purity_long <- ov_full %>%
  select(cluster, all_of(purity_vars)) %>%
  pivot_longer(cols = all_of(purity_vars),
               names_to = "purity_metric",
               values_to = "purity_value")

# global test (Kruskal-Wallis)
kw_tests <- ov_purity_long %>%
  group_by(purity_metric) %>%
  kruskal_test(purity_value ~ cluster)
print(kw_tests)

# Post-hoc Dunn test con correzione BH
posthoc_tests <- ov_purity_long %>%
  group_by(purity_metric) %>%
  dunn_test(purity_value ~ cluster, p.adjust.method = "BH")
print(posthoc_tests)

# Sava table
stat_tests_purity <- list(
  kruskal = kw_tests,
  posthoc = posthoc_tests
)

significant_pairs <- posthoc_tests %>%
  filter(p.adj < 0.05)

print(significant_pairs)

# Boxplot con facet_wrap
p <- ggboxplot(
  ov_purity_long, x = "cluster", y = "purity_value",
  fill = "cluster", add = "jitter"
) +
  #stat_compare_means(method = "kruskal.test", label.y = 1.05) +
  labs(
    title = "Tumor purity across clusters",
    x = "Cluster", y = "Purity"
  ) +
  facet_wrap(~purity_metric, scales = "free_y") +
  theme_minimal()

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/purity_clusters_boxplot.png", width=12, height=8, units = "in", res = 800)
p
dev.off()


print(stat_tests_purity$kruskal)
print(stat_tests_purity$posthoc)


p


#####  NA purity 


library(dplyr)
library(ggplot2)

purity_vars <- c("purity_hovernet","tcga_purity_mean","ESTIMATE",
                 "ABSOLUTE","LUMP","IHC","CPE")

#  %  NA
na_summary <- ov_full %>%
  summarise(across(all_of(purity_vars),
                   ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(cols = everything(),
               names_to = "purity_metric",
               values_to = "na_percent")

# Barplot delle % di NA
p_na <- ggplot(na_summary, aes(x = reorder(purity_metric, -na_percent),
                               y = na_percent, fill = purity_metric)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", na_percent)),
            vjust = -0.5, size = 3.5) +
  labs(
    title = "Percentage of missing values (NA) for purity",
    x = "Purity metric", y = "% NA"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/purity_na_barplot.png", width=12, height=8, units = "in", res = 800)
p_na
dev.off()

na_matrix <- ov_full %>%
  select(all_of(purity_vars)) %>%
  mutate(across(everything(), is.na))

library(ComplexUpset)

na_matrix <- ov_full %>%
  select(all_of(purity_vars)) %>%
  mutate(across(everything(), is.na)) %>%
  mutate(across(everything(), as.integer))

# UpSet plot
upset(
  na_matrix,
  intersect = purity_vars,
  name = "Missingness",
  width_ratio = 0.2,
  base_annotations = list(
    'Intersection size' = intersection_size(
      text_mapping=aes(label=after_stat(count)), # <-- fix
      text = list(size=3, vjust=-0.5)
    ),
    'Set size' = upset_set_size()
  )
)


library(tidyr)
library(ComplexUpset)


set_matrix <- ov_full %>%
  select(patientID, all_of(purity_vars)) %>%
  mutate(across(all_of(purity_vars), ~ !is.na(.))) %>%   
  mutate(across(all_of(purity_vars), as.integer))        


set_matrix_bin <- set_matrix %>% select(-patientID)

# UpSet plot
upset <- upset(
  set_matrix_bin,
  intersect = purity_vars,
  name = "Patients",
  width_ratio = 0.2,
  base_annotations = list(
    'Intersection size' = intersection_size()
  )
)

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/purity_na_upsetplot.png", width=12, height=8, units = "in", res = 800)
upset
dev.off()
