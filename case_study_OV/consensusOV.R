# consensusOV 

library(consensusOV)
library(SummarizedExperiment)
# load("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/se_OV_data.RData")

se_ov

# i need entrez id gene name 
library(org.Hs.eg.db)
library(AnnotationDbi)
columns(org.Hs.eg.db)

expr_counts <- assay(se_ov, "counts")

gene_ids <- gsub("\\..*", "", rownames(expr_counts))

expr_counts_collapsed <- rowsum(expr_counts, group = gene_ids)

#expr_counts_collapsed <- rowsum(expr_counts_collapsed, group = gene_ids)

ids <- mapIds(org.Hs.eg.db, rownames(expr_counts_collapsed), "ENTREZID", "ENSEMBL")

dim(expr_counts_collapsed) 

######### consensus #################################
consensus <- get.subtypes(as.matrix(expr_counts_collapsed), ids)

table(consensus$consensusOV.subtypes)

# ora che ho questi risultati che me ne faccio?
# IMR_consensus DIF_consensus PRO_consensus MES_consensus 
# 20            26            20            16 

# > table(consensus$consensusOV.subtypes)
# 
# IMR_consensus DIF_consensus PRO_consensus MES_consensus 
# 19            27            19            17 



# unirli ad db originale e confrontarli con i 4 gruppi del clustering degli embeddings 

sum(duplicated(ov_full$patientID))
ov_full$patientID[duplicated(ov_full$patientID)]

sum(duplicated(rownames(consensus[[2]])))
rownames(consensus[[2]])[duplicated(rownames(consensus[[2]]))]

cons_df <- consensus[[2]][!duplicated(rownames(consensus[[2]])), ]
cons_df <- data.frame(patientID = rownames(cons_df), cons_df)

cons_df$consensuOV <- consensus[[1]][!duplicated(rownames(consensus[[2]]))]
merged_OV <- merge(
  ov_full,
  cons_df,
  by = "patientID",
  all.x = TRUE
)

tab <- table(merged_OV$cluster, merged_OV$consensuOV)
tab
prop.table(tab, 1)
prop.table(tab, 2)
prop.table(tab)

addmargins(tab)

library(janitor)
tabyl(merged_OV, cluster, consensuOV) %>% adorn_percentages("row")

df_tab <- as.data.frame(tab)

cont_table <- ggplot(df_tab, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "cluster",
    y = "consensuOV",
    fill = "Counts",
    title = "Contingency Table Heatmap"
  ) +
  theme_minimal()

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/cont_table_consensuov.png", width=12, height=8, units = "in", res = 800)
cont_table
dev.off()

prop_tab <- prop.table(tab, margin = 1)

# convert to long data.frame for ggplot
df_tab <- as.data.frame(prop_tab)

ggplot(df_tab, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq, accuracy = 1)), color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "consensuOV",
    y = "cluster",
    fill = "Row %", 
    title = "Row-wise Percentage Heatmap"
  ) +
  theme_minimal()


######################################### prob #####################

library(dplyr)
library(tidyr)
library(ggplot2)

# put your real probability column names here:
prob_cols <- c("PRO_consensus", "IMR_consensus", "DIF_consensus", "MES_consensus")

# reshape into long format
merged_long <- merged_OV %>%
  dplyr::select(cluster, all_of(prob_cols)) %>%
  pivot_longer(cols = all_of(prob_cols),
               names_to = "subtype",
               values_to = "probability") %>%
  mutate(subtype = gsub("_prob", "", subtype))  # clean names

# boxplots of assignment probabilities by cluster
bp_cluster <- ggplot(merged_long, aes(x = cluster, y = probability, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~subtype, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of assignment consensusOV probabilities across clusters",
       x = "Cluster",
       y = "Probability")

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/boxplot_cluster_consensuov.png", width=12, height=8, units = "in", res = 800)
bp_cluster
dev.off()

# ggplot(merged_long, aes(x = cluster, y = probability, fill = cluster)) +
#   geom_violin(trim = FALSE, alpha = 0.6) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
#   facet_wrap(~subtype, scales = "free_y") +
#   theme_minimal() +
#   labs(title = "Assignment probability distributions across clusters",
#        x = "Cluster",
#        y = "Probability")

summary_table <- merged_long %>%
  group_by(cluster, subtype) %>%
  summarise(
    mean_prob = mean(probability, na.rm = TRUE),
    sd_prob   = sd(probability, na.rm = TRUE),
    n         = n(),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = sprintf("%.2f ± %.2f", mean_prob, sd_prob))

summary_table

library(rstatix)

stat_tests <- merged_long %>%
  group_by(subtype) %>%
  kruskal_test(probability ~ cluster) %>%
  ungroup()

print(stat_tests)

posthoc_tests <- merged_long %>%
  group_by(subtype) %>%
  dunn_test(probability ~ cluster, p.adjust.method = "BH") %>%
  ungroup()

print(posthoc_tests %>%
        filter(p.adj < 0.05))

######### consensus 2  #################################
# get.konecny.subtypes(as.matrix(expr_counts_clean), entrez.ids = ids_clean)
# table(consensus$consensusOV.subtypes)
# get.verhaak.subtypes(as.matrix(expr_counts_clean), entrez.ids = ids_clean)


valid_idx <- !is.na(ids)
expr_counts_clean <- expr_counts_collapsed[valid_idx, ]
ids_clean <- ids[valid_idx]
ids_clean <- make.unique(ifelse(is.na(ids), "NA", as.character(ids)))
expr_counts_clean <- expr_counts_collapsed
data(sigOvcAngiogenic)
get.bentink.subtypes(as.matrix(expr_counts_clean), entrez.ids = ids_clean)
get.helland.subtypes(as.matrix(expr_counts_clean), entrez.ids = ids_clean)

consensus <- get.subtypes(as.matrix(expr_counts_clean), ids_clean)

table(consensus$consensusOV.subtypes)

# IMR_consensus DIF_consensus PRO_consensus MES_consensus 
# 20            26            20            16 





sum(duplicated(ov_full$patientID))
ov_full$patientID[duplicated(ov_full$patientID)]

sum(duplicated(rownames(consensus[[2]])))
rownames(consensus[[2]])[duplicated(rownames(consensus[[2]]))]

cons_df <- consensus[[2]][!duplicated(rownames(consensus[[2]])), ]
cons_df <- data.frame(patientID = rownames(cons_df), cons_df)

cons_df$consensuOV <- consensus[[1]][!duplicated(rownames(consensus[[2]]))]
merged_OV <- merge(
  ov_full,
  cons_df,
  by = "patientID",
  all.x = TRUE
)

tab <- table(merged_OV$cluster, merged_OV$consensuOV)
tab
prop.table(tab, 1)
prop.table(tab, 2)
prop.table(tab)

addmargins(tab)

library(janitor)
tabyl(merged_OV, cluster, consensuOV) %>% adorn_percentages("row")

df_tab <- as.data.frame(tab)

ggplot(df_tab, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "cluster",
    y = "consensuOV",
    fill = "Counts",
    title = "Contingency Table Heatmap"
  ) +
  theme_minimal()


prop_tab <- prop.table(tab, margin = 1)

# convert to long data.frame for ggplot
df_tab <- as.data.frame(prop_tab)

ggplot(df_tab, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq, accuracy = 1)), color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "consensuOV",
    y = "cluster",
    fill = "Row %", 
    title = "Row-wise Percentage Heatmap"
  ) +
  theme_minimal()


######################################### prob #####################

library(dplyr)
library(tidyr)
library(ggplot2)

# put your real probability column names here:
prob_cols <- c("PRO_consensus", "IMR_consensus", "DIF_consensus", "MES_consensus")

# reshape into long format
merged_long <- merged_OV %>%
  select(cluster, all_of(prob_cols)) %>%
  pivot_longer(cols = all_of(prob_cols),
               names_to = "subtype",
               values_to = "probability") %>%
  mutate(subtype = gsub("_prob", "", subtype))  # clean names

# boxplots of assignment probabilities by cluster
ggplot(merged_long, aes(x = cluster, y = probability, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~subtype, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of assignment consensusOV probabilities across clusters",
       x = "Cluster",
       y = "Probability")

# ggplot(merged_long, aes(x = cluster, y = probability, fill = cluster)) +
#   geom_violin(trim = FALSE, alpha = 0.6) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
#   facet_wrap(~subtype, scales = "free_y") +
#   theme_minimal() +
#   labs(title = "Assignment probability distributions across clusters",
#        x = "Cluster",
#        y = "Probability")

summary_table <- merged_long %>%
  group_by(cluster, subtype) %>%
  summarise(
    mean_prob = mean(probability, na.rm = TRUE),
    sd_prob   = sd(probability, na.rm = TRUE),
    n         = n(),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = sprintf("%.2f ± %.2f", mean_prob, sd_prob))

summary_table

library(rstatix)

stat_tests <- merged_long %>%
  group_by(subtype) %>%
  kruskal_test(probability ~ cluster) %>%
  ungroup()

print(stat_tests)

posthoc_tests <- merged_long %>%
  group_by(subtype) %>%
  dunn_test(probability ~ cluster, p.adjust.method = "BH") %>%
  ungroup()

print(posthoc_tests)
