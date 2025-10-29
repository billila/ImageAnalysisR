# barplot nuclei hovernet output ###

library(ggplot2)
library(dplyr)
library(tidyr)


counts_total <- ov_full %>%
  select(inflammatory, benign_epithelial, necrotic, neoplastic, stromal, no_label) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(cols = everything(), names_to = "cell_type", values_to = "count")

fig <- ggplot(counts_total, aes(x = cell_type, y = count, fill = cell_type)) +
  geom_col() +
  theme_light() +
  labs(title = "Total distribution of nuclei by cell type",
       x = "Cell type",
       y = "Number of nuclei") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# pdf("/mnt/spca/run_spca_log/corr_check/heatmap_correlation_1.3M.pdf", width=16, height=9)
png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/ovarian_cancer_barplot.png", width=12, height=8, units = "in", res = 800)
fig
dev.off()


# barplot 
counts_img <- ov_full %>%
  pivot_longer(cols = c(inflammatory, benign_epithelial, necrotic, neoplastic, stromal, no_label),
               names_to = "cell_type", values_to = "count")

fig2 <- ggplot(counts_img, aes(x = factor(slide_name), y = count, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Distribution of nuclei by cell type and image",
       x = "Image",
       y = "Number of nuclei") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/ovarian_cancer_img_barplot.png", width=12, height=8, units = "in", res = 800)
fig2
dev.off()



counts_img_prop <- counts_img %>%
  group_by(slide_name) %>%
  mutate(prop = count / sum(count))

fig3 <- ggplot(counts_img_prop, aes(x = factor(slide_name), y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of nuclei by cell type and image",
       x = "Image",
       y = "Proportion of nuclei") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/ovarian_cancer_img_prop_barplot.png", width=12, height=8, units = "in", res = 800)
fig3
dev.off()



counts_img_prop <- counts_img %>%
  group_by(slide_name) %>%
  mutate(prop = count / sum(count))

counts_img_prop <- counts_img_prop %>%
  mutate(cell_type = factor(cell_type,
                            levels = c("neoplastic",     
                                       "benign_epithelial",
                                       "inflammatory",
                                       "necrotic",
                                       "stromal",
                                       "no_label")))


order_neoplastic <- counts_img_prop %>%
  filter(cell_type == "neoplastic") %>%
  arrange(desc(prop)) %>%
  pull(slide_name)

counts_img_prop <- counts_img_prop %>%
  mutate(slide_name = factor(slide_name, levels = order_neoplastic))

# 4. Plot
fig3 <- ggplot(counts_img_prop,
               aes(x = slide_name, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of nuclei by cell type and image",
       x = "Image",
       y = "Proportion of nuclei") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fig3

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/ovarian_cancer_img_prop_barplot_neoplastic_order.png", width=12, height=8, units = "in", res = 800)
fig3
dev.off()

