library(ConsensusTME)
library(SummarizedExperiment)

#load("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/se_OV_data.RData")

exprMatrix <- se_ov@assays@data@listData[["counts"]]
rownames(exprMatrix) <- gene_info$gene_name
colnames(exprMatrix) <- colnames(se_ov)

# xcell <- xCellAnalysis(exprMatrix, 
#                        signatures = signatures_UNIVOQUE_LauraRev_GeneSetCollection,
#                        cell.types.use = c("B_cells", "Dendritic_cells", 
#                                           "NK_cells", "T_cells", "Endothelial", 
#                                           "Mesenchymal", "CAFs", "Stromal_cells", 
#                                           "Cancer_cells", "MacroMono", "Granulocytes_EosNeu"))

library(GSVA)
library(tidyverse)
colnames(exprMatrix) <- make.unique(colnames(exprMatrix))
exprMatrix_clean <- exprMatrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  # ripara eventuali nomi duplicati delle colonne
  as_tibble(.name_repair = "unique") %>%
  group_by(gene) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  as.data.frame()

rownames(exprMatrix_clean) <- exprMatrix_clean$gene
exprMatrix_clean$gene <- NULL

exprMatrix_clean <- as.matrix(exprMatrix_clean)

####### consensus TME ##################
#bulkExpMatrix <- as.matrix(read.delim(bulkGeneExpression.txt, row.names = 1))

score_TME <- ConsensusTME::consensusTMEAnalysis(exprMatrix_clean, cancer = "OV", statMethod = "ssgsea")

# Producing ConsensusTME Estimates Using The Following Parameters:
#   Statistical Framework: "ssgsea"

# consensuTME
# Gene Sets For Cancer Type: "OV"
# Sample Size: 82
# ℹ GSVA version 2.2.0
# ! 4506 genes with constant values throughout the samples
# ℹ Calculating  ssGSEA scores for 19 gene sets
# ℹ Calculating ranks
# ℹ Calculating rank weights
# ℹ Normalizing ssGSEA scores
# ✔ Calculations finished

# xcell con signature laura 
# ℹ GSVA version 2.2.0
# ! 4506 genes with constant values throughout the samples
# ! No annotation metadata available in the input expression data object
# ! Attempting to directly match identifiers in expression data to gene sets
# ℹ Calculating  ssGSEA scores for 11 gene sets
# ℹ Calculating ranks
# ℹ Calculating rank weights
# ✔ Calculations finished


xcell <- score_TME


#### stackplot all cell types ######
xcell_t <- as.data.frame(t(xcell))
xcell_norm <- xcell_t / rowSums(xcell_t)

library(tidyr)
library(dplyr)

xcell_long <- xcell_norm %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "Proportion")

xcell_long$Sample <- substr(xcell_long$Sample, 1, 12)

library(ggplot2)

xcell_long <- xcell_long %>%
  group_by(Sample, CellType) %>%
  summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop")

ggplot(xcell_long, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +  # oppure `element_text(angle = 90)`
  labs(title = "consensusTME", x = "Samples", y = "Proportion") +
  guides(fill = guide_legend(ncol = 2))

library(pheatmap)

pheatmap(t(xcell), 
         scale = "row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Immune/stromal cell abundance (xCell)")


xcell_long_cluster <- xcell_long %>%
  left_join(ov_full %>% dplyr::select(patientID, cluster), 
            by = c("Sample" = "patientID"))

xcell_long_cluster <- xcell_long_cluster %>%
  group_by(Sample, CellType, cluster) %>%
  summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop")

TME_clust <- ggplot(xcell_long_cluster, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  facet_grid(~ cluster, scales = "free_x", space = "free_x") +
  labs(
    title = "consensusTME score per cluster",
    x = "Samples",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(ncol = 2))

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/consensusTME_clust.png", width=12, height=8, units = "in", res = 800)
TME_clust
dev.off()

boxplot_TME <- ggplot(xcell_long_cluster, aes(x = cluster, y = Proportion, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ CellType, scales = "free_y") +
  theme_minimal() +
  labs(title = "consensusTME scores across clusters",
       x = "Cluster", y = "Proportion")

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/TME_boxplot_celltype.png", width=12, height=8, units = "in", res = 800)
boxplot_TME
dev.off()

boxplot_TME_clust <- ggplot(xcell_long_cluster, aes(x = CellType, y = Proportion, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_minimal() +
  labs(title = "consensusTME scores across clusters",
       x = "Cluster", y = "Proportion")

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/TME_boxplot_clust.png", width=12, height=8, units = "in", res = 800)
boxplot_TME_clust
dev.off()


library(dplyr)

kruskal_results <- xcell_long_cluster %>%
  group_by(CellType) %>%
  summarise(
    kruskal = list(kruskal.test(Proportion ~ cluster)),
    .groups = "drop"
  ) %>%
  mutate(
    p.value = sapply(kruskal, function(x) x$p.value)
  )

kruskal_results %>% dplyr::select(CellType, p.value)


library(FSA)   # per dunnTest

dunn_results <- xcell_long_cluster %>%
  group_by(CellType) %>%
  do({
    test <- FSA::dunnTest(Proportion ~ cluster, data = ., method = "bh")
    data.frame(test$res)
  })

dunn_results %>%
  filter(P.adj < 0.05)

