# clustree
# choose the number of cluster 

library(factoextra)
library(cluster)
library(NbClust)

library(dplyr)
library(clustree)

set.seed(123)

# Subset embeddings
embeddings <- ov_full[, embedding_cols]


for (k in 2:6) {
  km <- kmeans(embeddings, centers = k, nstart = 25)
  ov_full[[paste0("cluster_k", k)]] <- factor(km$cluster)
}


clustree <- clustree(ov_full, prefix = "cluster_k") +
  ggplot2::labs(title = "Clustree for k-means clustering")

png("/mnt/tcga_images/purity_ploidy_OV/case_study_OV/plot/clustree_kmeans.png", width=12, height=8, units = "in", res = 800)
clustree
dev.off()
