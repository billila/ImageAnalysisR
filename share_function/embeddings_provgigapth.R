library(tidyverse)


#input_folder <- "/mnt/tcga_images/prov_gigapath/output/TCGA_OV/"
input_folder <- "/home/ilaria/Documents/provgigapath_leonardo/TCGA_OV/"
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

process_single_file <- function(file_path) {
  df <- read_csv(file_path)
  print(file_path)
  embedding <- df$last_layer_embed[1] %>%
    str_remove_all("tensor\\(\\[\\[|\\]\\]\\)") %>%
    str_remove_all("\n") %>%
    str_squish() %>%
    str_split(",\\s*") %>%
    unlist() %>%
    as.numeric()
  
  tibble(
    slide_name = df$slide_name[1],
    tumor_type = str_extract(input_folder, "TCGA[_-][A-Z]+"),
    !!!set_names(as.list(embedding), paste0("V", seq_along(embedding)))  # Colonne V1, V2, V3, ...
  )
  
}


all_embeddings <- map_dfr(csv_files, process_single_file)

write_csv(all_embeddings, "/home/ilaria/Documents/provgigapath_leonardo/emb_gigapth_OV_last_embed.csv")

###############################################################################








