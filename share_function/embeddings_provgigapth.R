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



getEmbeddingLayer <- function(fnames, participant_ids = NULL, layer = NULL) {
  
  single_file <- FALSE
  if (is.character(fnames) && length(fnames) == 1) {
    fnames <- list(fnames)
    single_file <- TRUE
    if (is.null(participant_ids)) participant_ids <- "sample1"
  }
  
  embeddings_list <- vector("list", length(fnames))
  

  if (!is.null(participant_ids)) {
    names(embeddings_list) <- make.unique(participant_ids)
  } else {
    names(embeddings_list) <- paste0("sample", seq_along(fnames))
  }
  

  for (i in seq_along(fnames)) {
    if (file.exists(fnames[i])) {
      # Load CSV
      x <- readr::read_csv(fnames[[i]], show_col_types = FALSE)
      if (nrow(x) == 0) {
        embeddings_list[[i]] <- NA
        next
      }
      
      layer_ind <- if (is.null(layer)) {
        ncol(x)
      } else if (is.character(layer)) {
        match(layer, colnames(x))
      } else {
        layer + 1
      }
      
      # Extract tensor string
      tensor_string <- x[[1, layer_ind]]
      # Clean string and convert to numeric
      clean_string <- gsub("tensor\\(\\[\\[|\\]\\]\\)", "", tensor_string)
      embeddings_list[[i]] <- as.numeric(unlist(strsplit(clean_string, ",\\s*")))
      
    } else {
      embeddings_list[[i]] <- NA
    }
  }
  

  embeddings_list <- embeddings_list[!sapply(embeddings_list, function(x) all(is.na(x)))]


  if (single_file) return(embeddings_list[[1]])


  max_len <- max(sapply(embeddings_list, length))
  tensor_matrix <- t(sapply(embeddings_list, function(x) {
    if (length(x) < max_len) c(x, rep(NA, max_len - length(x))) else x
  }))
  
  rownames(tensor_matrix) <- names(embeddings_list)
  
  return(tensor_matrix)
}




