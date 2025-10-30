json_to_SpatialExperiment <- function(json_path, include_contours = FALSE) {
  # Load required packages
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Please install the 'jsonlite' package.")
  if (!requireNamespace("SpatialExperiment", quietly = TRUE))
    stop("Please install the 'SpatialExperiment' package.")
  
  library(jsonlite)
  library(SpatialExperiment)
  
  # Read JSON
  message("Reading JSON file: ", json_path)
  data <- fromJSON(json_path)
  if (is.null(data$nuc)) stop("No 'nuc' field found in JSON.")
  nuclei <- data$nuc
  
  # Type mapping
  type_map <- data.frame(
    type = 0:5,
    label = c("nolabe", "neopla", "inflam", "connec", "necros", "no-neo"),
    R = c(0, 255, 0, 0, 255, 255),
    G = c(0, 0, 255, 0, 255, 165),
    B = c(0, 0, 0, 255, 0, 0)
  )
  
  # Extract nuclei info
  message("Extracting nuclei data...")
  cells <- lapply(names(nuclei), function(id) {
    n <- nuclei[[id]]
    data.frame(
      cell_id = id,
      x = n$centroid[1],
      y = n$centroid[2],
      type = n$type,
      type_prob = n$type_prob,
      stringsAsFactors = FALSE
    )
  }) |> do.call(rbind, args = _)
  
  # Join with labels/colors
  cells <- merge(cells, type_map, by = "type", all.x = TRUE)
  
  # Build SpatialExperiment
  assay_data <- matrix(0, nrow = 0, ncol = nrow(cells))
  
  message("Creating SpatialExperiment object...")
  spe <- SpatialExperiment(
    assays = list(counts = assay_data),
    colData = cells,
    spatialCoords = as.matrix(cells[, c("x", "y")])
  )
  
  if (include_contours) {
    message("Adding contour data to metadata...")
    contour_list <- lapply(nuclei, function(n) n$contour)
    metadata(spe)$contours <- contour_list
  }
  
  metadata(spe)$type_map <- type_map
  
  message("Done! Created SpatialExperiment with ", ncol(spe), " cells.")
  return(spe)
}


spe <- json_to_SpatialExperiment(
  json_path = "/home/ilaria/Desktop/sid_he/output/json/day_11_liver_sid_21.json",
  include_contours = TRUE
)

head(colData(spe))
metadata(spe)$type_map


sfe <- SpatialFeatureExperiment::toSpatialFeatureExperiment(spe)


gg <- data.frame(spatialCoords(spe), colData(spe))
pal <- RColorBrewer::brewer.pal(length(unique(gg$label)), "Paired")

pal <- c(
  connec = "#8DD3C7",   # verde acqua
  inflam = "#FB8072",   # rosso salmone
  necros = "#BEBADA",   # viola chiaro
  neopla = "#80B1D3",   # blu medio
  `no-neo` = "#FDB462", # arancio
  nolabe = "#B3B3B3"    # grigio neutro
)

pal <- c(
  connec = "#3B9AB2",   # blu-verde (stromale)
  inflam = "#F21A00",   # rosso (infiammazione)
  necros = "#EDE530",   # giallo (necrosi)
  neopla = "#7E1E9C",   # viola (neoplastico)
  `no-neo` = "#999999", # grigio (non neoplastico)
  nolabe = "#CCCCCC"    # grigio chiaro (non etichettato)
)
pal <- c(
  connec = "#1b9e77",
  inflam = "#d95f02",
  necros = "#7570b3",
  neopla = "#e7298a",
  `no-neo` = "#66a61e",
  nolabe = "#999999"
)
pal <- c(
  connec = "#3B9AB2",   # blu-verde (stromale)
  inflam = "#F21A00",   # rosso (infiammazione)
  necros = "black",   # giallo (necrosi)
  neopla = "#7E1E9C",   # viola (neoplastico)
  `no-neo` = "#999999", # grigio (non neoplastico)
  nolabe = "#CCCCCC"    # grigio chiaro (non etichettato)
)

p_feat <- ggplot(gg, aes(x, y, col = label)) +
  geom_point(size = 0.001) +
  scale_color_manual(values = pal, name = "Tissue label") +
  guides(col = guide_legend(override.aes = list(size = 2))) +
  theme_void(base_size = 16)

p_feat
# p_feat <- ggplot(gg, aes(x, y, col = label)) +
#   geom_point(size = 0.01) +
#   #scale_color_manual(values = pal) +
#   guides(col = guide_legend(override.aes = list(size = 2))) +
#   theme_void() 

p_feat_nolegend <- p_feat + 
  theme(legend.position = "none") +
  coord_fixed(ratio = 1)  # forza rapporto 1:1

png_file <- "/home/ilaria/Desktop/sid_he/output/thumb/day_11_liver_sid_21.png"
# --- immagine originale ---
img   <- magick::image_read(png_file)
img <- magick::image_flip(img)



p_img <- ggdraw() + draw_image(img)

# --- affianca Original + Features ---
# combo <- plot_grid(p_img, p_feat, ncol = 2,
#                    labels = c("Original", "Features"))

combo <- plot_grid(
  p_img, p_feat_nolegend, ncol = 1,
  align = "hv", axis = "tblr",
  rel_widths = c(1, 1.08)  # puoi giocare con questi valori
)
# Mostra nel report
print(combo)
