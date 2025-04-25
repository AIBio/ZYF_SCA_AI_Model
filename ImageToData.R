#
setwd("/home/yhw/bioinfo/project-zyf/Release")
ProcessImage <- function(file, sample.n = 1000){
  library("imager")
  library("dplyr")
  library("ggplot2")
  
  # - 1. load image
  file.type <- gsub(".*\\.", "", file)
  file.name <- basename(file)
  image <- load.image(file)
  
  # convert image to grey
  # gray_image <- grayscale(image)
  
  # - 2. compute image gradient
  edges <- imgradient(image, "xy") %>% enorm() %>% threshold("95%")
  
  # - 3. extract coordicates
  edge_coords <- which(edges > 0, arr.ind = TRUE) %>% as.data.frame()
  edge_coords <- subset(edge_coords, dim2 < (max(edge_coords$dim2) - 50))
  dim(edge_coords)
  edge_coords <- edge_coords[sample(1:nrow(edge_coords), sample.n, replace = T), ] %>% 
    dplyr::arrange(dim1)
  dim(edge_coords)

  # - 4. extract x and y
  x_coords <- edge_coords[, "dim1"]
  y_coords <- -edge_coords[, "dim2"] + max(edge_coords[, "dim2"])
  coords <- data.frame(x = x_coords, y = y_coords) %>% 
    rownames_to_column(var = "rank")
  
  # - 5. plot
  pdf(gsub(paste0(file.type, "$"), "pdf", file), height = 4, width = 5)
  print(coords %>% 
          ggplot(aes(x = x, y = y)) + 
          geom_point() + 
          geom_smooth(method = "gam", se = TRUE) +
          theme_bw() +
          theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
                axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
                axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
                axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
                panel.grid = element_blank()))
  dev.off()
  colnames(coords)[2:3] <- paste0(gsub(file.type, "", file.name), colnames(coords)[2:3]) 
  
  # - 6. return data
  #res <- list(coords)
  #names(res) <- file.name
  return(coords)
}
image.file <- list.files("/home/yhw/bioinfo/project-zyf/Release/clustering/primary_release/peaks/", pattern = ".*png", full.names = T)
image.file <- list.files("/home/yhw/bioinfo/project-zyf/Release/MODEL2/MODEL2/picture", pattern = ".*png", full.names = T)
image.file <- list.files("/home/yhw/bioinfo/project-zyf/Release/mouse release/mouse release/ACH", pattern = ".*png", full.names = T)
image.file <- list.files("/home/yhw/bioinfo/project-zyf/Release/mns/mns/MNs", pattern = ".*png", full.names = T)
pd <- list()
for (i in image.file) {
  pd[[basename(i)]] <- ProcessImage(file = i)
}
pd.merge <- do.call(cbind, pd)
pd.merge <- pd.merge[, grep("y$", colnames(pd.merge))]
pd.merge <- apply(pd.merge, 2, function(x){(x - min(x))/(max(x) - min(x))}) %>% t()
#pd.merge <- apply(pd.merge, 2, function(x){(x - mean(x)) / (sd(x))}) %>% t()
rownames(pd.merge) <- paste0("png-", gsub(".png.*", "", rownames(pd.merge)))

# set color of CSI values
csi_palette <- RColorBrewer::brewer.pal(n = 12, name = "Spectral")
cols.rev <- TRUE
if (isTRUE(cols.rev)) {
  csi_palette <- rev(csi_palette)
}
plt.min <- 0
plt.max <- 1
plot.color <- colorRamp2(seq(plt.min, plt.max, (plt.max - plt.min) / (length(csi_palette) - 1)), csi_palette)
library(ComplexHeatmap)
# plot heatmap
ht <- Heatmap(pd.merge, row_km = 25, 
              row_gap = unit(2, "mm"), column_gap = unit(0, "mm"), border = TRUE,
              col = plot.color, name = "Normalized Value",
              cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(3, "cm"),
              clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
              show_column_names = F, column_names_rot = 45, column_names_side = "bottom",
              cluster_columns = F, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
              clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
              show_row_names = F, row_names_rot = 0, row_names_side = "right",
              width = unit(10, "cm"), height = unit(20, "cm"))
print(ht)
ht.width <- 8
ht.height <- 40
res.out <- file.path(getwd(), "clustering")
file.name <- "all_peaks"
method <- "complete"
cols.pal <- "Spectral"
pdf(file.path(res.out, paste0(file.name, "_", method, "_", cols.pal, ".pdf")),
    width = ht.width, height = ht.height
)
print(ht)
dev.off()
# output mode
image.name <- rownames(ht@matrix)
ht <- draw(ht)
dev.off()
image.order <- row_order(ht)
for (i in 1:length(image.order)) {
  write.table(image.name[image.order[[i]]],
              file.path(res.out, paste0(file.name, "_", method, "_", cols.pal, "_Mode_", i, ".txt")),
              quote = F, row.names = F, col.names = F)
  geneset <- list(image.name[image.order[[paste0("Mode.", i)]]])
  names(geneset) <- paste0("Mode.", i)
}

# output data for AI model
pd.merge <- as.data.frame(pd.merge)
colnames(pd.merge) <- gsub("V", "s", colnames(pd.merge))
rownames(pd.merge) <- gsub("png-", "", rownames(pd.merge))
meta <- read.csv("mns/mns/release organoids.csv")
common.group <- intersect(rownames(pd.merge), meta$ID)
meta <- subset(meta, ID %in% common.group) %>% arrange(ID)
pd.merge <- pd.merge[meta$ID, ]
if (all(rownames(pd.merge) == meta$ID)) {
  write.csv(pd.merge, "mns/mns/bin1000_signal.csv")
  write.csv(meta, "mns/mns/bin1000_meta.csv")
}
table(meta$PEAK2)
colnames(meta)
table(meta$PEAK)
# 
save.image("ImageToData.RData")
load("ImageToData.RData")
