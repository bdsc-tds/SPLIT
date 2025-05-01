## ----libs, message=FALSE------------------------------------------------------
if(!requireNamespace("spacexr", quietly = TRUE)){
  remotes::install_github("dmcable/spacexr") ## or remotes::install_github("jpromeror/spacexr@HD") for implementation of the doublet mode.
}
library(spacexr)

if(!requireNamespace("SPLIT", quietly = TRUE)){
  remotes::install_github("bdsc-tds/SPLIT") 
}
library(SPLIT)

library(dplyr)
library(Seurat)
library(readxl)
library(SingleCellExperiment)
library(httr)
library(ggplot2)

## ----load-chormium-metadata---------------------------------------------------
# read metadata
url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-43458-x/MediaObjects/41467_2023_43458_MOESM4_ESM.xlsx"
temp_file <- tempfile(fileext = ".xlsx")
GET(url, write_disk(temp_file, overwrite = TRUE))

chrom_metadata <- read_excel(temp_file, sheet = 1) %>% as.data.frame()
rownames(chrom_metadata) <- chrom_metadata$Barcode

## ----load-chromium------------------------------------------------------------
# read Chromium 
url <- "https://cf.10xgenomics.com/samples/cell-exp/7.0.1/Chromium_FFPE_Human_Breast_Cancer_Chromium_FFPE_Human_Breast_Cancer/Chromium_FFPE_Human_Breast_Cancer_Chromium_FFPE_Human_Breast_Cancer_count_sample_filtered_feature_bc_matrix.h5"
temp_file <- tempfile(fileext = ".h5")
GET(url, write_disk(temp_file, overwrite = TRUE))

chrom_counts <- Read10X_h5(temp_file)
chrom <- CreateSeuratObject(counts = chrom_counts, assay = "RNA", meta.data = chrom_metadata)

chrom$QC <- !is.na(chrom$Annotation)
chrom$is_hybrid <- grepl("Hybrid", chrom$Annotation, ignore.case = TRUE)
chrom <- subset(chrom, subset = QC == TRUE & is_hybrid == FALSE) # remove not annotated cells that cells that have sign of doublets


## ----class-df-----------------------------------------------------------------
cell_type_to_class <- c(
  "B_Cells" = "B cell",
  "CD4+_T_Cells" = "T cell",
  "CD8+_T_Cells" = "T cell",
  "IRF7+_DCs" = "Myeloid",
  "LAMP3+_DCs" = "Myeloid",
  "Macrophages_1" = "Myeloid",
  "Macrophages_2" = "Myeloid",
  "Mast_Cells" = "Myeloid",
  "DCIS 1" = "Epithelial",
  "DCIS 2" = "Epithelial",
  "Invasive_Tumor" = "Epithelial",
  "Prolif_Invasive_Tumor" = "Epithelial",
  "Myoepi_ACTA2+" = "Myoepithelial",
  "Myoepi_KRT15+" = "Myoepithelial",
  "Stromal" = "Stromal",
  "Perivascular-Like" = "Stromal",
  "Endothelial" = "Endothelial"
)

class_df <- data.frame(class = cell_type_to_class)

# and define colors for reproducibility 
library(RColorBrewer)

cell_types <- unique(chrom$Annotation)
colors <- brewer.pal(n = max(3, min(length(cell_types), 12)), name = "Set3")
# Recycle colors if not enough
colors <- rep(colors, length.out = length(cell_types))
pal <- setNames(colors, cell_types)

## ----load-xenium--------------------------------------------------------------
if(!requireNamespace("STexampleData", quietly = TRUE))
  remotes::install_github("lmweber/STexampleData")
xe_full_seu <- STexampleData::Janesick_breastCancer_Xenium_rep1() 

## Convert to Seurat to stay consistent with chromium object
sp_coords <- spatialCoords(xe_full_seu)
colnames(sp_coords) <- c("ST_1", "ST_2")

xe_full <- CreateSeuratObject(
  counts = counts(xe_full_seu),
  assay = "Xenium",
  meta.data = as.data.frame(colData(xe_full_seu))
)

xe_full[["spatial"]] <- CreateDimReducObject(sp_coords, assay = "Xenium", key = "ST_")

xe_full$x <- sp_coords[,1]
xe_full$y <- sp_coords[,2]
rm(xe_full_seu)

## ----downsampling-------------------------------------------------------------
DO_subset_xe <- FALSE 
X_lim <- c(6000, Inf) # cropping area 
Y_lim <- c(4000, Inf) # cropping area 

if(DO_subset_xe){
  xe <- subset(xe_full, subset = x > min(X_lim) & x < max(X_lim) & y > min(Y_lim) & y < max(Y_lim))
} else {
  xe <- xe_full
}

## ----rctd---------------------------------------------------------------------
DO_run_RCTC <- FALSE # FALSE to load pre-computed results

common_genes <- intersect(rownames(xe), rownames(chrom))
ref_labels <- chrom$Annotation %>% as.factor()

ref.obj <- Reference(GetAssayData(chrom, "RNA", "counts")[common_genes, ],
                     cell_types = ref_labels, min_UMI = 10, require_int = TRUE)

test.obj <- SpatialRNA(coords = xe@reductions$spatial@cell.embeddings %>% as.data.frame(),
                       counts = GetAssayData(xe, assay = "Xenium", layer = "counts")[common_genes, ],
                       require_int = TRUE)

if(!exists("class_df")) 
  class_df <- NULL

RCTD <- create.RCTD(
  test.obj, 
  ref.obj, 
  UMI_min = 10, 
  counts_MIN = 10, 
  UMI_min_sigma = 100,
  max_cores = BiocParallel::multicoreWorkers() - 1, 
  CELL_MIN_INSTANCE = 25,
  class_df = class_df # highly recommended 
  )

if(DO_run_RCTC){
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  saveRDS(RCTD, "~/precomp_rctd_class_aware.rds")
} else {
  message("reading precomp RCTD results")
  
  # Install googledrive if you haven't already
  if (!requireNamespace("googledrive", quietly = TRUE)) {
    install.packages("googledrive")
  }
  library(googledrive)
  drive_deauth()
  # Define the file ID from the Google Drive link
  file_id <- "1pTUKq49JbUFwVk7vttjZIFqkx-AKznRF" #"1DCalFIZJywOvrSGBSPqrHh9QINeQp8aq"
  local_path <- tempfile(fileext = ".rds")
  drive_download(as_id(file_id), path = local_path, overwrite = TRUE)
  RCTD <- readRDS(local_path)
  
}


## ----post-rctd----------------------------------------------------------------
RCTD <- SPLIT::run_post_process_RCTD(RCTD)
xe <- AddMetaData(xe, RCTD@results$results_df)
xe <- subset(xe, subset = nCount_Xenium >= 10)

cat("Proprtion of spot classes")
(xe$spot_class %>% table())/ncol(xe)*100 

## ----plot-raw-xenium, fig.width=16, message=FALSE-----------------------------

xe <- xe %>% SCTransform(assay = "Xenium", verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:50, verbose = FALSE)
p1 <- UMAPPlot(xe, group.by = "first_type", label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "none")
p2 <- UMAPPlot(xe, group.by = "second_type", cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom")
p3 <- UMAPPlot(xe, group.by = "spot_class") + theme_void() + theme(aspect.ratio = 1, legend.position = "right")
 
p1 | p2 | p3


## ----spatial-plot, fig.width=12, warning=FALSE--------------------------------
DimPlot(xe, reduction = "spatial", group.by = "first_type", raster = TRUE, cols = pal) + coord_fixed()

## ----SPLIT, message=FALSE-----------------------------------------------------
# Run SPLIT purification
res_split <- SPLIT::purify(
  counts = GetAssayData(xe, assay = 'Xenium', layer = 'counts'), # or any gene x cells counts matrix
  rctd = RCTD,
  DO_purify_singlets = TRUE # Optional. If TRUE, singlets with an available secondary type are purified the same way as doublets_certain; otherwise, left unchanged.
  )

# Create a purified Seurat object
xe_purified <- CreateSeuratObject(
  counts = res_split$purified_counts,
  meta.data = res_split$cell_meta,
  assay = "Xenium"
)

# Optional: Filter, normalize and visualize
xe_purified <- subset(xe_purified, subset = nCount_Xenium > 5)
xe_purified <- xe_purified %>%
  SCTransform(assay = "Xenium", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

#UMAPPlot(xe_purified, group.by = c("first_type"), label = T, repel = T) + theme_void() + theme(aspect.ratio = 1, legend.position = "none")

## ----plot-raw-split-purified, fig.width=12, message=FALSE---------------------
p1 <- UMAPPlot(xe, group.by = c("first_type"), label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "none") + ggtitle("Raw Xenium data")

p2 <- UMAPPlot(xe_purified, group.by = c("first_type"), label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "right") + ggtitle("SPLIT-purified Xenium data")

p3 <- UMAPPlot(xe_purified, group.by = c("spot_class")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom") + ggtitle("SPLIT-purified Xenium data colored by spot class")
p4 <- UMAPPlot(xe_purified, group.by = c("purification_status")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom") + ggtitle("SPLIT-purified Xenium data colored by purification status")


(p1|p2) 
(p3|p4)


## ----spatial-nw---------------------------------------------------------------
sp_nw <- SPLIT::build_spatial_network(
  xe, 
  reduction = "spatial",
  dims = 1:2, 
  DO_prune = TRUE, 
  rad_pruning = 15, # remove connections further than 15um
  k_knn = 20
  )

sp_nw <- SPLIT::add_spatial_metric(spatial_neighborhood = sp_nw, rctd = RCTD)
sp_neigh_df <- SPLIT::neighborhood_analysis_to_metadata(sp_nw)

xe <- AddMetaData(xe, sp_neigh_df)

rm(sp_nw, sp_neigh_df)

## ----plot-neigh-weight-second-type, message=FALSE-----------------------------
# Plot magnitude of local diffusion on UMAP
FeaturePlot(xe, features = c("neighborhood_weights_second_type")) + theme_void() + theme(aspect.ratio = 1)

# Plot distribution of local diffusion value
hist(xe$neighborhood_weights_second_type)

# Plot distribution of local diffusion value per `spot_class`
xe@meta.data %>% filter(!is.na(spot_class)) %>% 
  ggplot(aes(x = spot_class, y = neighborhood_weights_second_type, color = spot_class)) + geom_boxplot() + labs(title = "Local neighbohood diffusion by spot class") + theme_minimal() 

## ----spatially-aware-split, fig.width=12, message=FALSE-----------------------
xe_purified_balanced_score <- SPLIT::balance_raw_and_purified_data_by_score(
  xe_raw = xe,
  xe_purified = xe_purified,
  default_assay = "Xenium", # should be param, but can wait 
  spot_class_key = "spot_class",
  threshold = 0.05, # lower -> more cells will be purified
  score_name = "neighborhood_weights_second_type"
)

# Optional: Filter, normalize and visualize
xe_purified_balanced_score <- subset(xe_purified_balanced_score, subset = nCount_Xenium > 5)
xe_purified_balanced_score <- xe_purified_balanced_score %>%
  SCTransform(assay = "Xenium", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

p5 <- UMAPPlot(xe_purified_balanced_score, group.by = c("first_type"), label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "none") + ggtitle("Spatially-aware SPLIT-purified Xenium data")
p6 <- UMAPPlot(xe_purified_balanced_score, group.by = c("purification_status")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom")

p5|p6

## ----transcriptomics-nw-------------------------------------------------------
tr_nw <- build_transcriptomics_network(
  xe,
  DO_prune = FALSE,
  k_knn = 100
)
tr_nw <- add_transcriptomics_metric(transcriptomics_neighborhood = tr_nw, rctd = RCTD) 
tr_neigh_df <- neighborhood_analysis_to_metadata(tr_nw)
xe <- AddMetaData(xe, tr_neigh_df)

rm(tr_nw, tr_neigh_df)

## ----run-split-shift----------------------------------------------------------
xe_split_shift <- SPLIT::balance_raw_and_purified_data_by_score(
  xe_raw = xe,
  xe_purified = xe_purified,
  default_assay = "Xenium",
  spot_class_key = "spot_class",
  threshold = 0.05, # to be consistent with spatially-aware SPLIT results
  score_name = "neighborhood_weights_second_type",
  DO_swap_lables = TRUE
)

# Optional: Filter, normalize and visualize
xe_split_shift <- subset(xe_split_shift, subset = nCount_Xenium > 5)
xe_split_shift <- xe_split_shift %>%
  SCTransform(assay = "Xenium", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

## ----vis-split-shift, fig.width=16, message=FALSE-----------------------------
p7 <- UMAPPlot(xe_split_shift, group.by = c("first_type"), label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "none") + ggtitle("SPLIT-shift-purified Xenium data")
p8 <- UMAPPlot(xe_split_shift, group.by = c("purification_status")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom")
p9 <- UMAPPlot(xe_split_shift, group.by = c("swap")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom")

p7|p8|p9

## ----vis-split-shift-swap-facet, fig.width=12, message=FALSE------------------
# Visualize results faceted by swapping status
p10 <- UMAPPlot(xe_split_shift, group.by = c("first_type"), split.by = "swap", raster = F, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom") + ggtitle("SPLIT-shift-purified Xenium data faceted by lable swapping status")
p10

## ----summary-plot, fig.width=18, fig.height=7, message=FALSE------------------
p1 | p2+theme(legend.position = "bottom") | p5 | p7

## ----fig.width=12-------------------------------------------------------------
pie_df <- get_pieplot_df(rctd = RCTD)
cois <- xe@meta.data %>% filter(neighborhood_weights_second_type > .2, first_type %in% c("CD4+_T_Cells", "CD8+_T_Cells"), spot_class != "reject") %>% rownames()

p <- plot_pie_around_cell(pie_df, cell_id = cois[1], radius = 50, cols = pal)
plot(p)

