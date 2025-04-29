## ----message=FALSE------------------------------------------------------------
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

## -----------------------------------------------------------------------------
# read metadata
url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-43458-x/MediaObjects/41467_2023_43458_MOESM4_ESM.xlsx"
temp_file <- tempfile(fileext = ".xlsx")
GET(url, write_disk(temp_file, overwrite = TRUE))

chrom_metadata <- read_excel(temp_file, sheet = 1) %>% as.data.frame()
rownames(chrom_metadata) <- chrom_metadata$Barcode

## -----------------------------------------------------------------------------
# read Chromium 
url <- "https://cf.10xgenomics.com/samples/cell-exp/7.0.1/Chromium_FFPE_Human_Breast_Cancer_Chromium_FFPE_Human_Breast_Cancer/Chromium_FFPE_Human_Breast_Cancer_Chromium_FFPE_Human_Breast_Cancer_count_sample_filtered_feature_bc_matrix.h5"
temp_file <- tempfile(fileext = ".h5")
GET(url, write_disk(temp_file, overwrite = TRUE))

chrom_counts <- Read10X_h5(temp_file)
chrom <- CreateSeuratObject(counts = chrom_counts, assay = "RNA", meta.data = chrom_metadata)

chrom$QC <- !is.na(chrom$Annotation)
chrom <- subset(chrom, subset = QC == TRUE)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
DO_subset_xe <- FALSE 
X_lim <- c(6000, Inf)
Y_lim <- c(4000, Inf)

if(DO_subset_xe){
  xe <- subset(xe_full, subset = x > min(X_lim) & x < max(X_lim) & y > min(Y_lim) & y < max(Y_lim))
} else {
  xe <- xe_full
}

## -----------------------------------------------------------------------------
DO_run_RCTC <- FALSE 

common_genes <- intersect(rownames(xe), rownames(chrom))
ref_labels <- chrom$Annotation %>% as.factor()

ref.obj <- Reference(GetAssayData(chrom, "RNA", "counts")[common_genes, ],
                     cell_types = ref_labels, min_UMI = 10, require_int = TRUE)

test.obj <- SpatialRNA(coords = xe@reductions$spatial@cell.embeddings %>% as.data.frame(),
                       counts = GetAssayData(xe, assay = "Xenium", layer = "counts")[common_genes, ],
                       require_int = TRUE)

RCTD <- create.RCTD(test.obj, ref.obj, UMI_min = 10, counts_MIN = 10, UMI_min_sigma = 100, max_cores = 20, CELL_MIN_INSTANCE = 25)

if(DO_run_RCTC){
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  saveRDS(RCTD, "~/precomp_rctd.rds")
} else {
  message("reading precomp RCTD results")
  
  # Install googledrive if you haven't already
  if (!requireNamespace("googledrive", quietly = TRUE)) {
    install.packages("googledrive")
  }
  library(googledrive)
  drive_deauth()
  # Define the file ID from the Google Drive link
  file_id <- "1DCalFIZJywOvrSGBSPqrHh9QINeQp8aq"
  local_path <- tempfile(fileext = ".rds")
  drive_download(as_id(file_id), path = local_path, overwrite = TRUE)
  RCTD <- readRDS(local_path)
  
}


## -----------------------------------------------------------------------------
RCTD <- SPLIT::run_post_process_RCTD(RCTD)
xe <- AddMetaData(xe, RCTD@results$results_df)
xe <- subset(xe, subset = nCount_Xenium >= 10)

## ----fig.width=16-------------------------------------------------------------

xe <- xe %>% SCTransform(assay = "Xenium") %>% RunPCA() %>% RunUMAP(dims = 1:50)
p1 <- UMAPPlot(xe, group.by = "first_type", label = T, repel = T) + theme_void() + theme(aspect.ratio = 1, legend.position = "none")
p2 <- UMAPPlot(xe, group.by = "second_type") + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom")
p3 <- UMAPPlot(xe, group.by = "spot_class") + theme_void() + theme(aspect.ratio = 1, legend.position = "right")
 
p1 | p2 | p3


## -----------------------------------------------------------------------------
# Run SPLIT purification
res_split <- SPLIT::purify(
  counts = GetAssayData(xe, assay = 'Xenium', layer = 'counts'), # or any gene x cells counts matrix
  rctd = RCTD,
  DO_purify_singlets = TRUE # optional
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
  SCTransform(assay = "Xenium") %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20)

UMAPPlot(xe_purified, group.by = c("first_type"), label = T, repel = T) + theme_void() + theme(aspect.ratio = 1, legend.position = "none")

## ----fig.height=8-------------------------------------------------------------
p1 <- UMAPPlot(xe, group.by = c("first_type"), label = T, repel = T) + theme_void() + theme(aspect.ratio = 1, legend.position = "none") + ggtitle("Raw Xenium data")

p2 <- UMAPPlot(xe_purified, group.by = c("first_type"), label = T, repel = T) + theme_void() + theme(aspect.ratio = 1, legend.position = "right") + ggtitle("SPLIT-purified Xenium data")

p1|p2


## ----fig.height=8-------------------------------------------------------------
DimPlot(xe, reduction = "spatial", group.by = "first_type") + coord_fixed()

