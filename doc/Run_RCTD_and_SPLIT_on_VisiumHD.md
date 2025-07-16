RCTD Annotation and SPLIT Purification for VisiumHD Sample
================

- [Introduction](#introduction)
- [Overview](#overview)
  - [Load Data](#load-data)
  - [Read Visium data (counts and cell
    coordinates)](#read-visium-data-counts-and-cell-coordinates)
- [RCDT annotation](#rcdt-annotation)
- [Purification](#purification)
  - [SPLIT (default)](#split-default)
  - [Spatially-aware SPLIT](#spatially-aware-split)
  - [SPLIT-shift](#split-shift)
- [Summary](#summary)

# Introduction

This vignette demonstrates how to annotate VisiumHD spatial
transcriptomics data using RCTD, followed by **purification with
SPLIT**.

⚠️ **Important Notice**

SPLIT currently requires **doublet-mode** RCTD results generated with
the original [spacexr GitHub
repository](https://github.com/dmcable/spacexr) or the faster [HD
fork](https://github.com/jpromeror/spacexr/tree/HD).  
🚧 **Compatibility with the newly released [Bioconductor
version](https://www.bioconductor.org/packages/release/bioc/html/spacexr.html)
of spacexr is under development.**

# Overview

In this vignette, we:

1.  **Run RCTD Annotation**
    - We begin by running RCTD annotation on a VisiumHD sample, using
      matched Chromium data from a public 10x Genomics dataset as the
      reference.
2.  **Apply Default SPLIT Purification**
    - The default SPLIT purification method is then applied to clean the
      annotated VisiumHD sample, refining the initial annotations.
3.  **Apply Spatially-Aware SPLIT**
    - We apply spatially-aware SPLIT, which purifies cells showing signs
      of contamination based on local spatial diffusion patterns.
4.  **Apply SPLIT-Shift**
    - Finally, SPLIT-shift is applied to swap primary and secondary
      labels based on transcriptomic neighborhood heterogeneity,
      improving the accuracy of cell type assignments.

This pipeline assumes that cell type assignments — originally derived
from the Chromium reference — are refined and reliable for downstream
analysis.

``` r
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
library(sf)
library(stringr)
```

## Load Data

For this vignette, we use a publicly available VisiumHD dataset from the
10x Genomics database, originating from:

> **Oliveira, M.F.d., Romero, J.P., Chung, M. et al.** *High-definition
> spatial transcriptomic profiling of immune cell populations in
> colorectal cancer.* *Nature Genetics* 57, 1512–1523 (2025).
> <https://doi.org/10.1038/s41588-025-02193-3>

This dataset provides high-resolution spatial transcriptomics data
suitable for downstream analysis with RCTD and SPLIT.

------------------------------------------------------------------------

### Load Chromium Dataset (Reference)

We manually load metadata for the Chromium single-cell dataset and
metadata from the same study, which will serve as the **reference** for
RCTD annotation.

``` r
url <- "https://raw.githubusercontent.com/10XGenomics/HumanColonCancer_VisiumHD/main/MetaData/SingleCell_MetaData.csv.gz"
temp_file <- tempfile(fileext = ".csv.gz")

GET(url, write_disk(temp_file, overwrite = TRUE))
```

    ## Response [https://raw.githubusercontent.com/10XGenomics/HumanColonCancer_VisiumHD/main/MetaData/SingleCell_MetaData.csv.gz]
    ##   Date: 2025-07-16 13:12
    ##   Status: 200
    ##   Content-Type: application/octet-stream
    ##   Size: 7.29 MB
    ## <ON DISK>  /tmp/RtmpdbVaPH/file149a857084b18.csv.gz

``` r
chrom_metadata <- read.csv(temp_file) #%>% as.data.frame()
rownames(chrom_metadata) <- chrom_metadata$Barcode
```

Manually load Chromium from 10x

``` r
# read Chromium 
url <- "https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5"

temp_file <- tempfile(fileext = ".h5")
GET(url, write_disk(temp_file, overwrite = TRUE))
```

    ## Response [https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5]
    ##   Date: 2025-07-16 13:12
    ##   Status: 200
    ##   Content-Type: application/x-hdf5
    ##   Size: 426 MB
    ## <ON DISK>  /tmp/RtmpdbVaPH/file149a876ed849d.h5

``` r
chrom_counts <- Read10X_h5(temp_file)
chrom <- CreateSeuratObject(counts = chrom_counts, assay = "RNA", meta.data = chrom_metadata)

chrom <- subset(chrom, subset = QCFilter == "Keep") # remove cell that did not pass QC 

rm(chrom_counts)
```

#### Choosing annotation level.

This reference dataset has 2 levels of annotation (ie., `Level1` being
the broader and Level2 being more specific). In this tutorial we will
use `Level1` for simplicity and acceleration. But if you are using more
fine annotation, we **highly recommended** providing `class_df`
parameter into RCTD for more robust RCTD annotation Providing
higher-level cell type classes improves RCTD accuracy and
**significantly** reduces the number of rejected cells, preserving more
cells (rejects are excluded from the downstream analyses and do not
undergo SPLIT purification).

``` r
chrom$Annotation <- chrom$Level1

# Since we use Level1 annotation, there is no need to provide broader level classes, but if you plan using more fine-grained anootation, rpovide their broadee level to enhance RCTD's robustness. For this dataset it would be:
if(FALSE){
  mat <- as.matrix(table(chrom$Level1, chrom$Level2))
  cell_type_to_class <- apply(mat, 2, function(col) rownames(mat)[which.max(col)])
  
  class_df <- data.frame(class = cell_type_to_class)
}
```

``` r
# and define colors for reproducibility 
library(RColorBrewer)

cell_types <- unique(chrom$Annotation)
colors <- brewer.pal(n = max(3, min(length(cell_types), 12)), name = "Set3")
# Recycle colors if not enough
colors <- rep(colors, length.out = length(cell_types))
pal <- setNames(colors, cell_types)
```

### Load VisiumHD Dataset

We begin by downloading the cell-segmented data from the 10x Genomics
website. **Note:** the file is approximately **9 GB**, so you may prefer
to download it manually and reference it from a permanent location. You
can do this via:

``` bash
curl -O https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_segmented_outputs.tar.gz
tar -xvzf Visium_HD_Human_Colon_Cancer_segmented_outputs.tar.gz
```

Or use a tremporary location.

``` r
# download data form 10x counts and coordinates
url <- "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_segmented_outputs.tar.gz"
temp_file <- tempfile(fileext = ".tar.gz")
output_dir <- tempfile()  # or set to a specific folder path

# Download the file
GET(url, write_disk(temp_file, overwrite = TRUE))
```

    ## Response [https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_segmented_outputs.tar.gz]
    ##   Date: 2025-07-16 13:13
    ##   Status: 200
    ##   Content-Type: application/x-tar
    ##   Size: 9.51 GB
    ## <ON DISK>  /tmp/RtmpdbVaPH/file149a887e0763.tar.gz

``` r
# Create output directory if needed
dir.create(output_dir, showWarnings = FALSE)

# Unzip the .tar.gz file
untar(temp_file, exdir = output_dir)

# Check extracted contents
list.files(output_dir, recursive = TRUE)
```

    ##  [1] "segmented_outputs/analysis/clustering/gene_expression_graphclust/clusters.csv"                    
    ##  [2] "segmented_outputs/analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv"            
    ##  [3] "segmented_outputs/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv"             
    ##  [4] "segmented_outputs/analysis/clustering/gene_expression_kmeans_3_clusters/clusters.csv"             
    ##  [5] "segmented_outputs/analysis/clustering/gene_expression_kmeans_4_clusters/clusters.csv"             
    ##  [6] "segmented_outputs/analysis/clustering/gene_expression_kmeans_5_clusters/clusters.csv"             
    ##  [7] "segmented_outputs/analysis/clustering/gene_expression_kmeans_6_clusters/clusters.csv"             
    ##  [8] "segmented_outputs/analysis/clustering/gene_expression_kmeans_7_clusters/clusters.csv"             
    ##  [9] "segmented_outputs/analysis/clustering/gene_expression_kmeans_8_clusters/clusters.csv"             
    ## [10] "segmented_outputs/analysis/clustering/gene_expression_kmeans_9_clusters/clusters.csv"             
    ## [11] "segmented_outputs/analysis/diffexp/gene_expression_graphclust/differential_expression.csv"        
    ## [12] "segmented_outputs/analysis/diffexp/gene_expression_kmeans_10_clusters/differential_expression.csv"
    ## [13] "segmented_outputs/analysis/diffexp/gene_expression_kmeans_2_clusters/differential_expression.csv" 
    ## [14] "segmented_outputs/analysis/diffexp/gene_expression_kmeans_3_clusters/differential_expression.csv" 
    ## [15] "segmented_outputs/analysis/diffexp/gene_expression_kmeans_4_clusters/differential_expression.csv" 
    ## [16] "segmented_outputs/analysis/diffexp/gene_expression_kmeans_5_clusters/differential_expression.csv" 
    ## [17] "segmented_outputs/analysis/diffexp/gene_expression_kmeans_6_clusters/differential_expression.csv" 
    ## [18] "segmented_outputs/analysis/diffexp/gene_expression_kmeans_7_clusters/differential_expression.csv" 
    ## [19] "segmented_outputs/analysis/diffexp/gene_expression_kmeans_8_clusters/differential_expression.csv" 
    ## [20] "segmented_outputs/analysis/diffexp/gene_expression_kmeans_9_clusters/differential_expression.csv" 
    ## [21] "segmented_outputs/analysis/pca/gene_expression_10_components/components.csv"                      
    ## [22] "segmented_outputs/analysis/pca/gene_expression_10_components/dispersion.csv"                      
    ## [23] "segmented_outputs/analysis/pca/gene_expression_10_components/features_selected.csv"               
    ## [24] "segmented_outputs/analysis/pca/gene_expression_10_components/projection.csv"                      
    ## [25] "segmented_outputs/analysis/pca/gene_expression_10_components/variance.csv"                        
    ## [26] "segmented_outputs/analysis/umap/gene_expression_2_components/projection.csv"                      
    ## [27] "segmented_outputs/cell_segmentations.geojson"                                                     
    ## [28] "segmented_outputs/cloupe.cloupe"                                                                  
    ## [29] "segmented_outputs/filtered_feature_cell_matrix.h5"                                                
    ## [30] "segmented_outputs/filtered_feature_cell_matrix/barcodes.tsv.gz"                                   
    ## [31] "segmented_outputs/filtered_feature_cell_matrix/features.tsv.gz"                                   
    ## [32] "segmented_outputs/filtered_feature_cell_matrix/matrix.mtx.gz"                                     
    ## [33] "segmented_outputs/graphclust_annotated_cell_segmentations.geojson"                                
    ## [34] "segmented_outputs/graphclust_annotated_nucleus_segmentations.geojson"                             
    ## [35] "segmented_outputs/nucleus_segmentations.geojson"                                                  
    ## [36] "segmented_outputs/raw_feature_cell_matrix.h5"                                                     
    ## [37] "segmented_outputs/raw_feature_cell_matrix/barcodes.tsv.gz"                                        
    ## [38] "segmented_outputs/raw_feature_cell_matrix/features.tsv.gz"                                        
    ## [39] "segmented_outputs/raw_feature_cell_matrix/matrix.mtx.gz"                                          
    ## [40] "segmented_outputs/spatial/aligned_fiducials.jpg"                                                  
    ## [41] "segmented_outputs/spatial/aligned_tissue_image.jpg"                                               
    ## [42] "segmented_outputs/spatial/cytassist_image.tiff"                                                   
    ## [43] "segmented_outputs/spatial/detected_tissue_image.jpg"                                              
    ## [44] "segmented_outputs/spatial/scalefactors_json.json"                                                 
    ## [45] "segmented_outputs/spatial/tissue_hires_image.png"                                                 
    ## [46] "segmented_outputs/spatial/tissue_lowres_image.png"

## Read Visium data (counts and cell coordinates)

``` r
vhd_counts <- Seurat::Read10X_h5(file.path(output_dir, "segmented_outputs/filtered_feature_cell_matrix.h5"))

# get VisiumHD cell coordinates (optionally, RCTD works w/o coordinates as well)
cell_segmentation_path <- file.path(output_dir, "segmented_outputs/cell_segmentations.geojson")
cell_segmentation <- st_read(cell_segmentation_path)
```

    ## Reading layer `cell_segmentations' from data source `/tmp/RtmpdbVaPH/file149a874e8d6ae/segmented_outputs/cell_segmentations.geojson' using driver `GeoJSON'
    ## Simple feature collection with 220704 features and 1 field
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 40587.5 ymin: -28.8 xmax: 65202.7 ymax: 22701
    ## Geodetic CRS:  WGS 84

``` r
st_crs(cell_segmentation) <- 32632 
```

    ## Warning: st_crs<- : replacing crs does not reproject data; use st_transform for that

``` r
cell_segmentation <- cell_segmentation %>%
  mutate(cell_name = str_glue("cellid_{str_pad(cell_id, 9, pad = '0')}-1"))
  
vhd_coords <- st_centroid(cell_segmentation)
```

    ## Warning: st_centroid assumes attributes are constant over geometries

``` r
rownames(vhd_coords) <- vhd_coords$cell_name
  
vhd_coords <- vhd_coords %>%
  mutate(ST_1 = st_coordinates(.)[,1],
         ST_2 = st_coordinates(.)[,2]) %>%
  st_set_geometry(NULL) %>% 
  select(ST_1, ST_2)
```

``` r
## Convert to Seurat to stay consistent with chromium object
vhd <- CreateSeuratObject(
  counts = vhd_counts,
  assay = "VHD",
  meta.data = vhd_coords
)

vhd[["spatial"]] <- CreateDimReducObject(vhd_coords[colnames(vhd),] %>% as.matrix(), assay = "VHD", key = "ST_")

vhd$x <- vhd$ST_1
vhd$y <- vhd$ST_2

vhd <- subset(vhd, subset = nCount_VHD > 100)
rm(vhd_counts, vhd_coords)
```

``` r
DimPlot(vhd, reduction = "spatial") # simple visualization
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggplot(cell_segmentation) +
  geom_sf(alpha = .5,) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Cell Segmentation") +
  scale_fill_manual(values = pal)
```

    ## Warning: No shared levels found between `names(values)` of the manual scale and the data's fill values.

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

#### Downsample VisiumHD Dataset for Faster RCTD (Optional)

Running RCTD can be time-consuming, especially on large datasets like
VisiumHD To speed up computation during this tutorial, we optionally
provide code to downsample the dataset. We recommend spatial cropping
(rather than random sampling) to preserve neighborhood structure, which
is important for downstream analysis.

That said, downsampling is **not required **, we provide precomputed
RCTD results below so you can skip running RCTD altogether if desired.

``` r
DO_subset_vhd <- TRUE 
X_lim <- c(45000, 60000) # cropping area 
Y_lim <- c(10000, 20000) # cropping area 

if(DO_subset_vhd){
  vhd <- subset(vhd, subset = x > min(X_lim) & x < max(X_lim) & y > min(Y_lim) & y < max(Y_lim))
} 
```

# RCDT annotation

Run RCTD Annotation on VisiumHD Running RCTD on large datasets can be
computationally intensive and may take several hours. To streamline the
workflow, we provide the full code for reproducibility. However, we
recommend loading a pre-computed RCTD object by setting
`DO_run_RCTD <- FALSE`.

``` r
DO_run_RCTC <- FALSE # FALSE to load pre-computed results

common_genes <- intersect(rownames(vhd), rownames(chrom))
ref_labels <- chrom$Annotation %>% as.factor()

ref.obj <- Reference(GetAssayData(chrom, "RNA", "counts")[common_genes, ],
                     cell_types = ref_labels, require_int = TRUE)
```

    ## Warning in Reference(GetAssayData(chrom, "RNA", "counts")[common_genes, : Reference: number of cells per cell type is 63900, larger than maximum
    ## allowable of 10000. Downsampling number of cells to: 10000

``` r
test.obj <- SpatialRNA(coords = vhd@reductions$spatial@cell.embeddings %>% as.data.frame(),
                       counts = GetAssayData(vhd, assay = "VHD", layer = "counts")[common_genes, ],
                       require_int = TRUE)

if(!exists("class_df")) 
  class_df <- NULL

rctd <- create.RCTD(
  test.obj, 
  ref.obj, 
  max_cores = 10, 
  class_df = class_df # highly recommended if annotation provided at the fine level
  )
```

    ## Begin: process_cell_type_info

    ## process_cell_type_info: number of cells in reference: 82237

    ## process_cell_type_info: number of genes in reference: 17696

    ## 
    ##               B cells           Endothelial            Fibroblast Intestinal Epithelial               Myeloid              Neuronal 
    ##                 10000                  7916                 10000                 10000                 10000                  4321 
    ##         Smooth Muscle               T cells                 Tumor 
    ##                 10000                 10000                 10000

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size 1.3 GiB

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size 1.0 GiB

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size 1.3 GiB
    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size 1.3 GiB
    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size 1.3 GiB
    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size 1.3 GiB
    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size 1.3 GiB
    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size 1.3 GiB

    ## End: process_cell_type_info

    ## create.RCTD: getting regression differentially expressed genes:

    ## get_de_genes: B cells found DE genes: 135

    ## get_de_genes: Endothelial found DE genes: 367

    ## get_de_genes: Fibroblast found DE genes: 257

    ## get_de_genes: Intestinal Epithelial found DE genes: 335

    ## get_de_genes: Myeloid found DE genes: 376

    ## get_de_genes: Neuronal found DE genes: 294

    ## get_de_genes: Smooth Muscle found DE genes: 347

    ## get_de_genes: T cells found DE genes: 460

    ## get_de_genes: Tumor found DE genes: 285

    ## get_de_genes: total DE genes: 2379

    ## create.RCTD: getting platform effect normalization differentially expressed genes:

    ## get_de_genes: B cells found DE genes: 271

    ## get_de_genes: Endothelial found DE genes: 708

    ## get_de_genes: Fibroblast found DE genes: 481

    ## get_de_genes: Intestinal Epithelial found DE genes: 714

    ## get_de_genes: Myeloid found DE genes: 729

    ## get_de_genes: Neuronal found DE genes: 655

    ## get_de_genes: Smooth Muscle found DE genes: 740

    ## get_de_genes: T cells found DE genes: 1042

    ## get_de_genes: Tumor found DE genes: 792

    ## get_de_genes: total DE genes: 4514

``` r
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
  file_id <- "1wgwYQVdYzJbsAt0bYW-Is9eHnCSTzNMh"
  local_path <- tempfile(fileext = ".rds")
  drive_download(as_id(file_id), path = local_path, overwrite = TRUE)
  RCTD <- readRDS(local_path)
}
```

    ## reading precomp RCTD results

    ## File downloaded:

    ## • 'precomp_rctd_visHD_lite.rds' <id: 1wgwYQVdYzJbsAt0bYW-Is9eHnCSTzNMh>

    ## Saved locally as:

    ## • '/tmp/RtmpdbVaPH/file149a84dbf91d5.rds'

``` r
rm(chrom)
gc()
```

    ##             used   (Mb) gc trigger    (Mb)   max used    (Mb)
    ## Ncells  15807651  844.3   26650366  1423.3   26650366  1423.3
    ## Vcells 374623409 2858.2 1636850072 12488.2 2534206190 19334.5

Visualize RCTD Annotation Post-process RCDT output and add results into
VisiumHD object

``` r
RCTD <- SPLIT::run_post_process_RCTD(RCTD)
```

    ## Correcting singlets ...

    ## Updating scores ...

    ## Add coordinates to results ...

    ## Computing alternative annotations ...

    ## Replacing results_df ...

``` r
vhd <- AddMetaData(vhd, RCTD@results$results_df)
vhd <- subset(vhd, subset = nCount_VHD >= 100)
```

``` r
cat("Proprtion of spot classes")
```

    ## Proprtion of spot classes

``` r
(vhd$spot_class %>% table())/ncol(vhd)*100 
```

    ## .
    ##            reject doublet_uncertain   doublet_certain           singlet 
    ##          1.439481          4.885608         23.076431         70.598479

``` r
vhd <- vhd %>% NormalizeData() 
vhd <- vhd %>% FindVariableFeatures()
vhd <- vhd %>% ScaleData()
vhd <- vhd %>% RunPCA()
vhd <- vhd %>% RunUMAP(dims = 1:50)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

``` r
p1 <- UMAPPlot(vhd, group.by = "first_type", label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "none")
p2 <- UMAPPlot(vhd, group.by = "second_type", cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom")
p3 <- UMAPPlot(vhd, group.by = "spot_class") + theme_void() + theme(aspect.ratio = 1, legend.position = "right")
 
p1 | p2 | p3
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/plot-raw-visiumHD-1.png)<!-- -->
Spatial Visualization

``` r
DimPlot(vhd, reduction = "spatial", group.by = "first_type", raster = TRUE, cols = pal) + coord_fixed()
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/spatial-plot-1.png)<!-- -->

``` r
cell_segmentation <- cell_segmentation %>%
  filter(cell_name %in% rownames(vhd@meta.data))

cell_segmentation <- cell_segmentation %>%
  left_join(
    vhd@meta.data %>%
      select(first_type, second_type, spot_class) %>%
      tibble::rownames_to_column("cell_name"),
    by = "cell_name"
  )

ggplot(cell_segmentation) +
  geom_sf(aes(fill = first_type, color = first_type), alpha = .5) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Cell Segmentation") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal)
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/spatial-plot-2.png)<!-- -->

# Purification

## SPLIT (default)

This section runs the default SPLIT purification and visualizes purified
data.

``` r
# Run SPLIT purification
res_split <- purify(
  counts = GetAssayData(vhd, assay = 'VHD', layer = 'counts'), # or any gene x cells counts matrix
  rctd = RCTD,
  DO_parallel = F,
  n_workers = NULL,
  chunk_size = 5000,
  DO_purify_singlets = T # Optional. If TRUE, singlets with an available secondary type are purified the same way as doublets_certain; otherwise, left unchanged.
)
```

    ## [1] 17696 62592
    ## N_genes =  17696Processing certain doublets...
    ## 58633 
    ## 0 %
    ## 9 %
    ## 17 %
    ## 26 %
    ## 34 %
    ## 43 %
    ## 51 %
    ## 60 %
    ## 68 %
    ## 77 %
    ## 85 %
    ## 94 %
    ## Processing uncertain doublets...
    ## 3058 
    ## 0 %
    ## Purification completed in  5.559927object.size(all_doublet_results): 552409040Combaning doublets results ...

    ## Processed 10000 / 61691
    ## Processed 20000 / 61691
    ## Processed 30000 / 61691
    ## Processed 40000 / 61691
    ## Processed 50000 / 61691
    ## Processed 60000 / 61691
    ## Processed 61691 / 61691
    ## Combining results completed in  3.320345object.size(purified): 466812768

``` r
# Create a purified Seurat object
vhd_purified <- CreateSeuratObject(
  counts = res_split$purified_counts,
  meta.data = res_split$cell_meta,
  assay = "VHD"
)

# Optional: Filter, normalize and visualize
vhd_purified <- subset(vhd_purified, subset = nCount_VHD > 100)
vhd_purified <- vhd_purified %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:50)
```

### Visually compare results of Raw and SPLIT-Purified data

``` r
p1 <- UMAPPlot(vhd, group.by = c("first_type"), label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "none") + ggtitle("Raw VisiumHD data")

p2 <- UMAPPlot(vhd_purified, group.by = c("first_type"), label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "right") + ggtitle("SPLIT-purified VisiumHD data")

p3 <- UMAPPlot(vhd_purified, group.by = c("spot_class")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom") + ggtitle("SPLIT-purified VisiumHD data colored by spot class")
p4 <- UMAPPlot(vhd_purified, group.by = c("purification_status")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom") + ggtitle("SPLIT-purified VisiumHD data colored by purification status")


(p1|p2) 
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/plot-raw-split-purified-1.png)<!-- -->

``` r
(p3|p4)
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/plot-raw-split-purified-2.png)<!-- -->

## Spatially-aware SPLIT

SPLIT can leverage spatial information to assess the abundance of
secondary signals in the local neighborhood (i.e., local diffusion
potential), enabling selective decomposition only when contamination is
likely. This spatially informed strategy helps prevent overcorrection of
phenotypes that may be underrepresented or absent in the reference.
Specifically, we first compute the spatial neighborhood for each cell,
then identify and purify cells that have sign of local diffusion of the
secondary cell type.

``` r
sp_nw <- SPLIT::build_spatial_network(
  vhd, 
  reduction = "spatial",
  dims = 1:2, 
  DO_prune = TRUE, 
  rad_pruning = 50, # remove connections further than 15um
  k_knn = 20
  )
```

    ## Computing nearest neighbors

    ## Only one graph name supplied, storing nearest-neighbor graph only

    ## N = 883146 ( 74 %) edges were pruned

``` r
sp_nw <- SPLIT::add_spatial_metric(spatial_neighborhood = sp_nw, rctd = RCTD)
sp_neigh_df <- SPLIT::neighborhood_analysis_to_metadata(sp_nw)

vhd <- AddMetaData(vhd, sp_neigh_df)

rm(sp_nw, sp_neigh_df)
```

### Visualize local diffusion of secondary cell type

The score `neighborhood_weights_second_type` corresponds to the average
weight of the secondary cell type in cell’s spatial neighborhood.

``` r
# Plot magnitude of local diffusion on UMAP
FeaturePlot(vhd, features = c("neighborhood_weights_second_type")) + theme_void() + theme(aspect.ratio = 1)
```

    ## Warning: The `slot` argument of `FetchData()` is deprecated as of SeuratObject 5.0.0.
    ## ℹ Please use the `layer` argument instead.
    ## ℹ The deprecated feature was likely used in the Seurat package.
    ##   Please report the issue at <https://github.com/satijalab/seurat/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/plot-neigh-weight-second-type-1.png)<!-- -->

``` r
# Plot distribution of local diffusion value
hist(vhd$neighborhood_weights_second_type)
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/plot-neigh-weight-second-type-2.png)<!-- -->

``` r
# Plot distribution of local diffusion value per `spot_class`
vhd@meta.data %>% filter(!is.na(spot_class)) %>% 
  ggplot(aes(x = spot_class, y = neighborhood_weights_second_type, color = spot_class)) + geom_boxplot() + labs(title = "Local neighbohood diffusion by spot class") + theme_minimal() 
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/plot-neigh-weight-second-type-3.png)<!-- -->

We now purify cells that have secondary signal in their spatial
neighborhood (e.g., `neighborhood_weights_second_type`) and keep other
cells unchanged

``` r
vhd_purified_balanced_score <- SPLIT::balance_raw_and_purified_data_by_score(
  xe_raw = vhd,
  xe_purified = vhd_purified,
  default_assay = "VHD", # 
  spot_class_key = "spot_class",
  threshold = 0.05, # lower -> more cells will be purified
  score_name = "neighborhood_weights_second_type"
)

# Optional: Filter, normalize and visualize
vhd_purified_balanced_score <- subset(vhd_purified_balanced_score, subset = nCount_VHD > 100)
vhd_purified_balanced_score <- vhd_purified_balanced_score %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:50, verbose = FALSE)

p5 <- UMAPPlot(vhd_purified_balanced_score, group.by = c("first_type"), label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "none") + ggtitle("Spatially-aware SPLIT-purified VisiumHD data")
p6 <- UMAPPlot(vhd_purified_balanced_score, group.by = c("purification_status")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom")

p5|p6
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/spatially-aware-split-1.png)<!-- -->

## SPLIT-shift

In some cases, the contamination signal is so strong that RCTD assigns
the cell to its secondary cell type. To address this, we introduce
SPLIT-shift—an approach that refines phenotype assignments by swapping
the primary and secondary cell type labels based on transcriptional
neighborhood homogeneity.

For this, we need to compute transcriptomics neighborhood

``` r
tr_nw <- build_transcriptomics_network(
  vhd,
  DO_prune = FALSE,
  k_knn = 100
)
```

    ## Computing nearest neighbors

    ## Only one graph name supplied, storing nearest-neighbor graph only

``` r
tr_nw <- add_transcriptomics_metric(transcriptomics_neighborhood = tr_nw, rctd = RCTD) 
tr_neigh_df <- neighborhood_analysis_to_metadata(tr_nw)
vhd <- AddMetaData(vhd, tr_neigh_df)

rm(tr_nw, tr_neigh_df)
```

And then, we set `DO_swap_lables = TRUE` to allow SPLIT-shift

``` r
vhd_split_shift <- SPLIT::balance_raw_and_purified_data_by_score(
  xe_raw = vhd,
  xe_purified = vhd_purified,
  default_assay = "VHD",
  spot_class_key = "spot_class",
  threshold = 0.05, # to be consistent with spatially-aware SPLIT results
  score_name = "neighborhood_weights_second_type",
  DO_swap_lables = TRUE
)

# Optional: Filter, normalize and visualize
vhd_split_shift <- subset(vhd_split_shift, subset = nCount_VHD > 100)
vhd_split_shift <- vhd_split_shift %>%
  SCTransform(assay = "VHD", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:50, verbose = FALSE)
```

    ## Warning: The `slot` argument of `GetAssayData()` is deprecated as of SeuratObject 5.0.0.
    ## ℹ Please use the `layer` argument instead.
    ## ℹ The deprecated feature was likely used in the Seurat package.
    ##   Please report the issue at <https://github.com/satijalab/seurat/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

    ## Warning: The `slot` argument of `SetAssayData()` is deprecated as of SeuratObject 5.0.0.
    ## ℹ Please use the `layer` argument instead.
    ## ℹ The deprecated feature was likely used in the Seurat package.
    ##   Please report the issue at <https://github.com/satijalab/seurat/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

### Visualize SPLIT-shift

``` r
p7 <- UMAPPlot(vhd_split_shift, group.by = c("first_type"), label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "none") + ggtitle("SPLIT-shift-purified VisiumHD data")
p8 <- UMAPPlot(vhd_split_shift, group.by = c("purification_status")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom")
p9 <- UMAPPlot(vhd_split_shift, group.by = c("swap")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom")

p7|p8|p9
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/vis-split-shift-1.png)<!-- -->

``` r
# Visualize results faceted by swapping status
p10 <- UMAPPlot(vhd_split_shift, group.by = c("first_type"), split.by = "swap", raster = F, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom") + ggtitle("SPLIT-shift-purified VisiumHD data faceted by lable swapping status")
p10
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/vis-split-shift-swap-facet-1.png)<!-- -->

# Summary

To run SPLIT, you need a single-cell reference with reliable cell type
labels, which is used to annotate VisiumHD data using RCTD in doublet
mode.  
We strongly recommend providing a broader-level mapping of the reference
cell types to higher-level classes. This helps RCTD produce more robust
results and reduces the number of rejected cells, which are excluded
from downstream analysis.

After annotation, SPLIT can be run in several **combinable** modes to
purify the data:

1.  **Default SPLIT** purifies all `doublets_certain`, all
    `doublets_uncertain`, and—if `DO_purify_singlets = TRUE` — singlets
    that show signs of a secondary cell type.  
    Rejected cells are always removed as unreliable.

2.  **Spatially-aware SPLIT** purifies any cells that show signs of
    contamination based on local spatial diffusion — i.e., having
    secondary signal in their spatial neighborhood.

3.  **SPLIT-shift** allows swapping the primary and secondary cell type
    labels based on transcriptional neighborhood homogeneity.

``` r
p1 | p2+theme(legend.position = "bottom") | p5 | p7
```

![](Run_RCTD_and_SPLIT_on_VisiumHD_files/figure-gfm/summary-plot-1.png)<!-- -->
