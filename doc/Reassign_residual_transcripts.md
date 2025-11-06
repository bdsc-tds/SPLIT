Reassigning residual trancripts to neighboring cells
================

- [Introduction](#introduction)
- [Purification](#purification)
  - [SPLIT (default)](#split-default)
- [Reassignment of filtered out
  transcripts](#reassignment-of-filtered-out-transcripts)
  - [Compute spatial network](#compute-spatial-network)
  - [Evaluate transcript
    re-assigment](#evaluate-transcript-re-assigment)
- [Exploting cells that cannot send their residual
  transcripts](#exploting-cells-that-cannot-send-their-residual-transcripts)
- [Do residuals even make sense?](#do-residuals-even-make-sense)

# Introduction

This vignette demonstrates how to reassign residual (filtered out)
transcripts to the neighboring cells. For this, we will use the results
of the [vignette on running SPLIT on
Xenium](https://github.com/bdsc-tds/SPLIT/blob/main/vignettes/Run_RCTD_and_SPLIT_on_Xenium.Rmd)

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
library(SingleCellExperiment)
library(ggplot2)
```

### Load Xenium Dataset

We load the Xenium spatial transcriptomics data using the
`STexampleData` package, which provides convenient access to example
spatial datasets for analysis and demonstration.

``` r
if(!requireNamespace("STexampleData", quietly = TRUE))
  remotes::install_github("lmweber/STexampleData")
xe_full_seu <- STexampleData::Janesick_breastCancer_Xenium_rep1() 
```

    ## see ?STexampleData and browseVignettes('STexampleData') for documentation

    ## loading from cache

``` r
## Convert to Seurat to stay consistent with chromium object
sp_coords <- spatialCoords(xe_full_seu)
colnames(sp_coords) <- c("ST_1", "ST_2")

xe <- CreateSeuratObject(
  counts = counts(xe_full_seu),
  assay = "Xenium",
  meta.data = as.data.frame(colData(xe_full_seu))
)

xe[["spatial"]] <- CreateDimReducObject(sp_coords, assay = "Xenium", key = "ST_")

xe$x <- sp_coords[,1]
xe$y <- sp_coords[,2]
rm(xe_full_seu)
```

You can optionally crop the slice to reduce computations, but we will
use the entire dataset.

``` r
DO_subset_xe <- FALSE 
X_lim <- c(6000, Inf) # cropping area 
Y_lim <- c(4000, Inf) # cropping area 

if(DO_subset_xe){
  xe <- subset(xe, subset = x > min(X_lim) & x < max(X_lim) & y > min(Y_lim) & y < max(Y_lim))
} 
```

Download pre-computed RCTD (ie., cell-type annotation) results

``` r
# Install googledrive if you haven't already
if (!requireNamespace("googledrive", quietly = TRUE)) {
  install.packages("googledrive")
}
library(googledrive)
drive_deauth()
# Define the file ID from the Google Drive link
file_id <- "1pTUKq49JbUFwVk7vttjZIFqkx-AKznRF" 
local_path <- tempfile(fileext = ".rds")
drive_download(as_id(file_id), path = local_path, overwrite = TRUE)
```

    ## File downloaded:

    ## • 'precomp_rctd_class_aware.rds' <id: 1pTUKq49JbUFwVk7vttjZIFqkx-AKznRF>

    ## Saved locally as:

    ## • '/tmp/Rtmpt49rM9/file2114f52c2b811.rds'

``` r
RCTD <- readRDS(local_path)
```

Visualize RCTD Annotation Post-process RCDT output and add results into
Xenium object

``` r
RCTD <- SPLIT::run_post_process_RCTD(RCTD)
```

    ## Correcting singlets ...

    ## Updating scores ...

    ## Add coordinates to results ...

    ## Computing alternative annotations ...

    ## Replacing results_df ...

``` r
xe <- AddMetaData(xe, RCTD@results$results_df)

xe <- subset(xe, subset = nCount_Xenium >= 10)

xe <- xe %>% SCTransform(assay = "Xenium", verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:50, verbose = FALSE)
```

``` r
# and define colors for reproducibility 
library(RColorBrewer)

cell_types <- unique(xe$first_type) %>% sort()
colors <- brewer.pal(n = max(3, min(length(cell_types), 12)), name = "Set3")
# Recycle colors if not enough
colors <- rep(colors, length.out = length(cell_types))
pal <- setNames(colors, cell_types)
```

``` r
p1 <- UMAPPlot(xe, group.by = "first_type", label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "none")
p2 <- UMAPPlot(xe, group.by = "second_type", cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom")
p3 <- UMAPPlot(xe, group.by = "spot_class") + theme_void() + theme(aspect.ratio = 1, legend.position = "right")
 
p1 | p2 | p3
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/plot-raw-xenium-1.png)<!-- -->

Spatial Visualization

``` r
DimPlot(xe, reduction = "spatial", group.by = "first_type", raster = F, cols = pal) + coord_fixed()
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/spatial-plot-1.png)<!-- -->

# Purification

## SPLIT (default)

This section runs the default SPLIT purification and visualizes purified
data.

``` r
# Run SPLIT purification
res_split <- SPLIT::purify(
  counts = GetAssayData(xe, assay = 'Xenium', layer = 'counts'), # or any gene x cells counts matrix
  rctd = RCTD,
  DO_purify_singlets = TRUE 
)
```

    ## [1]    307 163849
    ## N_genes = 307 
    ## Processing certain doublets...
    ## 158702 
    ## 0 %
    ## 6 %
    ## 13 %
    ## 19 %
    ## 25 %
    ## 32 %
    ## 38 %
    ## 44 %
    ## 50 %
    ## 57 %
    ## 63 %
    ## 69 %
    ## 76 %
    ## 82 %
    ## 88 %
    ## 95 %
    ## Processing uncertain doublets...
    ## 2830 
    ## 0 %
    ## Combaning doublets results ...

    ## Processed 10000 / 161532
    ## Processed 20000 / 161532
    ## Processed 30000 / 161532
    ## Processed 40000 / 161532
    ## Processed 50000 / 161532
    ## Processed 60000 / 161532
    ## Processed 70000 / 161532
    ## Processed 80000 / 161532
    ## Processed 90000 / 161532
    ## Processed 100000 / 161532
    ## Processed 110000 / 161532
    ## Processed 120000 / 161532
    ## Processed 130000 / 161532
    ## Processed 140000 / 161532
    ## Processed 150000 / 161532
    ## Processed 160000 / 161532
    ## Processed 161532 / 161532

``` r
# Create a purified Seurat object
xe_purified <- CreateSeuratObject(
  counts = res_split$purified_counts,
  meta.data = res_split$cell_meta,
  assay = "Xenium"
)
```

Running a standard analysis on purified data to get new (purified) UMAP

``` r
# Optional: Filter, normalize and visualize
xe_purified <- subset(xe_purified, subset = nCount_Xenium > 5)
xe_purified <- xe_purified %>%
  SCTransform(assay = "Xenium", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

#UMAPPlot(xe_purified, group.by = c("first_type"), label = T, repel = T) + theme_void() + theme(aspect.ratio = 1, legend.position = "none")
```

### Visually compare results of Raw and SPLIT-Purified data

``` r
p1 <- UMAPPlot(xe, group.by = c("first_type"), label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "none") + ggtitle("Raw Xenium data")

p2 <- UMAPPlot(xe_purified, group.by = c("first_type"), label = T, repel = T, cols = pal) + theme_void() + theme(aspect.ratio = 1, legend.position = "right") + ggtitle("SPLIT-purified Xenium data")

p3 <- UMAPPlot(xe_purified, group.by = c("spot_class")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom") + ggtitle("SPLIT-purified Xenium data colored by spot class")
p4 <- UMAPPlot(xe_purified, group.by = c("purification_status")) + theme_void() + theme(aspect.ratio = 1, legend.position = "bottom") + ggtitle("SPLIT-purified Xenium data colored by purification status")


(p1|p2) 
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/plot-raw-split-purified-1.png)<!-- -->

``` r
(p3|p4)
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/plot-raw-split-purified-2.png)<!-- -->

# Reassignment of filtered out transcripts

For this, we need to compute: - spatial network to define neighboring
cells to which the residual transcripts will be reassigned - reassign
residual counts (`raw_counts - purified_counts`) among cell’s
neighboring cells of the secondary cell type (done by
`reassign_residual_counts()`)

## Compute spatial network

``` r
sp_nw <- SPLIT::build_spatial_network(
  xe, 
  reduction = "spatial",
  dims = 1:2, 
  DO_prune = TRUE, 
  rad_pruning = 15, # remove connections further than 15um
  k_knn = 20
)
```

    ## Computing nearest neighbors

    ## Only one graph name supplied, storing nearest-neighbor graph only

    ## N = 2424044 ( 78 %) edges were pruned

``` r
sp_nw <- SPLIT::add_spatial_metric(spatial_neighborhood = sp_nw, rctd = RCTD)

sp_neigh_df <- SPLIT::neighborhood_analysis_to_metadata(sp_nw)
xe <- AddMetaData(xe, sp_neigh_df)
#rm(sp_nw, sp_neigh_df)
```

Compute reassignment

``` r
pur_genes <- rownames(xe_purified)
pur_cells <- colnames(xe_purified)
pur_counts <- GetAssayData(xe_purified, layer = "counts", assay = "Xenium") #[pur_genes, pur_cells]
raw_counts <- GetAssayData(xe, layer = "counts", assay = "Xenium") #[pur_genes, pur_cells]


reassignmet_results <- SPLIT::reassign_residual_counts(
  raw_counts = raw_counts, 
  corrected_counts = pur_counts, 
  spatial_network = sp_nw, 
  purification_status = xe_purified$purification_status,
  mode = "count_proportional", # "uniform" for uniform distribution among neighbors of contaminating cell type, or "count_proportional" for the reassignment proportional to contaminating neighbors' number of counts
  return_reassignment_operator = TRUE # optional, to explore which cells effectively send transcripts and which cells receive (see details below)
)

reassigned_counts <- reassignmet_results$corrected_counts
reassignment_operator <- reassignmet_results$reassignment_operator

#rm(reassignmet_results)
```

## Evaluate transcript re-assigment

We can now explore which cells sent transcripts and which cells received
them upon re-assignment.

``` r
reassignment_operator <- reassignment_operator[pur_cells, pur_cells]
cells_that_uneffectively_send_transccripts <- rownames(reassignment_operator)[rowSums(reassignment_operator) == 0]

receiving <- colSums(reassigned_counts)
sending <- colSums(raw_counts[pur_genes, pur_cells] - pur_counts[pur_genes, pur_cells])
sending_effectively <- sending 
sending_effectively[cells_that_uneffectively_send_transccripts] <- 0

sending_receiving <- data.frame(nCount_sent = sending, nCount_sent_effectively = sending_effectively, nCount_received = receiving)

plot(sending_receiving$nCount_sent, sending_receiving$nCount_sent_effectively)
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
cell.meta <- xe_purified@meta.data
xe_purified_balanced_score <- CreateSeuratObject(
  counts = pur_counts,
  assay = "Xenium"
)

xe_post_redistribution <- CreateSeuratObject(
  counts = reassigned_counts[pur_genes, pur_cells], #  no need to add purified counts as they have been already accounted in `reassigned_counts`
  meta.data = res_split$cell_meta,
  assay = "Xenium"
)

xe_post_redistribution <- AddMetaData(xe_post_redistribution, sending_receiving)

xe_post_redistribution$nCount_Xenium
```

    ##          1          2          8          9         10         11         14         15         16         17         18         20         23         26         27 
    ##  26.000000  37.840823  43.923451  11.197165  61.380218  39.272797  24.825798 103.621569  50.527331   7.930965  49.965610  12.230983  42.968316  26.132306  13.011576 
    ##         30         32         35         39         40         42         43         44         45         46         47         48         49         53         56 
    ##  48.491993  18.498101  19.586945  47.387264   7.698639  16.301487  55.903330  35.202101  13.435353  15.072056  35.000000  14.052568  23.579705  38.918512  44.418019 
    ##         57         61         65         66         67         68         69         70         71         72         73         74         75         76         77 
    ##  73.592892  36.725213  20.844868  38.482130  17.863005 149.079358 216.744259 105.726632  75.364887  97.509582  87.044351  67.233153 195.000000  45.204063 122.507556 
    ##         78         79         80         81         82         83         84         85         86         87         88         89         90         91         92 
    ##  87.508437  47.808918  60.590797 129.000000  96.434565  39.844596 185.775489 125.596246 138.794653 101.424106  42.775661 119.219207 134.860623  81.666493 146.560519 
    ##         93         94         95         96         97         98         99        100        101        102        103        104        105        106        107 
    ## 102.123416  31.162442 135.316753 181.206571  74.063706  92.673041  75.798360  75.781177 127.842279  84.135924 111.000000  81.569780 142.000000 158.653733 360.234542 
    ##        108        109        110        111        112        113        114        115        116        117        118        119        120        121        122 
    ##  97.704432 173.547417  61.254006 132.756377 144.125581 260.000000 172.104053  99.111578  73.540513 140.681654 357.609842 143.438113 141.514391 505.696625 112.192657 
    ##        123        124        125        126        127        128        129        130        131        132        133        134        135        136        137 
    ##  93.588892 376.445264 101.060722 296.692488  98.696523  46.000000 141.305920 257.000000 275.476125  23.135666 276.968012 239.236385 182.527385 135.286643 222.034981 
    ##        138        139        140        141        142        143        144        145        146        147        148        149        150        151        152 
    ## 208.207197 242.550305 246.232341 189.498077 237.979374 146.249962 625.490910 362.373102 340.698393 156.195678 331.972021  52.000000 105.783971 390.000000 204.505165 
    ##        153        154        155        156        157        158        159        160        161        162        163        164        165        166        167 
    ## 259.004871 301.515738 241.223506 103.039379 170.000000 276.000000 323.642933 218.137585 227.000000 399.000000 195.133555 112.685730 419.000000  93.000000 424.000000 
    ##        168        169        170        171        172        173        174        175        176        177        178        180        182        183        184 
    ##  92.017814 134.769681  84.628541 263.359431 115.563297 248.723333 187.000000  84.238650  12.549253 195.120394 485.694744 110.824510 111.372944  23.066955 134.000000 
    ##        185        186        188        189        191        192        193        195        196        197        198        199        201        202        203 
    ##  41.389732 212.000000 204.582919  24.909142  14.309316  45.442352   9.942865 178.256397  28.473766  76.267766 258.216992 224.000000 328.156363 262.928854  96.951178 
    ##        204        205        206        207        208        209        210        211        212        213        214        215        216        217        218 
    ## 180.000000  77.000000 182.064767 271.000000  88.727009 150.000000 167.613641 105.919530  90.000000 148.000000  97.846713 103.000000  63.000000  96.776923 109.828626 
    ##        219        220        221        222        223        224        225        226        227        228        229        230        231        232        233 
    ##  99.000000  80.000000 118.224876  64.044171 199.291033 238.000000 112.000000 268.000000  79.663337 141.964202 148.000000 214.000000 204.000000 144.890281  61.331023 
    ##        234        235        236        237        238        239        240        241        242        243        244        245        246        247        248 
    ## 354.505414 167.000000  52.294831 333.397957  96.883326 278.000000 265.844554  33.356377 277.000000 102.687521 392.673470  76.495428 192.000215 202.217635 155.028511 
    ##        249        250        251        252        253        254        255        256        257        258        259        260        261        262        263 
    ##  17.257505 167.184214 108.172169 156.000000  44.497309 199.659256 134.656063 214.656092 262.000000 364.344355 313.214260 260.000000 130.794467 160.538998  23.935644 
    ##        264        265        266        267        268        269        270        271        272        273        274        275        276        277        278 
    ## 192.000000 110.540622 135.440392 249.274948 218.767108 284.000000  76.168282  91.036806  53.320371 150.828527 283.757772 241.150618 157.855720  51.299793 165.154641 
    ##        279        280        282        283        284        285        286        287        288        289        290        291        292        293        294 
    ## 221.145005  83.379618  87.897677  97.229868 285.943114  95.843475 296.000000 232.717553 143.468829 177.811640 178.105897  49.672136 515.221967 266.866827  69.219771 
    ##        295        296        297        298        299        300        301        302        303        304        305        306        308        309        310 
    ##  53.986187 160.384190  95.376066  26.962950  63.563800  72.750791  28.913744  93.549062 198.073464 142.929470  23.536288  88.417332 123.000000 165.648785 245.409837 
    ##        311        312        313        314        315        316        317        318        320        321        322        323        324        325        326 
    ## 260.424402  85.331747 238.512832 112.202365 287.391479  96.354546 177.000000 305.137799  51.045445 611.554296 100.060226 127.488362  35.238247  27.532422  83.768602 
    ##        327        328        329        330        331        332        333        334        335        336        337        338        339        340        341 
    ## 226.666203  46.137382 116.000000 130.544458  57.259912 386.458637 233.322797 162.354496 417.000000 153.593313  57.298594 117.506577 187.489155 107.350651 100.861124 
    ##        342        343        344        345        346        347        348        349        350        351        352        353        354        355        356 
    ## 148.463363 183.000000 435.582668 138.413997  89.422297 124.058565  47.306978 514.693022  92.240091  45.063778  51.041218  56.864803  41.044538 118.060853 137.550156 
    ##        357        358        359        360        361        362        363        364        365        366        367        368        369        370        371 
    ##  56.785283  37.978901  84.558552  82.860258  94.057361  77.336771  39.003285  77.188848  65.879959 122.635068  75.221956  58.222803  77.249387  44.912574  41.625185 
    ##        372        373        374        375        376        377        378        379        380        381        382        383        384        385        386 
    ##  79.813001 119.822985  72.822695 209.821147  47.241499 126.887892  94.428526  29.377330 111.419872 156.231208  39.941244  45.185014  85.135693  28.381931  60.631934 
    ##        387        388        389        390        391        392        393        394        395        396        397        399        400        402        403 
    ## 118.394981  45.448051  51.582852  30.352939  46.538911 177.000000 192.912635  79.445505  16.355081  81.471894  93.519312  22.261456  94.561591  38.898442  63.806498 
    ##        404        405        406        407        408        409        410        411        412        413        414        415        416        417        418 
    ##  20.841090  23.012098  80.434096 126.987519 103.000000  45.375125 111.000000 142.000000 214.000000  99.000000  79.000000 376.000000 356.000000 403.000000 209.868933 
    ##        419        420        421        422        423        424        425        426        427        428        429        430        431        432        433 
    ## 186.000000 240.313980 261.513423 358.966822 108.472214 358.772480 488.000000 230.420460 285.155917 269.799877  52.694383 165.617645 227.794677 276.264774 134.547390 
    ##        434        435        436        437        438        439        440        441        443        444        446        447        448        449        450 
    ## 311.675879 245.205784 202.787334 316.455568 210.536236  22.325859  96.693328  24.827038  87.000000  94.463654 115.546154 330.854688 452.638227 197.195659 517.907394 
    ##        451        452        453        454        455        456        457        458        459        460        461        462        463        464        465 
    ## 237.046712 332.861550 207.074372 154.430628 333.927297 100.582497 153.002846 133.154360 142.209439  73.072703 248.171705 366.567563 136.349595  92.019414 115.969636 
    ##        466        467        468        469        470        471        472        473        474        475        477        478        479        480        481 
    ## 129.184980  86.637372 107.435437 106.206374 217.181458 122.446118 222.566011  86.226467 169.312499 255.921865 307.972249  20.011343 182.433370  42.000000  58.433644 
    ##        482        483        484        485        486        488        489        490        491        492        493        494        495        496        497 
    ##  33.000000  74.000000  16.011179  45.103741  29.825705  44.249204 176.825477 149.137413 129.522051  93.640108 191.207391  81.850313 132.438680  60.821728 126.726639 
    ##        498        499        500        501        502        503        504        505        506        507        508        509        511        512        513 
    ##  98.422107 105.094098  73.314272  77.485581 118.000000  97.392924 130.186589  85.258235 131.955296 145.656332 152.804633  82.778769  46.784720 195.359259  59.391070 
    ##        514        515        516        517        518        519        520        521        522        523        524        525        526        527        528 
    ##  74.724240  83.871695 153.399417 199.899399  36.479729 108.759400  50.473905  63.118283 125.040039 105.540230  61.576185 110.280912 356.395621  28.445030 106.284104 
    ##        529        530        531        532        533        534        535        536        537        538        539        540        541        542        543 
    ## 197.000000 167.557928  91.364087  43.529483 134.051099  69.483372  60.000000  16.384662  67.227766  48.806222 135.745085  37.000000  21.682092  37.053817 187.252259 
    ##        544        545        546        547        548        549        550        551        552        553        554        555        556        557        558 
    ##  66.884810 119.333092 140.753028 126.179204 139.205271  73.616137  66.033020  97.201506  66.424207  77.000000  75.230211  45.197046 274.106553 117.893447 164.632178 
    ##        559        560        561        562        563        564        565        566        567        568        569        570        571        572        574 
    ##  97.905652 107.687258 118.439408  53.816535 113.767975  13.147827 190.873725  60.912937 162.431045 117.591424  97.373725 110.507453  55.596080  66.388874 436.672052 
    ##        575        576        577        578        579        581        582        583        584        585        586        587        588        589        590 
    ##  26.385679 237.731800 115.319173  59.711885  43.721161  35.708921 199.449535 149.822650 137.153646  33.533333  33.243109 216.445866  22.777126  50.957285  71.506114 
    ##        592        593        594        595        596        597        598        601        603        604        605        606        608        609        610 
    ##  97.620684  31.841298  95.369594  70.294994  35.763399  81.313346  73.095040 121.470894  53.475170  59.895344  15.680013  36.856808  74.088215  78.000000 116.741334 
    ##        611        612        613        614        616        617        618        619        620        621        622        623        624        625        626 
    ##  61.000000 141.000000 119.662736 136.750475  24.842153  94.171769  23.432935   9.824396 102.423812  11.497239  59.878465 121.147249 124.516413  39.000000  80.149895 
    ##        627        628        629        630        631        632        633        634        635        638        641        642        643        644        645 
    ##  85.287852  70.717697 125.814218 157.790627 359.000000 229.685770  45.376245  74.000000  59.621422 100.864142  68.878478 123.963095 201.055098  55.251946  83.281033 
    ##        646        647        648        649        650        651        652        653        654        655        656        657        658        659        660 
    ##  82.000000 154.795451  61.026438 173.264091  73.852350 153.654306  26.288655  31.121881  71.401704 109.671106  38.290379  86.796375  60.645281 111.477216  15.000000 
    ##        661        662        663        664        665        666        667        668        669        670        671        672        673        675        676 
    ##  78.151925  71.268011  57.289255  34.367527  46.573552  40.528468 123.457507 161.467692  50.013498 159.577742  35.000000  43.579539  45.457806 131.205995 193.117264 
    ##        678        679        680        681        682        683        684        685        686        687        688        689        690        691        692 
    ## 201.463361  86.255828 114.478316 170.119547 161.372654 265.153880  99.704631  84.188992 138.927022  31.000000  63.385607  23.545271  90.566291 111.478427  43.725270 
    ##        693        694        695        696        697        698        699        700        701        702        703        704        705        706        707 
    ## 109.620385 249.303196 189.293855 104.595043 157.845446  72.906237  93.800562 155.707513  49.843993 165.895172 160.672438  31.137045  49.049375 194.181835  69.070250 
    ##        708        709        710        711        712        713        714        715        716        717        718        719        720        721        722 
    ##  50.450276  81.093398  69.388516  84.683409  68.051458 101.245671  39.780666 131.584429 107.605403 132.838304  98.283530  52.394925 126.482065  88.000000 129.593295 
    ##        723        724        725        726        727        728        730        731        732        733        734        735        736        737        738 
    ##  96.599907 141.360225 182.567998 159.000000 143.398429 160.236642  68.676272  92.925298  78.340390 161.549881 332.493567 118.775717 134.561337 148.117704 144.180587 
    ##        739        740        741        742        743        744        745        746        747        748        749        750        751        752        753 
    ## 159.378316 148.546101 262.847686 101.234643  46.925821 138.546816 198.897908  92.126703  96.988820 215.032788  97.094275  78.237895  13.468232  55.040732  89.265248 
    ##        754        755        756        757        758        759        760        761        762        763        764        765        766        767        768 
    ## 307.974922  60.959970 254.895000  41.237246 306.655897 114.082521 131.574063 154.014730 182.535992 406.297696 153.447693  60.951114 289.016315 572.000000 300.212702 
    ##        769        770        771        772        773        774        775        776        777        778        779        780        781        782        783 
    ## 274.000000 215.027223 158.000000 258.836879 267.000000 176.000000 275.475093 183.262254 289.000000 295.966452  82.921412 248.736808  91.000000 160.000000 329.199189 
    ##        784        785        786        787        788        789        790        791        792        793        794        795        796        797        798 
    ## 110.456926 170.053505 119.782997  94.000000 105.287808 147.000000  14.718077 115.924386 187.707907  60.171027  74.521827  74.221809  26.000000  72.000000  50.000000 
    ##        799        800        801        802        803        804        805        806        807        808        809        810        811        812        813 
    ##  68.021313  66.349656  47.240476 123.184574  88.888230  53.921471  78.390217 106.845013 121.713481  81.875936  88.755700  58.851296  99.541237  89.682012 175.753252 
    ##        814        815        816        817        818        819        820        821        822        823        824        825        826        827        829 
    ## 284.515149 127.016681  83.187129 198.998655  13.366547 155.232802  27.201926 164.751930  70.743562  56.528537 278.994502  84.902783  33.095403 110.717928 125.000000 
    ##        830        831        832        833        834        835        836        837        838        839        840        841        842        843        844 
    ##  84.987407  90.419119 171.000000  43.172863  75.768388 260.187617 109.000000  91.179927 162.380295  46.176504 190.570845 107.394720 101.029623 110.735745 128.555261 
    ##        845        846        847        848        849        850        851        852        853        854        855        856        857        858        859 
    ## 120.643482  58.715570  63.240469 210.111770 159.061214  86.766637 108.112853  55.708567  73.184461  65.723509  94.313516  81.749781  86.035168  77.664083  50.622291 
    ##        860        861        862        863        864        865        866        867        868        869        870        871        872        873        874 
    ##  58.182250 285.112779 300.697671  92.001626 238.280490 100.862303 110.921976 173.211651 273.401574 386.173515  91.201617 101.682996 189.292487  64.042412  64.000000 
    ##        875        876        877        879        880        881        882        883        884        885        886        887        888        889        890 
    ##  77.000000 277.840490 213.823496 186.799410  92.000000  27.346620 166.017925  82.325315 218.946055 227.000000  85.755015 109.131953  66.330845  82.689158 101.240552 
    ##        891        892        895        896        897        898        899        901        902        904        905        906        907        908        909 
    ##  82.753854  71.799979 204.770990 238.397044 425.440089  19.547919 120.811268 343.859798  66.150592  87.000000  61.191345  58.954963  51.431092  67.000000 125.517190 
    ##        910        911        912        913        914        915        916        918        919        920        921        922        923        924        925 
    ##  90.959862  64.269321  21.905526  28.173132  67.306827  28.036220  52.111173 102.889705  13.436796 137.793816 130.579798  21.000000  51.087580  72.000000  23.442042 
    ##        926        927        928        929        930        931        932        933        934        935        936        937        938        939        940 
    ## 190.921892 134.768705  91.087024 101.293070 203.050924 171.145622 168.658364 328.015694 270.068524 324.528912 194.836506 271.372681 341.000000 305.000000 344.066468 
    ##        941        942        943        944        945        946        947        948        949        950        951        952        953        954        955 
    ## 245.193763  24.454865  80.001331 149.563824  68.508936 143.748413 121.415654  57.243366 127.032060 126.208659  46.498885 156.000000 240.000000 143.062266 146.000000 
    ##        956        957        958        959        960        961        962        963        964        965        966        967        968        969        970 
    ## 103.679513  54.312588  27.772109 116.000000  46.088343 423.757051 126.391312 152.204748 122.000000 116.461030  79.000000 223.442983  49.000000 155.937563 170.797546 
    ##        971        973        974        975        976        977        978        979        980        981        982        983        984        985        986 
    ## 351.000000  76.804749 138.199471 187.714591  85.049979 195.455095 108.244183  47.265736 147.606450 129.939790 150.000000 186.000000 139.000000 225.000000 130.000000 
    ##        987        988        990        991        992        993        994        995        996        997        998        999       1000       1001       1002 
    ##  72.000000 111.162389  75.078137  99.120114  50.244289 255.733931 286.965909  77.332155  88.009418 162.615089 277.757649 247.714431 183.000000 158.393550 101.619977 
    ##       1003       1004       1005       1006       1007       1008       1009       1010       1011       1012       1014       1015       1016       1017       1018 
    ## 334.000000 121.265619 153.204855 176.756353 147.033997 271.069630 513.179869 401.817812  99.040557 118.934536 262.586419  87.596575 338.926256 138.899587 162.165294 
    ##       1019       1020       1021       1022       1023       1025       1026       1027       1028       1029       1030       1031       1032       1033       1034 
    ## 110.807957  95.614274  61.464891  55.473425  17.248604 150.713553 150.786030  51.496588  65.552761  27.494671 128.824613  67.279659 177.936869  27.201021  93.374587 
    ##       1035       1036       1037       1038       1039       1040       1041       1042       1043       1044       1045       1046       1047       1048       1049 
    ##  23.175051  65.422964  66.711500  99.788301 104.911451  58.019582  79.448644  36.049210 101.261757  33.634181  46.081043  58.168310  41.739614  35.608983  11.285010 
    ##       1050       1051       1052       1053       1054       1055       1056       1057       1058       1059       1060       1061       1064       1065       1067 
    ##  28.999981  18.424519  43.790397  25.633732  24.185646  38.715349  35.664776  71.783370  21.824296  35.275818  13.406364  35.773395  16.528016  55.898390  26.792074 
    ##       1068       1070       1072       1073       1074       1076       1079       1080       1081       1082 
    ##  23.013005  14.084435  23.607497  67.793197  53.978986  32.487972  84.239556  52.839983  51.970844  86.780737 
    ##  [ reached 'max' / getOption("max.print") -- omitted 160420 entries ]

``` r
xe_purified_balanced_score$nCount_Xenium
```

    ##          1          2          8          9         10         11         14         15         16         17         18         20         23         26         27 
    ##  26.000000  37.840823  38.000000  11.197165  61.380218  27.893810  24.825798  89.552534  33.458851   7.930965  23.654019  12.230983  22.473773  16.713080  13.011576 
    ##         30         32         35         39         40         42         43         44         45         46         47         48         49         53         56 
    ##  48.491993  18.498101   6.888432  47.387264   7.698639  16.301487  55.903330  13.603389   8.391121  13.334448  35.000000  14.052568  17.688409  12.149495   9.010911 
    ##         57         61         65         66         67         68         69         70         71         72         73         74         75         76         77 
    ##  73.592892  26.107001  20.844868  38.482130  17.863005 149.079358 216.744259  89.000000  45.499107  97.509582  87.044351  58.277504 195.000000  45.204063 116.177376 
    ##         78         79         80         81         82         83         84         85         86         87         88         89         90         91         92 
    ##  81.499353  47.808918  57.563319 129.000000  91.893349  39.844596 185.775489 120.000000 132.372731  95.000000  42.775661  53.000000 117.327900  69.430117 146.560519 
    ##         93         94         95         96         97         98         99        100        101        102        103        104        105        106        107 
    ##  92.204904  22.683693 120.589265 181.206571  74.063706  92.673041  75.798360  75.781177 127.842279  84.135924 111.000000  81.569780 142.000000 135.111812 360.234542 
    ##        108        109        110        111        112        113        114        115        116        117        118        119        120        121        122 
    ##  97.704432 173.547417  61.254006 132.756377 144.125581 260.000000 172.104053  99.111578  73.540513 122.000000  58.961021 143.438113 141.514391 228.296780 112.192657 
    ##        123        124        125        126        127        128        129        130        131        132        133        134        135        136        137 
    ##  93.588892 376.445264 101.060722 296.692488  98.696523  46.000000 141.305920 257.000000 275.476125  23.135666 276.968012 239.236385 178.350792  76.644213 174.625136 
    ##        138        139        140        141        142        143        144        145        146        147        148        149        150        151        152 
    ## 203.459277 197.040963 208.208393 160.750811 183.205337 146.249962 555.000000 362.373102 340.698393 156.195678 331.972021  52.000000 105.783971 390.000000 204.505165 
    ##        153        154        155        156        157        158        159        160        161        162        163        164        165        166        167 
    ## 259.004871 301.515738 241.223506 103.039379 170.000000 276.000000 323.642933 198.000000 227.000000 399.000000 195.133555  99.000000 419.000000  93.000000 424.000000 
    ##        168        169        170        171        172        173        174        175        176        177        178        180        182        183        184 
    ##  78.332085 134.769681  84.628541 236.000000 115.563297 230.268538 187.000000  84.238650  12.549253 195.120394 427.000000  91.000000 111.372944  17.870993 134.000000 
    ##        185        186        188        189        191        192        193        195        196        197        198        199        201        202        203 
    ##  41.389732 212.000000 170.704767  24.909142  12.000000  43.847151   9.942865 178.256397  28.473766  76.267766 206.959596 224.000000 225.274495 231.000000  87.708023 
    ##        204        205        206        207        208        209        210        211        212        213        214        215        216        217        218 
    ## 180.000000  77.000000 182.064767 271.000000  88.727009 150.000000 167.613641 105.919530  90.000000 148.000000  97.846713 103.000000  63.000000  96.776923 109.828626 
    ##        219        220        221        222        223        224        225        226        227        228        229        230        231        232        233 
    ##  99.000000  80.000000 118.224876  64.044171 176.472587 238.000000 112.000000 268.000000  79.013014 141.964202 148.000000 214.000000 204.000000 144.890281  61.331023 
    ##        234        235        236        237        238        239        240        241        242        243        244        245        246        247        248 
    ## 354.505414 167.000000  52.294831 299.000000  96.883326 278.000000 212.532757  33.356377 277.000000 102.687521 331.092188  53.968015 151.532515 202.217635 133.572360 
    ##        249        250        251        252        253        254        255        256        257        258        259        260        261        262        263 
    ##  17.257505 167.184214  46.356383 156.000000  44.497309 170.003164 134.656063 185.000000 262.000000 364.344355 313.214260 260.000000 130.794467 110.749589  23.935644 
    ##        264        265        266        267        268        269        270        271        272        273        274        275        276        277        278 
    ## 192.000000 110.540622 135.440392 177.000000 218.767108 284.000000  76.168282  91.036806  53.320371 150.828527 283.757772 241.150618 157.855720  51.299793 165.154641 
    ##        279        280        282        283        284        285        286        287        288        289        290        291        292        293        294 
    ## 191.221198  83.379618  66.127544  97.229868  70.000000  95.843475 296.000000 232.717553 143.468829 177.811640 154.000000  49.672136 444.000000 196.909984  69.219771 
    ##        295        296        297        298        299        300        301        302        303        304        305        306        308        309        310 
    ##  24.432701  36.703662  95.376066  26.962950  54.213304  72.750791  28.913744  93.549062 198.073464 142.929470  23.536288  88.417332 123.000000 165.648785 245.409837 
    ##        311        312        313        314        315        316        317        318        320        321        322        323        324        325        326 
    ## 260.424402  65.913947 238.512832  97.708942 287.391479  96.354546 177.000000 269.000000  51.045445 611.554296  82.794072 127.488362  35.238247  27.532422  83.768602 
    ##        327        328        329        330        331        332        333        334        335        336        337        338        339        340        341 
    ## 138.527328  46.137382 116.000000  90.244691  46.942547 168.422770 233.322797 111.368361 417.000000  65.448155  18.424179 117.506577 159.000000 107.350651 100.861124 
    ##        342        343        344        345        346        347        348        349        350        351        352        353        354        355        356 
    ## 148.463363 183.000000 246.000000 138.413997  89.422297 124.058565  47.306978 486.000000  92.240091  37.546970  51.041218  37.555518  35.021224  72.000000  53.493774 
    ##        357        358        359        360        361        362        363        364        365        366        367        368        369        370        371 
    ##  41.443659  37.978901  76.025525  77.709883  42.858269  54.059825  39.003285  57.305301  65.879959  81.050764  70.898198  23.953827  77.249387  44.912574  41.625185 
    ##        372        373        374        375        376        377        378        379        380        381        382        383        384        385        386 
    ##  52.106916  45.929901  60.508421 162.000000  38.731024 107.178853  81.846068  29.377330 105.824020 141.290961  39.941244  45.185014  80.222262  28.381931  60.631934 
    ##        387        388        389        390        391        392        393        394        395        396        397        399        400        402        403 
    ##  97.161879  45.448051  51.582852  30.352939  46.538911 177.000000 147.000000  58.000000  16.355081  46.377943  61.000000  16.947813  94.561591  38.898442  63.806498 
    ##        404        405        406        407        408        409        410        411        412        413        414        415        416        417        418 
    ##  20.841090  23.012098  80.434096 126.987519 103.000000  45.375125 111.000000 142.000000 214.000000  99.000000  79.000000 376.000000 356.000000 403.000000 177.900880 
    ##        419        420        421        422        423        424        425        426        427        428        429        430        431        432        433 
    ## 186.000000 240.313980 261.513423 311.031947 101.791007 338.000000 488.000000 230.420460 285.155917 223.183911  52.694383  62.213575 227.794677 208.201473  39.502458 
    ##        434        435        436        437        438        439        440        441        443        444        446        447        448        449        450 
    ## 311.675879 224.845111 184.753595 316.455568 177.674169  22.325859  96.693328  22.955767  87.000000  94.463654 115.546154 164.427569 282.522023 165.000000 274.814695 
    ##        451        452        453        454        455        456        457        458        459        460        461        462        463        464        465 
    ## 144.397612 332.861550 207.074372 122.054739 313.000000 100.582497 153.002846 105.613631 124.930604  73.072703 248.171705 163.794030 105.015514  92.019414 115.969636 
    ##        466        467        468        469        470        471        472        473        474        475        477        478        479        480        481 
    ## 129.184980  86.637372  88.332850 106.206374 217.181458  80.428661 222.566011  86.226467 169.312499 255.921865 291.000000  20.011343 163.838002  42.000000  58.433644 
    ##        482        483        484        485        486        488        489        490        491        492        493        494        495        496        497 
    ##  33.000000  74.000000  16.011179  45.103741  29.825705  44.249204 176.825477 108.484942 129.522051  93.640108 191.207391  57.863537 132.438680  41.013224  61.771259 
    ##        498        499        500        501        502        503        504        505        506        507        508        509        511        512        513 
    ##  98.422107  66.084027  53.000000  64.961123 118.000000  67.475542 118.093572  85.258235 100.249819 145.656332 135.529494  82.778769  46.784720 150.354611  59.391070 
    ##        514        515        516        517        518        519        520        521        522        523        524        525        526        527        528 
    ##  32.995352  83.871695 153.399417 169.030471  36.479729  93.000000  50.473905  55.423837 110.310671  73.797354  61.576185 110.280912 258.840651  28.445030  83.065528 
    ##        529        530        531        532        533        534        535        536        537        538        539        540        541        542        543 
    ## 197.000000 145.105615  91.364087  35.146512 134.051099  64.743073  60.000000  16.384662  36.992880  48.806222 110.781932  37.000000  21.682092  37.053817 149.000000 
    ##        544        545        546        547        548        549        550        551        552        553        554        555        556        557        558 
    ##  66.884810 119.333092  47.346808  93.036847 122.406777  73.616137  66.033020  97.201506  66.424207  77.000000  75.230211  45.197046 223.107309  97.000756 155.201537 
    ##        559        560        561        562        563        564        565        566        567        568        569        570        571        572        574 
    ##  97.905652 107.687258 116.569359  50.750481 111.942100  13.147827 153.000000  60.000000 158.186466 115.588852  97.373725  78.979673  55.596080  26.663391 139.550294 
    ##        575        576        577        578        579        581        582        583        584        585        586        587        588        589        590 
    ##  26.385679 237.731800 115.319173  59.711885  43.721161  35.708921 199.449535 149.822650 137.153646  33.533333  33.243109 212.153057  22.777126  50.957285  71.506114 
    ##        592        593        594        595        596        597        598        601        603        604        605        606        608        609        610 
    ##  65.034446  31.841298  60.610940  70.294994  35.763399  81.313346  52.455950 121.470894  53.475170  41.005329  15.680013  36.856808  74.088215  78.000000 116.741334 
    ##        611        612        613        614        616        617        618        619        620        621        622        623        624        625        626 
    ##  61.000000 141.000000  79.001162 115.832448  24.842153  79.073043  23.432935   9.824396  86.853251  11.497239  59.878465 121.147249  95.542126  39.000000  70.309130 
    ##        627        628        629        630        631        632        633        634        635        638        641        642        643        644        645 
    ##  50.596982  70.717697  88.700391 157.790627 359.000000  90.431052  45.376245  74.000000  59.621422 100.864142  68.878478 123.963095 201.055098  42.854084  77.000000 
    ##        646        647        648        649        650        651        652        653        654        655        656        657        658        659        660 
    ##  82.000000  93.012070  55.550665 129.276161  73.852350  81.944684  26.288655  26.066565  46.093374  62.426679  38.290379  86.796375  60.645281  74.629583  15.000000 
    ##        661        662        663        664        665        666        667        668        669        670        671        672        673        675        676 
    ##  78.151925  71.268011  57.289255  30.806143  46.573552  40.528468  54.794393  68.456546  50.013498 118.886093  35.000000  43.579539  45.457806  86.164400 172.178620 
    ##        678        679        680        681        682        683        684        685        686        687        688        689        690        691        692 
    ## 201.463361  86.255828 114.478316  34.408045 161.372654 265.153880  99.704631  67.404951  71.478280  31.000000  63.385607  19.490899  72.828611  95.918406  43.725270 
    ##        693        694        695        696        697        698        699        700        701        702        703        704        705        706        707 
    ##  56.518624 157.898239 189.293855 104.595043 157.845446  72.906237  93.800562  98.068268  49.843993 165.895172 153.929108  15.495911  31.358866 194.181835  69.070250 
    ##        708        709        710        711        712        713        714        715        716        717        718        719        720        721        722 
    ##  50.450276  81.093398  50.224879  84.683409  49.114997 101.245671  39.780666  86.515979  63.864439 113.228513  77.366741  50.900402 111.595088  88.000000  92.037129 
    ##        723        724        725        726        727        728        730        731        732        733        734        735        736        737        738 
    ##  86.008401  91.444218 125.521133 159.000000 103.412533 160.236642  36.357193  65.268387  51.920101  84.988311 293.246751 103.010844 127.824128  59.668524 144.180587 
    ##        739        740        741        742        743        744        745        746        747        748        749        750        751        752        753 
    ##  82.845051  86.247724 227.955877  68.920197  23.735931 106.735759 156.000000  50.882368  79.129374 138.390070  88.226918  62.000000  13.468232  41.429813  57.677359 
    ##        754        755        756        757        758        759        760        761        762        763        764        765        766        767        768 
    ## 293.420464  58.135971 207.712127  36.327832 243.150398  84.593226 115.362220  43.026927 179.468005 361.794596  74.876031  60.951114 284.960333 572.000000 300.212702 
    ##        769        770        771        772        773        774        775        776        777        778        779        780        781        782        783 
    ## 274.000000 215.027223 158.000000 258.836879 267.000000 176.000000 275.475093 183.262254 289.000000 295.966452  82.921412 248.736808  91.000000 160.000000 301.547408 
    ##        784        785        786        787        788        789        790        791        792        793        794        795        796        797        798 
    ##  71.535047  98.444482  94.348219  94.000000  99.390977 147.000000  14.718077 109.446213  75.412900  60.171027  74.521827  74.221809  26.000000  72.000000  50.000000 
    ##        799        800        801        802        803        804        805        806        807        808        809        810        811        812        813 
    ##  68.021313  66.349656  47.240476  86.425050  88.888230  44.427439  45.817656  66.173523  85.525918  59.049518  88.755700  58.851296  85.994470  49.877365 175.753252 
    ##        814        815        816        817        818        819        820        821        822        823        824        825        826        827        829 
    ## 242.630212 112.641535  48.855622 198.998655  13.366547 155.232802  27.201926 140.870639  70.743562  56.528537 258.000000  78.104709  28.526137 110.717928 125.000000 
    ##        830        831        832        833        834        835        836        837        838        839        840        841        842        843        844 
    ##  20.858046  78.972414 171.000000  43.172863  75.768388 260.187617 109.000000  91.179927 162.380295  46.176504 190.570845 102.774968  15.380248  72.350625 128.555261 
    ##        845        846        847        848        849        850        851        852        853        854        855        856        857        858        859 
    ##  62.198743  37.423681  63.240469 122.000000 159.061214  86.766637 108.112853  55.708567  53.470835  65.723509  94.313516  81.749781  86.035168  77.664083  50.622291 
    ##        860        861        862        863        864        865        866        867        868        869        870        871        872        873        874 
    ##  58.182250 285.112779 256.680411  74.452289 204.345140 100.862303  37.064650 173.211651 202.443986 328.653858  91.201617 101.682996  97.360755  64.042412  64.000000 
    ##        875        876        877        879        880        881        882        883        884        885        886        887        888        889        890 
    ##  77.000000 230.000000  77.000000  72.947431  92.000000  27.346620 143.177697  82.325315 218.946055 227.000000  85.755015  29.217814  66.330845  82.689158  42.000000 
    ##        891        892        895        896        897        898        899        901        902        904        905        906        907        908        909 
    ##  82.753854  71.799979  96.445544 218.550738 152.851753  19.547919 110.651850 308.000000  43.275352  87.000000  59.182312  58.954963  50.000000  67.000000  58.000000 
    ##        910        911        912        913        914        915        916        918        919        920        921        922        923        924        925 
    ##  88.097677  62.287808  21.905526  25.239919  60.527847  23.000000  43.280755 100.000000  12.996460 123.000000 130.579798  21.000000  40.877575  72.000000  23.442042 
    ##        926        927        928        929        930        931        932        933        934        935        936        937        938        939        940 
    ## 142.000000 134.768705  91.087024 101.293070 203.050924 171.145622 168.658364 328.015694 270.068524 324.528912 194.836506 271.372681 341.000000 305.000000 344.066468 
    ##        941        942        943        944        945        946        947        948        949        950        951        952        953        954        955 
    ## 245.193763  24.454865  80.001331 131.855428  31.159772 143.748413  70.783261  39.399805 101.387942 103.919772  34.138684 156.000000 240.000000 143.062266 146.000000 
    ##        956        957        958        959        960        961        962        963        964        965        966        967        968        969        970 
    ## 103.679513  54.312588  27.772109 116.000000  44.674397 423.757051  71.720694 152.204748 122.000000 116.461030  79.000000 178.543817  49.000000 155.937563 167.624024 
    ##        971        973        974        975        976        977        978        979        980        981        982        983        984        985        986 
    ## 351.000000  76.804749 136.000000 185.000000  70.972130 182.741522 108.244183  47.265736  86.000000 129.939790 150.000000 186.000000 139.000000 225.000000 130.000000 
    ##        987        988        990        991        992        993        994        995        996        997        998        999       1000       1001       1002 
    ##  72.000000 111.162389  71.000000  99.120114  46.537395 255.733931 286.965909  34.087469  88.009418 128.479675 261.157908 247.714431 183.000000 158.393550  70.353923 
    ##       1003       1004       1005       1006       1007       1008       1009       1010       1011       1012       1014       1015       1016       1017       1018 
    ## 334.000000 121.265619 118.878707 176.756353 147.033997 271.069630 210.733946 401.817812  77.673852 118.934536 262.586419  87.596575 241.821171  91.179945 149.509948 
    ##       1019       1020       1021       1022       1023       1025       1026       1027       1028       1029       1030       1031       1032       1033       1034 
    ## 110.807957  63.364037  61.464891  55.473425  17.248604 139.749763 115.000000  51.496588  60.776072  24.719997  63.740545  42.582839 168.000000  18.452866  86.641983 
    ##       1035       1036       1037       1038       1039       1040       1041       1042       1043       1044       1045       1046       1047       1048       1049 
    ##  20.619155  65.422964  66.711500  90.150405  96.168502  58.019582  79.448644  36.049210  91.241657  33.634181  21.322700  54.157215  26.542752  35.608983  11.285010 
    ##       1050       1051       1052       1053       1054       1055       1056       1057       1058       1059       1060       1061       1064       1065       1067 
    ##  28.999981  17.587124  33.662050  24.656772  24.185646  25.152650  27.677537  53.412719  21.824296  29.642109  13.406364  33.737010  16.528016  24.160765  16.952839 
    ##       1068       1070       1072       1073       1074       1076       1079       1080       1081       1082 
    ##  23.013005  14.084435  23.607497  28.365605  50.103711  32.487972  77.224977  45.825405  51.970844  86.780737 
    ##  [ reached 'max' / getOption("max.print") -- omitted 160420 entries ]

# Exploting cells that cannot send their residual transcripts

What are the cells that have signal of contamination, but do not have
neighbors of the secondary cell type around?

``` r
xe@meta.data <- xe@meta.data %>%
  mutate(purified_with_no_secondary_type_sp_neighbors = 
           (rownames(xe@meta.data) %in% cells_that_uneffectively_send_transccripts) &
           (rownames(xe@meta.data) %in% rownames(res_split$cell_meta)[res_split$cell_meta$purification_status == "purified"])
         )

UMAPPlot(xe, group.by = "purified_with_no_secondary_type_sp_neighbors") + theme_void() + theme(aspect.ratio = 1)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
nb_weight_features <- c("neighborhood_weights_second_type", "neighborhood_weights_first_type")
nb_count_features <- c("second_type_neighbors_no_reject_N", "first_type_neighbors_N")
FeaturePlot(xe, features = nb_weight_features) & theme_void() & theme(aspect.ratio = 1)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
FeaturePlot(xe, features = nb_count_features) & theme_void() & theme(aspect.ratio = 1)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
FeatureScatter(xe, feature1 = nb_weight_features[1], feature2 = nb_weight_features[2], group.by = "purified_with_no_secondary_type_sp_neighbors")
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
FeatureScatter(xe, feature1 = nb_weight_features[1], feature2 = nb_weight_features[2], group.by = "first_type", cols = pal)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

``` r
VlnPlot(xe, features = nb_count_features, group.by = "purified_with_no_secondary_type_sp_neighbors", pt.size = 0)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->

``` r
df <- data.frame(
  n_counts_pur = xe_purified_balanced_score$nCount_Xenium, 
  n_counts_reassigned = xe_post_redistribution$nCount_Xenium
)
plot(df$n_counts_pur, df$n_counts_reassigned)
abline(0,1)
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
summary(df$n_counts_pur)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##    5.006   68.136  122.230  153.620  210.784 1240.000

``` r
summary(df$n_counts_reassigned)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##    5.017   79.907  140.134  170.442  231.000 1381.225

``` r
xe_post_redistribution <- xe_post_redistribution %>%
  SCTransform(assay = "Xenium", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)
```

``` r
p_pr <- UMAPPlot(xe_post_redistribution, group.by = c("first_type"), cols = pal, label = T, repel = T) + ggtitle("Fully purified + counts reassignment") & theme(aspect.ratio = 1)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

``` r
p_pr_st <- UMAPPlot(xe_post_redistribution, group.by = c("second_type"), cols = pal) + ggtitle("Second type") & theme(aspect.ratio = 1)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

``` r
p_pur <- UMAPPlot(xe_purified, group.by = "first_type", cols = pal, label = T, repel = T) + ggtitle("Fully purified") & theme(aspect.ratio = 1) & NoLegend()
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

``` r
p_pur + p_pr
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
p_pr_st
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
xe_post_redistribution$nCount_sent_uneffectively <- xe_post_redistribution$nCount_sent - xe_post_redistribution$nCount_sent_effectively


FeaturePlot(xe_post_redistribution, features = c("nCount_sent", "nCount_sent_effectively", "nCount_received")) & theme(aspect.ratio = 1)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
FeaturePlot(xe_post_redistribution, features =  c("nCount_sent", "nCount_received", "nCount_sent_uneffectively", "nCount_sent_effectively"), raster = F) &
  scale_color_gradientn(
    colors = c("gray", "blue"),
    #colours = viridis::viridis(100),
    limits = c(NA, 150)
  ) &
  theme(aspect.ratio = 1)
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
VlnPlot(xe_post_redistribution, features = c("nCount_sent", "nCount_received", "nCount_sent_uneffectively", "nCount_sent_effectively"), ncol = 2,  group.by = "first_type", cols = pal, pt.size = 0, same.y.lims = T)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
VlnPlot(xe_post_redistribution, features = c("nCount_sent", "nCount_received", "nCount_sent_uneffectively", "nCount_sent_effectively"), ncol = 2,  group.by = "second_type", cols = pal, pt.size = 0, same.y.lims = T)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
plot(xe_post_redistribution$nCount_sent_effectively, xe_post_redistribution$nCount_received)
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
total_transcripts_raw <- rowSums(GetAssayData(xe, layer = "counts", assay = "Xenium"))
total_transcripts_raw_common_genes <- total_transcripts_raw[pur_genes]
total_transcripts_purified <- rowSums(GetAssayData(xe_purified_balanced_score, layer = "counts", assay = "Xenium"))
total_transcripts_purified_and_reassigned <- rowSums(GetAssayData(xe_post_redistribution, layer = "counts", assay = "Xenium"))

df_total_transcripts <- data.frame(
  total_transcripts_raw_common_genes = total_transcripts_raw_common_genes,
  total_transcripts_purified = total_transcripts_purified,
  total_transcripts_purified_and_reassigned = total_transcripts_purified_and_reassigned
)

df_total_transcripts %>% 
  ggplot(aes(x = total_transcripts_raw_common_genes, y = total_transcripts_purified_and_reassigned)) +
  geom_point() + geom_abline() +
  scale_x_log10() + scale_y_log10() + theme_classic() + theme(aspect.ratio = 1)
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
df_plot <- data.frame(
  category = factor(
    c("Raw common genes", "Purified", "Purified + Reassigned"),
    levels = c("Raw common genes", "Purified", "Purified + Reassigned")  # enforce order
  ),
  value = c(
    sum(df_total_transcripts$total_transcripts_raw_common_genes),
    sum(df_total_transcripts$total_transcripts_purified),
    sum(df_total_transcripts$total_transcripts_purified_and_reassigned)
  )
)

ggplot(df_plot, aes(x = category, y = value, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(x = NULL, y = "Total transcripts") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = .5))
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
plot(xe_post_redistribution$nCount_sent, xe_post_redistribution$nCount_received)
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
xe_post_redistribution$nCount_sent %>% hist()
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
VlnPlot(xe_post_redistribution, features = c("nCount_sent", "nCount_received"), pt.size= 0)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`
    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

``` r
sum(xe_post_redistribution$nCount_received) - sum(xe_post_redistribution$nCount_sent)
```

    ## [1] 21002555

# Do residuals even make sense?

``` r
residual_counts <- raw_counts[pur_genes, pur_cells] -  res_split$purified_counts[pur_genes, pur_cells]
  
xe_residuals <- CreateSeuratObject(
  counts = residual_counts, 
  meta.data = res_split$cell_meta,
  assay = "Xenium"
)

xe_residuals <- subset(xe_residuals, nCount_Xenium > 1)

xe_residuals <- xe_residuals %>%
  SCTransform(assay = "Xenium", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)
```

Yes, the residuals group by the secondary (ie, residual) cell type

``` r
UMAPPlot(xe_residuals, group.by = c("first_type"), cols = pal, label = F, repel = F)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
UMAPPlot(xe_residuals, group.by = c("second_type"), cols = pal, label = T, repel = T)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
hist(xe_residuals$nCount_Xenium)
```

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

``` r
VlnPlot(xe_residuals, features = "nCount_Xenium", group.by = "second_type", cols = pal, pt.size = 0)
```

    ## Rasterizing points since number of points exceeds 100,000.
    ## To disable this behavior set `raster=FALSE`

![](/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/SPLIT/doc/Reassign_residual_transcripts_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->
