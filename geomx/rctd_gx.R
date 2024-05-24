# scripts for performing cell-deconvolutions in Nanostring GeoMx data using RCTD method described in https://doi.org/10.1038/s41587-021-00830-w

library(spacexr)
library(Seurat)
library(Matrix)
library(doParallel)
library(ggplot2)
library(plyr)
library(data.table)
library(msigdbr)
library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(EnhancedVolcano)
library(viridis)
library(rstatix)
library(ComplexHeatmap)
library(readxl)
register(SerialParam())

############################################################################################
#### Inputs
############################################################################################

analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
rdata = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata"
results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results_run180124'
contrast_path <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/contrasts_matrix.xlsx"
snrnaseq_ref <- "/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/merged_kam.sil.web_seurat.Rdata"
setwd(analysis_dir)

############################################################################################
###### Part 2:  Format geomx spatial data
############################################################################################
# setwd
setwd(results_folder)

# load geomx normalised data
load(rdata)
# select samples
#keep_index <- gxdat_s$Diagnosis == "CTR"
keep_index <- gxdat_s$Diagnosis != "NTC"
# load in counts matrix
#counts <- gxdat_s@assays$GeoMx@counts[,keep_index]  # load in counts matrix
counts <- gxdat_s@assays$RNA@counts[,keep_index]  # load in counts matrix
# create metadata matrix 
meta <- gxdat_s@meta.data[keep_index,]
# confirm names
table(colnames(counts) == row.names(meta))
# load in coordinates
coords = as.data.frame(gxdat_s@reductions$umap@cell.embeddings[keep_index,]) # we are using UMAP coordinates as dummy variables until we obtain DSP ROI locations
colnames(coords) <- c("imagerow","imagecol")
# process counts
nUMI <- colSums(counts) # In this case, total counts per spot

############################################################################################
###### Part 3: RCTD to obtain cell proportions
############################################################################################
# load snrnaseq reference
load(snrnaseq_ref)

### i) deconvolute using Kamath n= 63 referene
set.seed(123)
# construct ST to deconvolute for SN tissue
group2 <- meta$ROI %in% c("SND","SNL","SNM","SNV","VTA","LC","RN")
puck <- SpatialRNA(coords[group2,], counts[,group2], nUMI[group2])
barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')

# construct references
# NOTE; taking all cells from Kamath
toMatch <- c(names(table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Kamath"])))
group1 <- merge.combined.sct@meta.data$cell_type_merge %in% toMatch
sc_obj <- merge.combined.sct[ ,group1]

# remove cell-types with < 25 cells AND/OR cells of interes
cells_keep <- names(table(sc_obj@meta.data$cell_type_merge)[table(sc_obj@meta.data$cell_type_merge) > 25])
sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% cells_keep]
# get count data and restrict to only genes within GeoMx
counts_sc <- sc_obj[["RNA"]]$counts
counts_sc <- counts_sc[row.names(counts_sc) %in% row.names(puck@counts),]

# set cell types and quant nUMI
cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
cell_types <- as.factor(cell_types)
nUMIsc <- colSums(counts_sc)
## create reference
reference <- Reference(counts_sc, cell_types, nUMIsc)

##  run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# add deconvolution to meta-data and save
table(row.names(normalize_weights(myRCTD@results$weights)) == row.names(meta))
meta_rctd <- cbind(meta,normalize_weights(myRCTD@results$weights))
#save(meta_rctd,file = "meta_rctd_Kamath_n63.Rdata")



#### ii) deconvolute using Webber LC reference
# construct ST to deconvolute for SN tissue
group2 <- meta$ROI %in% c("SND","SNL","SNM","SNV","VTA","LC","RN")
puck <- SpatialRNA(coords[group2,], counts[,group2], nUMI[group2])
barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')

# construct references
toMatch <- c(names(table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Webber"])))
group1 <- merge.combined.sct@meta.data$cell_type_merge %in% toMatch
sc_obj <- merge.combined.sct[ ,group1]

# remove cell-types with < 25 cells
cells_keep <- names(table(sc_obj@meta.data$cell_type_merge)[table(sc_obj@meta.data$cell_type_merge) > 25])
sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% cells_keep]
# get count data and restrict to only genes within GeoMx
counts_sc <- round(10^sc_obj[["RNA"]]$counts,0)
counts_sc <- counts_sc[row.names(counts_sc) %in% row.names(puck@counts),]

# set cell types and quant nUMI
cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
cell_types <- as.factor(cell_types)
nUMIsc <- colSums(counts_sc)
## create reference
reference <- Reference(counts_sc, cell_types, nUMIsc)

##  run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# add deconvolution to meta-data and save
table(row.names(normalize_weights(myRCTD@results$weights)) == row.names(meta))
meta_rctd <- cbind(meta,normalize_weights(myRCTD@results$weights))
#save(meta_rctd,file = "meta_rctd_Webber_n10.Rdata")


