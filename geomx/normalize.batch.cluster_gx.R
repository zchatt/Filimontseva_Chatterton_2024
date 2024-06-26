# Normalization of GeoMx using quantile normalization
# Please refer to "normalization_review_gx.R" for methods in selecting quantile normalization

# Libraries
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/giotto_env.R")
library(preprocessCore)
library(Giotto)
library(Seurat)
library(ggpubr)
library(harmony)
library(terra)

#########
## INPUT ##
#########

# slide.name <- "hu_brain_001"
# analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/analysis"

# run_name = "geomx_sep2023"
# analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis"
# rdata <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis/geomx_sep2023_gt0.gs0.01_qc.gx.Rdata"

run_name = "geomx_oct2023"
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
rdata <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_gt0.01.gs0.01_qc.gx.Rdata"

# run_name = "geomx_oct2023_min"
# analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min"
# rdata <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min/geomx_oct2023_min_gt0.gs0.01_qc.gx.Rdata"
# 

############################################################################################
#### Part 1 : Convert to Seurat object and quantile normalise
############################################################################################
setwd(analysis_dir)
# Load qc data generated from geomx_lowlevel.R, should be in "/analysis" of the name "slide.name_qc.gx.Rdata"
load(rdata)

# normalise expression data and store in "data" slot - we have changed order since "as.Seurat" no longer working
# gxdat_s@assays$RNA holds normalised data
raw <- exprs(gxdat)
norm_quant <- as.data.frame(normalize.quantiles(raw))
colnames(norm_quant) <- colnames(raw)
row.names(norm_quant) <- row.names(raw)

# convert to Seurat object 
#gxdat_s <- as.Seurat(gxdat,ident = "SampleID", normData = "exprs", forceRaw=TRUE) # no longer working
gxdat_s <- CreateSeuratObject(as.matrix(norm_quant), meta.data = sData(gxdat))
gxdat_s@assays$RNA$counts <- raw
gxdat_s@assays$RNA$data <- norm_quant

# add metadata
tmp <- gxdat_s@meta.data
table(row.names(tmp) == row.names(sData(gxdat)))
#tmp["ROI"] <- sData(gxdat)$ROI
tmp <- cbind(tmp,sData(gxdat))
gxdat_s@meta.data <- tmp


############################################################################################
#### Part 2 : Pre-batch clustering
############################################################################################
gxdat_s <- FindVariableFeatures(object = gxdat_s )
gxdat_s <- ScaleData(object = gxdat_s)
gxdat_s <- RunPCA(gxdat_s, verbose = FALSE)

## plot PC1/2
p1 <- DimPlot(object = gxdat_s, reduction = "pca", pt.size = .5, group.by = "slide name")
p2 <- VlnPlot(object = gxdat_s, features = "PC_1", pt.size = .1, group.by = "slide name")

ggsave("bclus_1_pca12.png",
       ggarrange(p1, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

ggsave("bclus_2_pca.exp.png",
       ggarrange(p2, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")


############################################################################################
#### Part 3 : Batch correction
############################################################################################
## Batch correct using Harmony
gxdat_s <- gxdat_s %>% 
  RunHarmony("slide name", plot_convergence = TRUE)

## plot PC1/2
p1 <- DimPlot(object = gxdat_s, reduction = "harmony", pt.size = 2, group.by = "slide name")
p2 <- VlnPlot(object = gxdat_s, features = "harmony_1", pt.size = .1, group.by = "slide name")

ggsave("bclus_3_pca12.harmony.png",
       ggarrange(p1, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

ggsave("bclus_4_pca.exp.harmony.png",
       ggarrange(p2, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")


###########################################################################################
#### Part 3 : Dimension reduction
############################################################################################

gxdat_s <- gxdat_s %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

## plot
p1 <- DimPlot(gxdat_s, reduction = "umap", group.by = "slide name", pt.size = 1, split.by = 'slide name')
p2 <- DimPlot(gxdat_s, reduction = "umap", label = TRUE)

ggsave("bclus_5_runname.umap.png",
       ggarrange(p1, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

ggsave("bclus_6_clust.umap.png",
       ggarrange(p2, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

# variables of interest to overlay
vars <- c("Diagnosis","Brainregion","Brainbank_ID","Sex","DetectionThreshold","segment")
px1 <- DimPlot(gxdat_s, reduction = "umap",  group.by=vars[1])
px2 <- DimPlot(gxdat_s, reduction = "umap",  group.by=vars[2])
px3 <- DimPlot(gxdat_s, reduction = "umap",  group.by=vars[3])
px4 <- DimPlot(gxdat_s, reduction = "umap",  group.by=vars[4])
px5 <- DimPlot(gxdat_s, reduction = "umap",  group.by=vars[5])
px6 <- DimPlot(gxdat_s, reduction = "umap",  group.by=vars[6])

ggsave("bclus_Dx.umap.png",
       ggarrange(px1,px2,px4,px5, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")
ggsave("bclus_Brainregion.umap.png",
       ggarrange(px2, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")
ggsave("bclus_Brainbank_ID.umap.png",
       ggarrange(px3, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")
ggsave("bclus_Sex.umap.png",
       ggarrange(px4, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")
ggsave("bclus_DetectionThreshold.umap.png",
       ggarrange(px5, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")
ggsave("bclus_x.umap.png",
       ggarrange(p2,px1,px2,px5, ncol=2,nrow=2) + bgcolor("white"),
       device = "png")
ggsave("bclus_x2.umap.png",
       ggarrange(p2,px6,px4,px4, ncol=2,nrow=2) + bgcolor("white"),
       device = "png")

# save data
save(gxdat_s, file = paste0(run_name,"_seurat.Rdata"))
