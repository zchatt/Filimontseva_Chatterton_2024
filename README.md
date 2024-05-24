# Filimontseva_Chatterton_2024
Scripts used in the data analysis desribed in Filimontseva & Chatterton et al, "Changes in the production and type of neuromelanin are early features in neurons that selectively degenerate in Parkinson’s disease", 2024.


## geomx
1) Low level processing and quality control of Nanostring GeoMx data [quality_control_gx.R](quality_control_gx.R)
2) Alternative low level processing and quality control of Nanostring GeoMx data without applying LOQ thresholds [quality_control_min_gx.R](quality_control_min_gx.R)
3) Review of normalization methods Nanostring GeoMx data [normalization_review_gx.R](normalization_review_gx.R)
4) Quantile Normalization, batch correction and clustering of Nanostring GeoMx data [normalize.batch.cluster_gx.R](normalize.batch.cluster_gx.R)
5) Deconvolution of cell-type proportions in Nanostring GeoMx data using RCTD method [rctd_gx.R](rctd_gx.R)
6) Evaluation of cell-type proportions by brain region in GeoMx datat [diff_celltype.gx.R](diff_celltype.gx.R)
7) Analysis of Nanostring GeoMx data to identify DEGs associated with NM granule classes [deg_inm_gx.R](deg_inm_gx.R)
8) Evaluation of gene expression of genes of interest within GeoMx data [geneexp_gx.R](geneexp_gx.R)

## nm
Scripts used for the downstream analysis of Neuromelanin, following intraneuronal Neuromelanin quentification described in https://github.com/zchatt/Neuromelanin

1) Analysis of intraneuronal NM granules using linear regression and k-means clustering to define granule classes [regression_cluster_inm.R](nm/regression_cluster_inm.R)


## snranseq
1) Scripts evaluating published single-cell RNA sequencing datasets for gene expression related to Neuromelanin biosynthesis [snrnaseq_nm.R](snrnaseq_nm.R)

