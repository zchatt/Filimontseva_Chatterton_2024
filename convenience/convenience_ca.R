# convenience scripts for CA degeneration manuscript

# libs
library(viridis)
library(readxl)
library(RColorBrewer)

#---------- colour palettes ----------
Diagnosis_col = c("CTR"= "grey",
                  "ILBD" = "#00AFBB", 
                  "ePD" = "#E7B800",
                  "lPD" = "red")

age_group_col = magma(4)

ROI_col = c("SNL" = "darkorchid1",
            "SNV" = "purple",
            "SND" = "purple3",
            "SNM" = "purple4",
            "VTA" = "forestgreen",
            "LC" = "yellow",
            "SNV_young" = "pink")

Cell_col = c("SOX6_AGTR1" = brewer.pal(n = 9, name = "YlGnBu")[7] ,
             "SOX6_DDT" = brewer.pal(n = 9, name = "YlGnBu")[6] ,
             "SOX6_PART1" = brewer.pal(n = 9, name = "YlGnBu")[5],
             "SOX6_GFRA2" = brewer.pal(n = 9, name = "YlGnBu")[4],
             "CALB1_CALCR" = brewer.pal(n = 9, name = "YlOrRd")[2],
             "CALB1_TRHR" = brewer.pal(n = 9, name = "YlOrRd")[3],
             "CALB1_PPP1R17" = brewer.pal(n = 9, name = "YlOrRd")[4],
             "CALB1_CRYM_CCDC68" = brewer.pal(n = 9, name = "YlOrRd")[5],
             "CALB1_GEM" = brewer.pal(n = 9, name = "YlOrRd")[6] ,
             "CALB1_RBP4" = brewer.pal(n = 9, name = "YlOrRd")[7],
             "NE" = "yellow")

Cell_col2 = c("SOX6_AGTR1_TYRpos" = "red",
              "SOX6_AGTR1_TYRneg" = "purple",
              "NE" = "yellow")


#---------- Genes of Interest related to NM ----------
skin_melanin_enzymes = c("TYR", "TYRP1","DCT") # note DCT = "TYRP2"
alt_PD = c("ALDH1A1", "ALDH3A1", "DDT", "CRABP1", "MAGED2","RBP4")
CA_precursor_NM_stock_trans_metab_A9.A10 = c("TH","DDC","COMT","ALDH","MAO","AR","ADH")
CA_precursor_NM_stock_trans_metab_A6 = c("TH", "DDC", "DBH", "MAO", "COMT", "MHPG")
CA_functional_DA = c("SLC18A2", "SLC6A3","SLC18A1")
CA_functional_NE = c("SLC18A2", "SLC18A1", "SLC6A2")
Stress_granules.free_radical_scavenging = toupper(c("DDX6","DDX1", "DDX3a", "DDX17","DDX3X", "PABPC1", "eIF3A", "eIF3B", "eIF4G1",
                                                    "G3BP1","G3BP2", "RAB33A", "MAP1L3CB2","MAP1LC3B","NDUFS7", "ACOT8", "COA7", "PRDX2",
                                                    "PRDX1","PRDX2","PRDX5","TUBA1B","GPN1","GPNMB","EGFR","BLVTB","BLVRB","GSTT1", "SIRT5"))
Lysosome_pathways = c("CTSD", "LAMP2", "MAP1LC3B", "SCARB2", "UBA52", "GBA1","GBA","FBXO16", "ATG5", "HSAPA4L", "HSAPA4","HSPA9", "UCHL1")
upstream_euNM = c("TH","DDC","DBH","DDT") #enzymes upstream of euNM, DAQ-DAC-DHI (enzyme involved in these two steps and products, eg DDT)
genes_link_skinpig.PD = c("GCH1", "GPNMB", "HERC2", "LRRK2","MC1R","PRKN","SNCA", "TPCN2","TYR", "TRPM7", "VPS35") # Genes linking skin pigmentation and Parkinson’s disease

# yuhong list genes interest
yf_genes <- c("DDT",
              "MC1R","TYR","TYRP1",
              "OCA2","ALDH1A1","ALDH1A3","CRABP1","MAGED2","RBP4",
              "TH","DDC","DBH","MAO","MAOA","COMT","MHPG","VMAT2","DAT","SLC6A3","VMAT1","DDX6","DDX1","DDX3a",
              "DDX17","PABPC1","elF3A","elF3B","elF4G1","G3BP1","G3BP2","RAB33A","MAP1L3CB2","NDUFS7",
              "ACOT8","COA7","PRDX2","PRDX1","PRDX5","TUBA1B","GPN1","GPNMB","EGFR","BLVTB","GSTT1",
              "SIRT5","CTSD","LAMP2","MAP1LC3B","SCARB2","UBA52","GBA1","FOXO16","FOXO1","ATG5","HSPA4L",
              "HSAPA4","HSPA4","HSPA9","UCHL1")

# pigmentation network
pigmentation_network <- c("ADAM17","ADAMTS20","BRCA1","ED1","EDA","EDN3","EDNRB","EGFR",
                          "FGFR2","IKBKG","KIT","KITLG","KRT2A","KRT2","LMX1A","MCOLN3","MITF","PAX3","SFXN1","SNAI2",
                          "SOX10","SOX18","WNT1","WNT3A","DCT","GPNMB","MATP","SLC45A2","SLC45A2","RAB38","SILV","PMEL","TYR","TYRP1",
                          "AP3B1","AP3D1","VPS33A","CNO","BLOC1S4","HPS1","HPS3","HPS4","HPS5","HPS6","LYST","MU","GSTM1",
                          "OA1","GPR143","P","PLDN","BLOC1S6","RABGGTA","MLPH","MYO5A","MYO7A","RAB27A","ASIP","ATRN","GGT1","GL","MC1R",
                          "MGRN1","POMC","ATP7A","ATP7B","BCL2","ERCC2","DCX","GSR","ITG2B","ITGA2B",
                          "ITGB3","MAP3K14","PH","PPY","C4A","C4B","U2AF1","ZFP362","ZNF362","MECP2","PITX2","RS1","SCO2","TYMP","MLC","MLC1")

# yuhong + pigmentation genes
yf_pigment_genes <- unlist(read.delim(file="/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min/pigmentation+YuHong.genes_fixed.txt", header = F)[,1])
yf_pigment_genes <-  toupper(yf_pigment_genes)
yf_pigment_genes <- unique(yf_pigment_genes )

# other undefines
other <- c("ABCB6","ANKRD27","GLS","LAMP1","RAB32","RAB9A","MYEF2","SGSM2","FOXO1")

# create list 
gene_interest_list <- list(yf_pigment_genes, skin_melanin_enzymes,alt_PD,CA_precursor_NM_stock_trans_metab_A9.A10,CA_precursor_NM_stock_trans_metab_A6,
                           CA_functional_DA,CA_functional_NE,Stress_granules.free_radical_scavenging,Lysosome_pathways,upstream_euNM,
                           genes_link_skinpig.PD,pigmentation_network,other)

names(gene_interest_list) <- c("Genes of interest literature Review by Y.F. & A.F.","Skin Melanin Enzyme","Altered in PD","CA precursor / NM stock trans metab (A9.A10)","CA precursor / NM stock trans metab (A6)",
                               "CA functional (DA)","CA functional (NE)","Stress granules/ Free radical scavenging","Lysosome pathways","Upstream NM",
                               "Genes linking skin pigmentation and PD","Pigmentation network","Other genes of interest")

# tabulate
gene_df <- do.call(rbind, lapply(names(gene_interest_list), function(list_name) {
  data.frame(
    Gene = unlist(gene_interest_list[[list_name]]),
    Source = list_name
  )
}))

combined_df <- aggregate(Source ~ Gene, data = gene_df, FUN = function(x) paste(unique(x), collapse = ", "))


# write to table
write.table(combined_df,file = "NM_genes_interest.txt", sep="\t", quote = F)


# PD genetics gene lists
g4pd_cnv_path <- "/Users/zacc/USyd/spatial_transcriptomics/data/ref_genesets/t_cnv20200820.txt"
g4pd_cv_path <- "/Users/zacc/USyd/spatial_transcriptomics/data/ref_genesets/t_common_variant20200820.txt"
g4pd_rg_path <- "/Users/zacc/USyd/spatial_transcriptomics/data/ref_genesets/t_rare_gene20200820.txt"
g4pd_rv_path <- "/Users/zacc/USyd/spatial_transcriptomics/data/ref_genesets/t_rare_variant20200820.txt"
g4pd_deg_path <- "/Users/zacc/USyd/spatial_transcriptomics/data/ref_genesets/t_gene_expression20200820.txt"

tmp <- read.delim(g4pd_rg_path, sep='\t', row.names=NULL, header=T)
g4pd_rare_gene <- gsub(" ","",unique(unlist(tmp$Gene)))
tmp <- read.delim(g4pd_cnv_path, sep='\t', row.names=NULL, header=T)
g4pd_cnv <- unlist(strsplit(gsub(" ","",unique(unlist(tmp$Gene))),";|,|/"))
tmp <- read.delim(g4pd_cv_path, sep='\t', row.names=NULL, header=T)
g4pd_common_variant <- unlist(strsplit(gsub(" ","",unique(unlist(tmp$Gene))),";|,|/"))
tmp <- read.delim(g4pd_rv_path, sep='\t', row.names=NULL, header=T)
g4pd_rare_variant <- unlist(strsplit(gsub(" ","",unique(unlist(tmp$Gene))),";|,|/"))
tmp <- read.delim(g4pd_deg_path)
colnames(tmp) <- paste0("Published_PD_DEG_evidence_",colnames(tmp))
g4pd_DEG <- gsub(" ","",unique(unlist(tmp$Published_PD_DEG_evidence_Gene)))

gene_interest_list <- list(alt_PD,pigmentation_network,
                           g4pd_rare_gene, g4pd_cnv,g4pd_common_variant,g4pd_rare_variant,g4pd_DEG)

names(gene_interest_list) <- c("alt_PD","pigmentation_network",
                               "g4pd_rare_gene", "g4pd_cnv","g4pd_common_variant","g4pd_rare_variant","g4pd_DEG")
