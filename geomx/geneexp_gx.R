## Analysis of genes or gene-sets of interest from Spatial Transcriptomics (ST) 

library(dplyr)
library(viridis)
library(data.table)
library(ggpubr)
library(rstatix)
library(spacexr)
library(Seurat)
library(edgeR)
library(ComplexHeatmap)
library(tidyr)


###############
### Input
###############
# analysis dir
analysis_dir = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
setwd(analysis_dir)

# colour palettes
Diagnosis_col = c("CTR"= "grey","ILBD" = "#00AFBB", "ePD" = "#E7B800","lPD" = "red")
age_group_col = magma(4)
ROI_col = c("SNL" = "darkorchid1",
            "SNV" = "purple",
            "SND" = "purple3",
            "SNM" = "purple4",
            "VTA" = "forestgreen",
            "LC" = "yellow",
            "SNV_young" = "pink")


############################################################################################
### # 1a. Evaluate expression of published DA neuronal scRNA markers of PD in GeoMx data as % of controls
############################################################################################
# define  genes interest
genes_to_plot <- c("TH","SOX6","AGTR1","ALDH1A1","KCNJ6","ANXA1","RIT2","SATB1","CALB1","SLC17A6") # markers of vulnerable populations and resistance_markers "CALB1","SLC17A6"

# ## evaluate genes of interest in non-thresheld data
# load("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min/geomx_oct2023_min_gt0.gs0.01_qc.gx.Rdata")
# goi <- genes_to_plot
# goi_df <- data.frame(
#   Gene = goi,
#   Number_of_Segments = fData(gxdat)[goi, "DetectedSegments"],
#   DetectionRate = percent(fData(gxdat)[goi, "DetectionRate"]))
# ss <- tableGrob(goi_df)
# 
# ggsave("table.relevantgenes.png",ss,device = "png")
# # unexpectedly SOX6 is below gene detection threshold


# ST - GeoMx thresheld data
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata")
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat <- as.data.frame(gxdat_s@meta.data)

# define df and meta data 
genes_to_plot <- genes_to_plot[genes_to_plot %in% row.names(exp_dat)]

meta_dat <- meta_dat[meta_dat$ROI %in% c("SNM") &
                       meta_dat$segment == "TH" &
                       meta_dat$Diagnosis %in% c("CTR","ePD","lPD"),]
meta_dat$condt <- meta_dat$Diagnosis != "CTR"
exp_dat <- exp_dat[genes_to_plot,row.names(meta_dat)]

# make midbrain ROI
meta_dat$ROI2 <- "midbrain"


# DEG analysis between between PD and Control subjects within each ROI
# # format covariates
meta_dat$Sex <- as.factor(meta_dat$Sex)
meta_dat$condt <- as.factor(meta_dat$condt)
meta_dat$Age  <- as.numeric(meta_dat$Age)
meta_dat$DV200  <- as.numeric(meta_dat$DV200)
meta_dat$Brainbank_ID <- as.factor(meta_dat$Brainbank_ID)

# limma voom
unique_var <- unique(meta_dat$ROI2)
res <- list()
for (val in 1:length(unique_var)){  
  print(unique_var[val])
  
  # subset for var of interest
  meta_dat2 <- meta_dat[meta_dat$ROI2 == unique_var[val],]
  exp_dat2 <- exp_dat[,row.names(meta_dat2)]
  
  dge <- DGEList(exp_dat2, group = meta_dat2$condt)
  dge <- calcNormFactors(dge)
  design <- model.matrix( ~ condt + DV200 + Age + Sex, data=meta_dat2)
  vm <- voom(dge, design = design, plot = TRUE)
  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)
  tt <- topTable(fit, n = Inf, adjust.method = "BH")
  
  tt$contrast <- unique_var[val]
  res[[val]] <- tt
}

# name results
names(res) <- unique_cell_type
res <- lapply(res, function(x) {
  x$Gene <- row.names(x)
  row.names(x) <- NULL
  return(x)
})

# write summary to file
stat.summary <- do.call(rbind,res)
df <- as.data.frame(apply(stat.summary,2,as.character))
#write.table(df, file="deg_res.vun.genes_ctr.pd_TH.stat.summary.txt", sep="\t",row.names = F, quote = F)
#write.table(df, file="deg_res.vun.genes_ctr.pd_Full.stat.summary.txt", sep="\t",row.names = F, quote = F)

# plot valcano
dplot <- df
dplot$adj.P.Val_log10 <- -log10(as.numeric(dplot$adj.P.Val))
dplot$logFC <- as.numeric(dplot$condtTRUE)
dplot$adj.P.Val <- as.numeric(dplot$adj.P.Val)

sig_genes <- dplot[dplot$adj.P.Val < 0.05,]
tmp <- dplot[dplot$adj.P.Val < 0.05 & dplot$logFC > 0.5,]
up_genes_expected <- tmp[tmp$Gene %in% c("CALB1","SLC17A6"),]
up_genes_unexpected <- tmp[!tmp$Gene %in% c("CALB1","SLC17A6"),]

tmp <- dplot[dplot$adj.P.Val < 0.05 & dplot$logFC < -0.5,]
down_genes_expected <- tmp[tmp$Gene %in% c("TH","SOX6","AGTR1","ALDH1A1","KCNJ6","ANXA1","RIT2","SATB1"), ] 
down_genes_unexpected <- tmp[!tmp$Gene %in% c("TH","SOX6","AGTR1","ALDH1A1","KCNJ6","ANXA1","RIT2","SATB1"), ] 


g1 <- ggplot(dplot, aes(x = logFC, y = adj.P.Val_log10)) +
  geom_point(aes(colour = contrast), 
             alpha = 0.8, 
             shape = 16,
             size = 4) + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", size=0.2) + 
  geom_vline(xintercept = c(-0.5,0.5),
             linetype = "dashed",size=0.2) +
  geom_vline(xintercept = c(0),size=0.2) +
  scale_colour_manual(values = ROI_col) + 
  scale_x_continuous(breaks = c(seq(-5, 5, 1)),     
                     limits = c(-5, 5)) +
  geom_label_repel(data = up_genes_unexpected, # Add labels last to appear as the top layer  
                   aes(label = Gene), color = "grey",
                   force = 2,size = 2,
                   nudge_x = 2,nudge_y = 5, max.overlaps= 100) +
  geom_label_repel(data = down_genes_unexpected, # Add labels last to appear as the top layer  
                   aes(label = Gene), color = "black",
                   force = 2,size = 2,
                   nudge_x = -1, max.overlaps= 100) +
  geom_label_repel(data = up_genes_expected, # Add labels last to appear as the top layer  
                   aes(label = Gene), color = "dodgerblue",
                   force = 2, size = 2,
                   nudge_y = 10,nudge_x = 2, max.overlaps= 100) +
  geom_label_repel(data = down_genes_expected, # Add labels last to appear as the top layer  
                   aes(label = Gene), color = "red",
                   force = 2, size = 2,
                   nudge_y = 1,nudge_x = -1, max.overlaps= 100) +
  labs(title = "TH+ mask",
       x = "logFC",
       y = "-log10(adjusted P-value)",
       colour = "Brain Region") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + ylim(0,15)
  
# save plot
arrange <- ggarrange(plotlist=list(g1), nrow=2, ncol=2)
ggsave("volcano_res.vun.genes_ctr.pd_TH.pdf", arrange,width = 8, height = 6)
#ggsave("volcano_res.vun.genes_ctr.pd_Full.pdf", arrange,width = 8, height = 6)


#### plot % of controls

# define groupings and collection list
grp_list <- list(c("SND", "SNM", "SNV","VTA"),
                 c("SNV"),c("SNM"),c("VTA"))
br_list <- list(c("midbrain"),
                c("SNV"),c("SNM"),c("VTA"))
col_list <- list()

# collect mean/ sd for each group
for (i in 1:length(grp_list)){
  # ST - GeoMx thresheld data
  load("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata")
  exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
  meta_dat <- as.data.frame(gxdat_s@meta.data)
  
  # define df and meta data 
  genes_to_plot <- genes_to_plot[genes_to_plot %in% row.names(exp_dat)]
  
  meta_dat <- meta_dat[meta_dat$ROI %in% grp_list[[i]] &
                         meta_dat$segment == "TH" &
                         meta_dat$Diagnosis %in% c("CTR","ePD","lPD"),]
  meta_dat$condt <- meta_dat$Diagnosis != "CTR"
  exp_dat <- exp_dat[genes_to_plot,row.names(meta_dat)]
  
  
  # combine data
  data_table <- cbind(meta_dat,t(exp_dat))
  
  # select columns
  columns_to_calculate <- genes_to_plot
  
  # group by ROI and Diagnosis, then calculate the mean for each group
  df <- data_table
  group_means <- df %>%
    group_by(condt) %>%
    summarise(across(all_of(columns_to_calculate), mean, na.rm = TRUE), .groups = 'drop')
  
  # filter for CTR means
  ctr_means <- group_means %>%
    filter(condt == FALSE)
  
  # join the CTR means back with the original dataset to have control means for comparison
  ctr_means$orig.ident <- "SeuratProject" # create common variable to left join
  
  df_with_ctr_means <- data_table %>%
    left_join(ctr_means, by = "orig.ident", suffix = c("", "_ctr_mean"))
  
  # calculate the percentage of the control mean for each Diagnosis within each segment and ROI
  df_with_pct_of_ctr <- df_with_ctr_means %>%
    mutate(across(all_of(columns_to_calculate),
                  ~(.x / get(paste0(cur_column(), "_ctr_mean")) * 100),
                  .names = "pct_of_ctr_{.col}"))
  
  # test correct prc have been calculated 
  mean(df_with_pct_of_ctr$pct_of_ctr_TH[df_with_pct_of_ctr$condt == FALSE])
  
  # # define df and meta data 
  # genes_to_plot <- genes_to_plot[genes_to_plot %in% row.names(exp_dat)]
  # meta_dat$condt <- meta_dat$Diagnosis != "CTR"
  # exp_dat <- exp_dat[genes_to_plot,row.names(meta_dat)]
  
  # assign to data_table
  data_table <- df_with_pct_of_ctr
  
  # convert gene expression to long form
  data_long <- gather(data_table, condition, measurement, pct_of_ctr_TH:pct_of_ctr_SLC17A6, factor_key=TRUE)
  data_long$condition <- gsub("pct_of_ctr_","",data_long$condition)
  
  ## plot in single violin plot
  tmp <- aggregate(measurement ~ condition, data_long[data_long$condt == TRUE,], mean)
  fact_lvls <- tmp[order(-tmp[,"measurement"]),][,"condition"]
  data_long[,"condition"] <- factor(data_long[,"condition"], levels = fact_lvls)
  data_long$condt <- as.factor(data_long$condt)
  
  # ##  remove outliers for plotting
  # x <- data_long$measurement
  # outliers <- !x %in% boxplot.stats(x)$out
  # dplot <- data_long[outliers,]
  
  # g1 <- ggboxplot(
  #   data_long, x = "condition", y = "measurement",
  #   color = "condt",
  #   palette = c("grey","red"),
  #   add = "mean_sd", size = 0.2) +
  #   theme(legend.position = "none") +
  #   geom_hline(yintercept=100, linetype="dashed", color = "grey") +
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10),
  #         plot.title = element_text(size = 10)) +
  #   ylab("Norm. Count (% of CTR's)") + xlab("") +
  #   ggtitle("") 
  
  
  # group data to calculate mean and sd
  grouped <- group_by(data_long, condition, condt)
  df2 <- summarise(grouped, mean=mean(measurement), sd=sd(measurement))
  
  # assign name
  df2$BrainRegion <- br_list[[i]]

  # assign to collection list
  col_list[[i]] <- df2
  
  
}

# assign names
df3 <- do.call("rbind", col_list)

# select only PD
df3 <- df3[df3$condt == "TRUE",]

# plot barplots
p<- ggplot(df3, aes(x=condition, y=mean, fill=BrainRegion)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  labs( x="", y = "Norm. Count (% of CTR's)") +
  theme_classic() + 
  scale_fill_manual(values=c("grey","purple","purple4","forestgreen")) + 
  #ylim(0,330) +
  geom_hline(yintercept = 100, lty = 2)
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), legend.position = "none")

arrange <- ggarrange(plotlist=list(p), nrow=2, ncol=1, widths = c(2,2))
ggsave("bar_res.vun.cellmarkers_prc.ctr_CTRPD_midbrainregions_gx.pdf", arrange,width = 8, height = 6)
#ggsave("bar_res.vun.cellmarkers_prc.ctr_CTRPD_gx.pdf", arrange,width = 8, height = 6)
#ggsave("Vln_res.vun.cellmarkers_prc.ctr_CTRPD.pdf", arrange,width = 8, height = 6)
#ggsave("Vln_iNM_n_per.mm2.ROI.pdf", arrange)




############################################################################################
### # 1b. Evaluate expression of genes of interest between ROIs of Control GeoMx data
###########################################################################################

# ST - GeoMx thresheld data
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata")
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat <- as.data.frame(gxdat_s@meta.data)

# define df and meta data 
genes_to_plot <- genes_to_plot[genes_to_plot %in% row.names(exp_dat)]

meta_dat <- meta_dat[meta_dat$Brainregion %in% c("A9","A6") &
                       meta_dat$segment != "TH" &
                       meta_dat$Diagnosis %in% c("CTR"),]

# define condition to contrast
meta_dat$condt <- meta_dat$Brainregion == "A9"
exp_dat <- exp_dat[,row.names(meta_dat)]

# # make midbrain ROI
# meta_dat$ROI2 <- "midbrain"


# DEG analysis between between PD and Control subjects within each ROI
# # format covariates
meta_dat$Sex <- as.factor(meta_dat$Sex)
meta_dat$condt <- as.factor(meta_dat$condt)
meta_dat$Age  <- as.numeric(meta_dat$Age)
meta_dat$DV200  <- as.numeric(meta_dat$DV200)
meta_dat$Brainbank_ID <- as.factor(meta_dat$Brainbank_ID)

# limma voom
meta_dat2 <- meta_dat
exp_dat2 <- exp_dat

dge <- DGEList(exp_dat2, group = meta_dat2$condt)
dge <- calcNormFactors(dge)
design <- model.matrix( ~ condt + DV200 + Age + Sex, data=meta_dat2)
vm <- voom(dge, design = design, plot = TRUE)
fit <- lmFit(vm, design = design)
fit <- eBayes(fit)
tt <- topTable(fit, n = Inf, adjust.method = "BH")

## plot barplots
# genes_to_plot <- c("SLC31A1","ATP7A","ATP7B","ATOX1","MT1H","MMP14","MT2A","MT3","MT4")
# genes_to_plot <- c("MAOA","SLC18A2","ALDH1A1","TH")
# genes_to_plot <- c("FABP7","EIF4G1","HSPB1","HSPH1")
# genes_to_plot <- c("ABCA5","ABCA6","ABCA9","ABCA10")
# genes_to_plot <- c("ABCA5","HSP90AA1", "HSPA4L","HSPA4")
#genes_to_plot <- c("HSPA1A", "HSPA4L","HSPD1","HSPA4","HSPA9","HSP90AA1")

genes_to_plot <- genes_to_plot[genes_to_plot %in% row.names(exp_dat)]
exp_dat <- exp_dat[genes_to_plot,row.names(meta_dat)]
tt[genes_to_plot,]

dplot <- cbind(meta_dat,t(exp_dat))
data_long <- gather(dplot, condition, measurement, row.names(exp_dat)[1]:row.names(exp_dat)[length(row.names(exp_dat))], factor_key=TRUE)


# Calculate mean and standard deviation for each group
df1 <- data_long[,c("condt","condition","measurement")]
df_summary <- df1 %>%
  group_by(condt, condition) %>%
  summarise(
    mean_measurement = mean(measurement),
    sd_measurement = sd(measurement)
  )


df_summary$Brainregion <- as.character(df_summary$condt)
df_summary$Brainregion[df_summary$Brainregion == "TRUE"] <- "VTA"
df_summary$Brainregion[df_summary$Brainregion == "FALSE"] <- "SNC"
df_summary$Brainregion <- as.factor(df_summary$Brainregion)

# Create the bar plot with error bars extending upwards only
p1a <- ggplot(df_summary, aes(x=Brainregion, y=mean_measurement, fill=Brainregion)) +
  geom_bar(stat='identity', position=position_dodge(width=0.9), width=0.8, colour="black") + 
  geom_linerange(aes(ymin=mean_measurement, ymax=mean_measurement + sd_measurement), 
                 position=position_dodge(width=0.9)) +
  facet_wrap(~ condition, scales='free_y', ncol=4, strip.position='top') +
  scale_fill_manual(values=c("purple","forestgreen")) +
  ylab("Norm. Counts") + 
  theme_classic() + ylim(0,50) +
  theme(
    axis.text.x = element_text(angle=45, vjust=0.5, hjust=1, size=10),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = 'outside',
    strip.background = element_blank()
  )

arrange <- ggarrange(plotlist=list(p1a), nrow=2, ncol=1, widths = c(2,2))
#ggsave("barplot_a6a9_ctr_genesinterest_200524_9.pdf", arrange,width = 8, height = 6)
  



############################################################################################
### # 2. Evaluate expression of 130 pigmentation/ melanin genes of interest
############################################################################################

## 3) Plot heatmap of 130 genes
library(ComplexHeatmap)

# cells
y_list <- unique(meta_dat$cell_type_merge)

# aggregate data
x = t(plot_mat)
y = meta_dat
dplot <- as.data.frame(x) %>%
  group_by(y$cell_type_merge) %>% 
  summarise_all("mean")

# # plot data frames
dplot1 <- t(dplot[,2:ncol(dplot)])
dplot1 <- scale(dplot1)
#dplot1 <- log10(dplot1)

# column annotations
anno_df = data.frame(
  Cell = dplot$`y$cell_type_merge`
)

ha = HeatmapAnnotation(df = anno_df,
                       col = list(Cell = Cell_col))

# row labels
tmp <- row.names(dplot1)
tmp[tmp %in%  other] <- "Genes interest"
tmp[tmp %in%  pigmentation_network] <- "Pigmentation Network"
tmp[tmp %in%  upstream_euNM] <- "Upstream euNM"
tmp[tmp %in%  c(CA_precursor_NM_stock_trans_metab_A9.A10,CA_precursor_NM_stock_trans_metab_A6)] <- "CA precursor"
tmp[tmp %in%  c(CA_functional_DA,CA_functional_NE)] <- "CA Functional"
tmp[tmp %in%  Stress_granules.free_radical_scavenging] <- "Stress Granule/ Scavenge"
tmp[tmp %in%  Lysosome_pathways] <- "Lysosome pathways"
tmp[tmp %in%  skin_melanin_enzymes] <- "Melanin enzymes"
tmp[tmp %in% yf_genes] <- "Genes interest"


# draw heatmap
pdf("heatmap_kamath.webber_n130.pdf",width = 10, height = 10)
hm <- Heatmap(dplot1, cluster_columns = TRUE,
              top_annotation = ha,
              row_split = tmp,
              heatmap_legend_param = list(title = "Log10(mean(Cell Prop.))"),
              column_km = 4, row_names_gp = gpar(fontsize = 4),
              width = ncol(dplot1)*unit(5, "mm"), 
              height = nrow(dplot1)*unit(1.2, "mm")
)
draw(hm)
dev.off()




############################################################################################
### # 3. Evaluate the expression of TYR in non-threshold GeoMx data
############################################################################################
## 1) load and format data
# load normalised GeoMx data with no threshold
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min/geomx_oct2023_min_seurat.Rdata")

# get the expression data 
dat <- gxdat_s@assays$RNA$data
meta <- gxdat_s@meta.data

# Summarize the total counts for each gene
gene_summaries <- dat %>%
  as.data.frame() %>%  # Convert the matrix to a data frame
  rownames_to_column(var = "Gene") %>%  # Move gene names to a column
  mutate(total_count = rowSums(.[-1]))  # Calculate total counts for each gene

gene_summaries <- gene_summaries[,c(1,904)]

# Name of the gene you want to highlight
gene_of_interest <- "TYR"

# Extract the total count for the gene of interest
tyr_count <- gene_summaries %>%
  filter(Gene == gene_of_interest) %>%
  pull(total_count)

## 2) Create a histogram using ggplot2
g1 <- ggplot(gene_summaries[1:18675,], aes(x = total_count )) +
  geom_density(aes(y=..density..)) +
  geom_segment(aes(x = tyr_count, xend = tyr_count, y = 0, yend = 7), linetype = "dashed", color = "grey20") +
  labs(title = "Normalized Counts Ctr TH+ neurons",
       x = "Normalized Counts (nt)",
       y = "Frequency") + theme_bw() + scale_x_log10() +
  geom_label(label=gene_of_interest, x=log10(tyr_count),
             y=7, label.padding = unit(0.55, "lines"),  label.size = 0.35,color = "grey20",
             fill="white") + ylim(0,8)


## 3) Violin plots between ROI for CTR 
gene_name <- "TYR"
gene <- dat[gene_name,]
dplot <- cbind(meta,gene = unlist(gene))
data_table <- dplot[dplot$Diagnosis == "CTR" & dplot$segment == "TH",]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[,c("ROI","gene","Brainbank_ID", "Age","Sex","PMD.hs","Plate_ID","DV200","DeduplicatedReads","GenesDetected")]

# linear mixed-effects model 
y_list <- c("gene")
# define variables
x_variable = "ROI"
y_variable = y_list[1]
x_lab = ""
y_lab = "TYR Normalized Counts (nt)"
colour_palette = ROI_col

# Fit a linear mixed-effects model
model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID) + (1 | Age) + (1 | Sex) + (1 | PMD.hs) + (1 | DV200) + (1 | Plate_ID)")), data = data_table)

# perform post-hoc tests to compare different regions using Tukey method
posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
summary(posthoc) # no sig diff between ROIs

# format for plotting
tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)

# Perform F-test
roi_cont <- unique(data_table$ROI)
for (val in 1:length(roi_cont)){
  f_test_result <- var.test(data_table[data_table$ROI == roi_cont[val],"gene"], data_table[data_table$ROI != roi_cont[val],"gene"])
  print(roi_cont[val])
  print(f_test_result)
  
}
# sig diff between SNV and LC vs else

# make violin plot 
bxp <- ggviolin(
  data_table, x = x_variable, y = y_variable, 
  fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab(y_lab) + xlab(x_lab) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                   geom="pointrange", color="black") + ylim(0,70)

# save plot
arrange <- ggarrange(plotlist=list(bxp), nrow=2, ncol=2)
ggsave(paste("violin_TYR.CTR.ROI.pdf"), arrange,  width = 8, height = 6)


## 4) Test for depth and gene detection differences between TYR high and low
high_threshold <- quantile(data_table$gene, 0.75)
low_threshold <- quantile(data_table$gene, 0.25)
# Classify as high, medium, or low
data_table$tyr_expression_class <- ifelse(data_table$gene > high_threshold, "High", 
                                          ifelse(data_table$gene < low_threshold, "Low", "A"))

t.test(data_table$DeduplicatedReads[data_table$tyr_expression_class == "High"],
       data_table$DeduplicatedReads[data_table$tyr_expression_class == "Low"])

t.test(data_table$GenesDetected[data_table$tyr_expression_class == "High"],
       data_table$GenesDetected[data_table$tyr_expression_class == "Low"])

# ii) Violin plots between ROIs x Dx 
# compare diagnosis within each ROI
gene_name <- "CCL3"
gene <- dat[gene_name,]
dplot <- cbind(meta,gene = unlist(gene))
data_table <- dplot[dplot$segment != "TH" ,]
colnames(data_table) <- make.names(colnames(data_table))
data_table$Diagnosis <- factor(data_table$Diagnosis,levels=c("CTR","ILBD", "ePD","lPD"))
data_table <- data_table[,c("ROI","gene","Brainbank_ID", "Age","Sex","PMD.hs","Plate_ID","DV200","DeduplicatedReads","Diagnosis","GenesDetected")]

y_list <- c("gene")
res <- list()
roi_uniq <- unique(data_table$ROI)
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$ROI == roi_cont,]
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "Diagnosis"
    y_variable = y_list[i]
    x_lab = "Diagnosis"
    y_lab = gene_name
    colour_palette = Diagnosis_col
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID) + (1 | Age) + (1 | Sex) + (1 | PMD.hs) + (1 | DV200) + (1 | Plate_ID)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(roi_cont) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(0,(max(data_table2[, y_variable] * 1.4)))
    
    # set x and y in data_table2 for plotting
    data_table2$x <- data_table2[, x_variable]
    data_table2$y <- data_table2[, y_variable]
    
    # create stat.test df
    summary_posthoc <- summary(posthoc)
    group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
    group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
    p.adj = as.numeric(summary_posthoc$test$pvalues)
    stat.test <- as.data.frame(cbind(group1,group2,p.adj))
    stat.test$p.adj <- as.numeric(stat.test$p.adj)
    stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", " "))
    print(stat.test)
    # get sig values
    if(any(stat.test$p.adj < 0.05)){
      print("test")
      # add p-values to plot and save the result
      stat.test <- stat.test[stat.test$p.adj < 0.05, ]
      bxp <- bxp + stat_pvalue_manual(stat.test,
                                      y.position = max(data_table2[, y_variable]) * 1.2, step.increase = 0.1,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
    # ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

# save plot
arrange <- ggarrange(plotlist=res, nrow=2, ncol=3)
ggsave(paste("violin_CCL3.Full.Diagnosis.ROI.pdf"), arrange)


# Perform F-test
roi_cont <- unique(data_table$ROI)
for (val in 1:length(roi_cont)){
  roi_cont <- roi_uniq[val]
  data_table2 <- data_table[data_table$ROI == roi_cont,]
  f_test_result <- var.test(data_table2[data_table2$Diagnosis == "CTR","gene"],
                            data_table2[data_table2$Diagnosis != "CTR","gene"])
  print(roi_cont)
  print(f_test_result)
  
}




y_list <- c("gene")
res <- list()
roi_uniq <- unique(data_table$Diagnosis)
res <- list()
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "Diagnosis"
  y_variable = y_list[i]
  x_lab = "Diagnosis"
  y_lab = y_labs[i]
  colour_palette = Diagnosis_col 
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID) + (1 | Age) + (1 | Sex) + (1 | PMD.hs) + (1 | DV200) + (1 | Plate_ID)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
  summary(posthoc)
  
  # format for plotting
  tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  
  # make violin plot 
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab(y_lab) + xlab(x_lab) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                     geom="pointrange", color="black") + ylim(0,700)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  
  # get sig values
  if(length(stat.test$p.adj < 0.05) > 1){
    # add p-values to plot and save the result
    stat.test <- stat.test[stat.test$p.adj < 0.05, ]
    bxp <- bxp + stat_pvalue_manual(stat.test,
                                    y.position = 250, step.increase = 0.16,
                                    label = "p.adj.signif") 
  }
  
  # save plot
  #ggsave(paste("violin_iNM.CTR",y_variable,".png"), bxp)
  
  # assign plot
  res[[i]] <- bxp
}


ggsave(paste("Vln_TYR.TH.Dx",y_variable,".",x_lab,".",roi_cont,".png"), arrange)



### Classification of high and low expressing ROIs ###
gene_name <- "TYR"
gene <- unlist(dat[gene_name,])
meta <- cbind(meta,gene = gene)

# Define the high and low thresholds based on quantiles
high_threshold <- quantile(gene, 0.75)
low_threshold <- quantile(gene, 0.25)

# Classify as high, medium, or low
meta$tyr_expression_class <- ifelse(gene > high_threshold, "High", 
                                    ifelse(gene < low_threshold, "Low", "A"))

t.test(meta$GeneDetectionRate[meta$tyr_expression_class == "High"],
       meta$GeneDetectionRate[meta$tyr_expression_class == "Low"])

table(meta$tyr_expression_class,meta$ROI)

# Prop test 
# Counts of high gene expressions and total observations in each group
high_count_vta <- 100 # replace with your actual count for VTA
total_vta <- 500 # replace with your total count for VTA
high_count_snv <- 120 # replace with your actual count for SNV
total_snv <- 500 # replace with your total count for SNV

# Two-proportion z-test
result <- prop.test(x = c(high_count_vta, high_count_snv), 
                    n = c(total_vta, total_snv),
                    correct = FALSE) # 'correct = FALSE' disables Yates' continuity correction

# View the result
print(result)

### DEG analysis between TYR high and low
# extract meta data
meta_dat <- meta[meta$segment == "TH" & meta$Diagnosis == "CTR",]
exp_dat  <- gxdat_s@assays$RNA$counts[,row.names(meta_dat)]
table(row.names(meta_dat) == colnames(exp_dat))

# format numerical variables
meta_dat$DV200 <- as.numeric(meta_dat$DV200)
meta_dat$Age <- as.numeric(meta_dat$Age)
meta_dat$n_per.µm2.iNM <- as.numeric(meta_dat$n_per.µm2.iNM)
meta_dat$n_per.µm2.eNM <- as.numeric(meta_dat$n_per.µm2.eNM)
meta_dat$Diagnosis_stage <- as.character(meta_dat$Diagnosis)
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "CTR"] <- 0
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "ILBD"] <- 1
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "ePD"] <- 2
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "lPD"] <- 3
meta_dat$Diagnosis_stage <- as.numeric(meta_dat$Diagnosis_stage)

# format factors
factor_format <- c("Brainbank_ID","Sex","Diagnosis","Brainregion","ROI")
for (i in 1:length(factor_format)){
  meta_dat[,factor_format[i]] <- as.factor(meta_dat[,factor_format[i]])
}
meta_dat$Diagnosis <- factor(meta_dat$Diagnosis, levels=c('ePD','CTR','lPD','ILBD'))

# create design matrix
targ <- meta_dat
options(na.action='na.omit')
design <- model.matrix(reformulate("Brainbank_ID + DV200 + Age + ROI + Sex + tyr_expression_class + GeneDetectionRate"),
                       data=targ, 
                       drop = FALSE)

# make names
colnames(design) <- make.names(colnames(design))

# create expression df & targ with design row.names
targ <- targ[row.names(design),]
y <- exp_dat[,row.names(targ)]

# fit design
v <- voom(y,design)
vfit <- lmFit(v)

# Perform LIMMA contrasts
cont.matrix <- makeContrasts(A="tyr_expression_classHigh - tyr_expression_classLow",levels=design)
cont.matrix <- makeContrasts(A="gene",levels=design)


fit2 <- contrasts.fit(vfit, cont.matrix)
vfit2 <- eBayes(fit2)
options(digits=3)

# Select significant DEGs and assign to list
tmp <- topTable(vfit2,number = Inf,sort.by="P")
res <- tmp[tmp$P.Value < 0.05,]


