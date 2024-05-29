# scripts to evaulate cell-type proportion differences following RCTD cell-type deconvolution

# libraries
library(edgeR)
library(dplyr)
library(plyr)
library(EnhancedVolcano)
library(ggplot2)
library(ggridges)
library(MASS)
library(ggplot2)
library(viridis)
library(ggpubr)
library(EnvStats)
library(rstatix)
library(ggridges)
library(lme4)
library(multcomp)
library(randomForest)
library(e1071)
library(caret)
library(nnet)
library(forcats)
library(ggcorrplot)
library(plotly)
library(tidyverse)
library(ComplexHeatmap)

#----------- inputs ----------- 
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results_run180124/")

# Kamath RCTD results (n=63)
rctd_meta_data = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results_run180124/meta_rctd_Kamath_n63.Rdata"
load(rctd_meta_data)

ROI_col = c("SNV" = "purple",
            "SNM" = viridis(6)[2],
            "SND" = viridis(6)[3],
            "SNL" = viridis(6)[4],
            "VTA" = viridis(6)[5])

Kamath_DA_col = c("CALB1_CALCR" = "purple",
                  "CALB1_CRYM_CCDC68" = "purple",
                  "CALB1_GEM" = "purple",   
                  "CALB1_PPP1R17" = "purple",
                  "CALB1_RBP4" = "purple",
                  "CALB1_TRHR" = "purple",      
                  "SOX6_AGTR1" = "dodgerblue",
                  "SOX6_DDT" = "dodgerblue",
                  "SOX6_GFRA2" = "dodgerblue",    
                  "SOX6_PART1" = "dodgerblue")


##############################
### 2) regression analysis - Kamath et al.
##############################
# i) Evaluate each cell-type by region; Kamath n=63
# get control data for ROIs
meta_dat <- meta_rctd[meta_rctd$Diagnosis %in% "CTR" & meta_rctd$segment == "TH" & meta_rctd$ROI %in% c("RN","SND","SNM","SNV","VTA","SNL"),]
dim(meta_dat)

# format numerical variables
meta_dat$DV200 <- as.numeric(meta_dat$DV200)
meta_dat$Age <- as.numeric(meta_dat$Age)
meta_dat$GenesDetected <- as.numeric(meta_dat$GenesDetected)
colnames(meta_dat) <- make.names(colnames(meta_dat))

# select cols of interest
y_list <- c("CALB1_CALCR","CALB1_CRYM_CCDC68","CALB1_GEM", "CALB1_PPP1R17","CALB1_RBP4","CALB1_TRHR",
            "SOX6_AGTR1","SOX6_DDT","SOX6_GFRA2", "SOX6_PART1")
data_table <- meta_dat[,c(y_list,"ROI","Brainbank_ID","Age","Sex","PMD.hs","DV200")]

y_labs <- y_list
res <- list()
for(i in 1:length(y_list)){
  set.seed(123)
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = y_variable
  main_lab = "Kamath, All cells (n=63)"
  colour_palette = ROI_col 
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID) + (1 | Age) + (1 | Sex) + (1 | PMD.hs) + (1 | DV200)")), data = data_table)
  print(model)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
  summary(posthoc)
  
  # format for plotting
  set.seed(123)
  tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  print(tmp)
  fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # make violin plot 
  bxp <- ggboxplot(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), plot.title = element_text(size = 10)) + 
    ylab(y_lab) + xlab("") + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black") + ggtitle(main_lab) +
    ylim(0,(max(data_table$y)+(max(data_table$y) * 0.9)))
  
  
  # create stat.test df
  set.seed(123)
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
  if(sum(stat.test$p.adj < 0.05) > 1){
    # add p-values to plot and save the result
    stat.test <- stat.test[stat.test$p.adj < 0.05, ]
    bxp <- bxp + stat_pvalue_manual(stat.test,
                                    y.position = (max(data_table$y)+(max(data_table$y) * 0.2)), step.increase = 0.15,
                                    label = "p.adj.signif")
  }
  # assign plot
  res[[i]] <- bxp
}

# save
arrange <- ggarrange(plotlist=res[1:6], nrow=2, ncol=3, widths = c(2,2))
ggsave("Bxp_CellType.CTR.TH_n1.6_kamath_n63.pdf", arrange, width = 8, height = 6)

arrange <- ggarrange(plotlist=res[7:12], nrow=2, ncol=3, widths = c(2,2))
ggsave("Bxp_CellType.CTR.TH_n7.12_kamath_n63.pdf", arrange, width = 8, height = 6)
