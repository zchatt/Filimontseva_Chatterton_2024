# regression analysis of NM variables and patient age

# Part 1: Extract NM quantification's
# Part 2: NM Area x Age
# Part 3: NM Density x Age
# Part 4: NM class x Age



library(ggpp)
library(readxl)
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
library(tidyr)

theme_set(theme_minimal())

## inputs
analysis_dir <- "/Users/zacc/USyd/NM_analysis"
setwd(analysis_dir)
quant_data = "/Users/zacc/USyd/NM_analysis/df_iNMeNM_230124.Rdata"
meta_data = "/Users/zacc/USyd/NM_analysis/NM_data_170124/ASAP1-IHC cohort.xlsx"

##############################################################
### Part 1: Extract NM quantification's  ###
##############################################################
# load NM quant data
load(quant_data)

# load metadata
meta <- read_xlsx(meta_data,1)
meta$Diagnosis[meta$Diagnosis == "Ct"] <- "CTR"

# # merge and format
df_agg_merged <- merge(df_agg,meta,by="Brainbank_ID",all.x=TRUE)
df_agg_merged$Diagnosis <- factor(df_agg_merged$Diagnosis,levels=c("CTR","ILBD", "ePD","lPD"))
df_agg_merged <- df_agg_merged[df_agg_merged$ROI != "ROI 1",]
df_agg_merged$log10_Area <- log10(df_agg_merged$Area..µm..)
df_agg_merged$Sex[df_agg_merged$Sex == "1"] <- "M"
df_agg_merged$Sex[df_agg_merged$Sex == "2"] <- "F"
df_agg_merged$Diagnosis_stage <- as.character(df_agg_merged$Diagnosis)
df_agg_merged$Diagnosis_stage[df_agg_merged$Diagnosis_stage == "CTR"] <- 0
df_agg_merged$Diagnosis_stage[df_agg_merged$Diagnosis_stage == "ILBD"] <- 1
df_agg_merged$Diagnosis_stage[df_agg_merged$Diagnosis_stage == "ePD"] <- 2
df_agg_merged$Diagnosis_stage[df_agg_merged$Diagnosis_stage == "lPD"] <- 3
df_agg_merged$Brainregion <- df_agg_merged$ROI
df_agg_merged$Brainregion[df_agg_merged$Brainregion %in% c("LC")] <- "A6"
df_agg_merged$Brainregion[df_agg_merged$Brainregion %in% c("SND","SNM","SNL","SNV")] <- "A9"
df_agg_merged$Brainregion[df_agg_merged$Brainregion %in% c("VTA")] <- "A10"
df_agg_merged$Direction.Red <- df_agg_merged$`Mean (Red)`/ (df_agg_merged$`Mean (Red)` + df_agg_merged$`Mean (Green)` + df_agg_merged$`Mean (Blue)`)*100
df_agg_merged$Direction.Green <- df_agg_merged$`Mean (Green)`/ (df_agg_merged$`Mean (Red)` + df_agg_merged$`Mean (Green)` + df_agg_merged$`Mean (Blue)`)*100
df_agg_merged$Direction.Blue <- df_agg_merged$`Mean (Blue)`/ (df_agg_merged$`Mean (Red)` + df_agg_merged$`Mean (Green)` + df_agg_merged$`Mean (Blue)`)*100
df_agg_merged <- df_agg_merged[complete.cases(df_agg_merged$log10_Area),]
df_agg_merged$Estimated.gray.value <- (0.299*df_agg_merged$`Mean (Red)`) + (0.587*df_agg_merged$`Mean (Green)`) + (0.114 * df_agg_merged$`Mean (Blue)`)
colnames(df_agg_merged) <- make.names(colnames(df_agg_merged))

# format age group
df_agg_merged$Age <- as.numeric(df_agg_merged$Age)
dplot <- df_agg_merged
dplot <- dplot %>%
  mutate(
    # Create categories
    age_group = dplyr::case_when(
      Age <= 60 ~ "<60",
      Age > 60 & Age <= 80  ~ "60-80",
      Age > 80 & Age <= 90 ~ "80-90",
      Age > 90             ~ "> 90"
    ),
    # Convert to factor
    age_group = factor(
      age_group,
      level = c("<60","60-80","80-90", "> 90")
    )
  )
df_agg_merged <- dplot


# colour palettes
Diagnosis_col = c("CTR"= "white","ILBD" = "#00AFBB", "ePD" = "#E7B800","lPD" = "red")
age_group_col = magma(4)
ROI_col = c("SNL" = "darkorchid1",
            "SNV" = "purple",
            "SND" = "purple3",
            "SNM" = "purple4",
            "VTA" = "forestgreen",
            "LC" = "yellow",
            "SNV_young" = "pink")

RGB_col = c("Red" = "red",
            "Blue" = "blue",
            "Green" = "green")




#####################################################
### Part 2: NM Area x Age
#####################################################
data_table <- df_agg_merged[df_agg_merged$Diagnosis == "CTR" ,]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]

# linear mixed-effects model of iNM size between diagnosis in each Brain region; area
y_list <- c("log10_Area")
res <- list()
roi_uniq <- c("VTA","SNV","LC")
main_lab <- c("VTA","SNV","LC")
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$ROI == roi_cont,]
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "age_group"
    y_variable = y_list[i]
    x_lab = "Age"
    y_lab = "Log10(Area [µm²])"
    colour_palette = Diagnosis_col
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Sex) + (1 | PMD)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(age_group = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(main_lab[z]) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(1.5,4)
    
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
                                      y.position = 3.8, step.increase = 0.1,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
    #ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM.Brainregion.pdf", arrange, width = 8, height = 6)
#ggsave("Vln_iNM.ROI.pdf", arrange)


# linear mixed-effects model of iNM size between diagnosis in each Brain region; n/um2
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57),]
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID_recode, ROI, Brainregion,Age,Sex,PMD,ROI.Area..µm.., sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$n_per.mm2 <- (dplot$n / as.numeric(dplot$ROI.Area..µm..)) * 1000000
data_table <- dplot[complete.cases(dplot$Brainregion),]
colnames(data_table) <- make.names(colnames(data_table))

y_list <- c("n_per.mm2")
res <- list()
roi_uniq <- c("A10","A9","A6")
main_lab <- c("VTA","SNpc","LC")

roi_uniq <- unique(data_table$ROI)
main_lab <- roi_uniq
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$Brainregion == roi_cont,]
  
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "Diagnosis"
    y_variable = y_list[i]
    x_lab = "Diagnosis"
    y_lab = "iNM (No./mm²)"
    colour_palette = Diagnosis_col
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table2)
    
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
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(main_lab[z]) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(0,200)
    
    
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
      # add p-values to plot and save the result
      stat.test <- stat.test[stat.test$p.adj < 0.05, ]
      bxp <- bxp + stat_pvalue_manual(stat.test,
                                      y.position = 250, step.increase = 0.1,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
    #ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM_n_per.mm2.Brainregion.pdf", arrange,width = 8, height = 6)
#ggsave("Vln_iNM_n_per.mm2.ROI.pdf", arrange)