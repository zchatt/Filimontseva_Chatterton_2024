# scripts to evaulate cell-type proportion differences in GeoMx following RCTD cell-type deconvolution for iNM manuscript

# libraries
library(edgeR)
library(dplyr)
library(purrr)
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
library(ComplexHeatmap)


##############################
### Input 
##############################

setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results_run180124/")
rctd_kamath = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results_run180124/meta_rctd_Kamath_n63.Rdata"
rctd_webber = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results_run180124/meta_rctd_Webber_n10.Rdata"

ROI_col = c("SNV" = "purple",
            "SNM" = viridis(6)[2],
            "SND" = viridis(6)[3],
            "SNL" = viridis(6)[4],
            "VTA" = viridis(6)[5],
            "LC" = viridis(6)[6],
            "SNV_young" = "dodgerblue",
            "RN" = "red4")

CA_col = c("CALB1_CALCR" = "purple",
                  "CALB1_CRYM_CCDC68" = "purple",
                  "CALB1_GEM" = "purple",   
                  "CALB1_PPP1R17" = "purple",
                  "CALB1_RBP4" = "purple",
                  "CALB1_TRHR" = "purple",      
                  "SOX6_AGTR1" = "dodgerblue",
                  "SOX6_DDT" = "dodgerblue",
                  "SOX6_GFRA2" = "dodgerblue",    
                  "SOX6_PART1" = "dodgerblue",
                  "NE" = "forestgreen")


##############################
### 1) Format data 
##############################
## combine data from Kamath decon SN/VTA/RN and Webber decon LC
load(rctd_webber)
m_lc <- meta_rctd
load(rctd_kamath)
m_k <- meta_rctd
meta_dat <- merge(m_lc,m_k,by="row.names",keep = TRUE)
row.names(meta_dat) <- meta_dat$Row.names
meta_dat <- meta_dat[,-grep("*\\.y$", colnames(meta_dat))]
colnames(meta_dat) <- gsub("*\\.x$","",colnames(meta_dat))

## format diagnosis
meta_dat$Diagnosis_2 <- meta_dat$Diagnosis
meta_dat$Diagnosis_2[meta_dat$Diagnosis_2 %in% c("ePD",'lPD')] <- "PD"

# select columns
df <- meta_dat
columns_to_calculate <- colnames(df)[67:139]

# Step 1: Group by segment, ROI, and Diagnosis, then calculate the mean for each group
group_means <- df %>%
  group_by(segment, ROI, Diagnosis_2) %>%
  summarise(across(all_of(columns_to_calculate), mean, na.rm = TRUE), .groups = 'drop')

# Step 2: Filter for CTR means
ctr_means <- group_means %>%
  filter(Diagnosis_2 == "CTR")

# Step 3: Join the CTR means back with the original dataset to have control means for comparison
df_with_ctr_means <- df %>%
  left_join(ctr_means, by = c("segment", "ROI"), suffix = c("", "_ctr_mean"))

# Step 4: Calculate the percentage of the control mean for each Diagnosis within each segment and ROI
df_with_pct_of_ctr <- df_with_ctr_means %>%
  mutate(across(all_of(columns_to_calculate),
                ~(.x / get(paste0(cur_column(), "_ctr_mean")) * 100),
                .names = "pct_of_ctr_{.col}"))

meta_dat_total <- df_with_pct_of_ctr

# test correct prc have been calculated 
mean(meta_dat_total$pct_of_ctr_NE[meta_dat_total$ROI == "SNV" & meta_dat_total$Diagnosis_2 == "CTR" & meta_dat_total$segment == "TH"])
mean(meta_dat_total$pct_of_ctr_SOX6_DDT[meta_dat_total$ROI == "SNV" & meta_dat_total$Diagnosis_2 == "CTR" & meta_dat_total$segment == "TH"])



##############################
### 2) Plot Heatmaps of percentage changes
##############################
Diagnosis_col = c("CTR"= "grey", "PD" = "red")

# get control data for ROIs
meta_dat <- meta_dat_total[meta_dat_total$ROI %in% c("RN","SND","SNL","SNM","SNV","VTA") & meta_dat_total$Diagnosis_2 %in% c("CTR","PD"),]
dim(meta_dat)

# change segment annota
meta_dat$segment[meta_dat$segment == "TH"] <- "TH+ mask"
meta_dat$segment[meta_dat$segment == "Full ROI"] <- "No mask"

# cells
y_list <- c("pct_of_ctr_CALB1_CALCR","pct_of_ctr_CALB1_CRYM_CCDC68","pct_of_ctr_CALB1_GEM", 
            "pct_of_ctr_CALB1_PPP1R17","pct_of_ctr_CALB1_RBP4","pct_of_ctr_CALB1_TRHR","pct_of_ctr_SOX6_AGTR1",
            "pct_of_ctr_SOX6_DDT","pct_of_ctr_SOX6_GFRA2", "pct_of_ctr_SOX6_PART1")

# y_list <- c("CALB1_CALCR","CALB1_CRYM_CCDC68","CALB1_GEM", 
#             "CALB1_PPP1R17","CALB1_RBP4","CALB1_TRHR","SOX6_AGTR1",
#             "SOX6_DDT","SOX6_GFRA2", "SOX6_PART1")

### Heatmap plots
# aggregate data
x = meta_dat[,y_list ]
y = meta_dat
ct_agg_mean <- as.data.frame(x) %>%
  group_by(y$ROI,y$segment,y$Diagnosis_2) %>% 
  summarise_all("mean")

# # plot data frames
dplot <- ct_agg_mean
dplot1 <- t(dplot[,4:ncol(dplot)])
dplot1 <- log10(dplot1) # convert to log10 for visualisation of the cell prop
row.names(dplot1) <- gsub("pct_of_ctr_","",row.names(dplot1))

# column annotations
anno_df = data.frame(
  Region = dplot$`y$ROI` ,
  Diagnosis = dplot$`y$Diagnosis_2`,
  Masking = dplot$`y$segment` 
)


ha = HeatmapAnnotation(df = anno_df,
                       col = list(Masking = c("No mask"= "#00AFBB","TH+ mask" = "pink"),
                                  Diagnosis = Diagnosis_col,
                                  Region = ROI_col)
)
# row colors
tmp <- row.names(dplot1)
tmp[grep("SOX6",tmp)] <- "dodgerblue" 
tmp[grep("CALB1",tmp)] <- "purple"
tmp[grep("NE",tmp)] <- "forestgreen"
row_annotation = c(col = tmp)


# draw heatmap
pdf("heatmap_rctd_prc_control_DA.pdf",width = 8, height = 4)
hm <- Heatmap(dplot1, cluster_columns = FALSE,
              col = mako(100),
              top_annotation = ha,
              column_split=sort(as.factor((dplot$`y$ROI`))),
              row_names_gp = gpar(col = row_annotation, fontsize = 7),
              heatmap_legend_param = list(title = "log10(% of CTR)")) 
draw(hm,
     column_title = "Kamath (n=63)",
     column_title_gp=grid::gpar(fontsize=16))
dev.off()

##############################
### 3) Contrast percentage changes between PD and Control
##############################
# select masking
segment <- "Full ROI"

#segment <-"Full ROI"
y_list <- c("pct_of_ctr_CALB1_CALCR","pct_of_ctr_CALB1_CRYM_CCDC68","pct_of_ctr_CALB1_GEM", 
            "pct_of_ctr_CALB1_PPP1R17","pct_of_ctr_CALB1_RBP4","pct_of_ctr_CALB1_TRHR","pct_of_ctr_SOX6_AGTR1",
            "pct_of_ctr_SOX6_DDT","pct_of_ctr_SOX6_GFRA2", "pct_of_ctr_SOX6_PART1")

roi_uniq <- c("SND","SNM","SNV","VTA","SNL")
main_lab = segment

# format numerical variables
meta_dat <- meta_dat_total[meta_dat_total$segment == segment & meta_dat_total$Diagnosis_2 %in% c("CTR","PD"),]
meta_dat$DV200 <- as.numeric(meta_dat$DV200)
meta_dat$Age <- as.numeric(meta_dat$Age)
meta_dat$GenesDetected <- as.numeric(meta_dat$GenesDetected)
colnames(meta_dat) <- make.names(colnames(meta_dat))
data_table <- meta_dat[,c(y_list,"ROI","Brainbank_ID","Age","Sex","PMD.hs","DV200","Diagnosis_2")]

y_labs <- y_list
res_stats <- list()
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$ROI == roi_cont,]
  
  for(i in 1:length(y_list)){
    tryCatch({
      set.seed(123)
      print(i)
      # define variables
      x_variable = "Diagnosis_2"
      y_variable = y_list[i]
      x_lab = "Diagnosis_2"
      y_lab = y_variable
      
      # Fit a linear mixed-effects model
      model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID) + (1 | Age) + (1 | Sex) + (1 | PMD.hs) + (1 | DV200)")), data = data_table2)
      print(model)
      
      # perform post-hoc tests to compare different regions using Tukey method
      posthoc <- glht(model, linfct = mcp(Diagnosis_2 = "Tukey"))
      summary(posthoc)
      
      # create stat.test df
      set.seed(123)
      summary_posthoc <- summary(posthoc)
      group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
      group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
      p.adj = as.numeric(summary_posthoc$test$pvalues)
      stat.test <- as.data.frame(cbind(group1,group2,p.adj))
      stat.test$Estimate <- summary_posthoc$test$coefficients
      stat.test$p.adj <- as.numeric(stat.test$p.adj)
      stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                       symbols = c("***", "**", "*", ".", " "))
      
      stat.test$ROI <- roi_uniq[z]
      stat.test$CellType <- y_list[i]
      stat.test$Reference <- main_lab
      stat.test$segment <- segment
      
      print(stat.test)
      
      # apply to list 
      res_stats[[length(res_stats) + 1]] <- stat.test
      
    }, error = function(e) {
      message(paste("Error in iteration with ROI", roi_cont, "and variable", y_list[i], ":", e$message))
      # Optionally, you can add code here to handle the error, e.g., logging it or assigning a default value
    })
  }
}

stat.summary <- do.call("rbind", res_stats)

# write summary to file
df <- as.data.frame(apply(stat.summary,2,as.character))
#write.table(df, file="Percentage_DA_ctr_TH.stat.summary.txt", sep="\t",row.names = F, quote = F)
write.table(df, file="Percentage_DA_ctr_Full.stat.summary.txt", sep="\t",row.names = F, quote = F)


### Plots of 



# prepare data for plotting
data_table_long <- meta_dat_total[ meta_dat_total$Diagnosis_2 %in% c("CTR","PD") & 
                                    meta_dat_total$ROI == "SNV", ] %>%
  pivot_longer( cols = y_list,  names_to = "Condition", values_to = "Value")
data_table_long$Condition <- as.factor(data_table_long$Condition)


# select data
dplot <- data_table_long[data_table_long$segment == "TH" & data_table_long$ROI == "SNV",]

# set theme
theme_set(
  theme_classic() +
    theme(legend.position = "top")
)

# select data
dplot <- data_table_long[data_table_long$segment == "TH" & data_table_long$ROI == "SNV",]


# summarise data mean and sd
df.summary <- dplot %>%
  group_by(ROI,segment,Diagnosis_2, Condition) %>%
  summarise(
    sd = sd(Value, na.rm = TRUE),
    len = mean(Value)
  )
df.summary


# order cells by mean
tmp <- aggregate(Value ~ Condition, dplot, mean)
fact_lvls <- tmp[order(-tmp[,"Value"]),][,"Condition"]
dplot$Condition <- factor(dplot$Condition, levels = fact_lvls)

# Line plots with jittered points
dodge_width <- 0.9
ggplot(dplot, aes(x = Condition, color = Diagnosis_2)) +
  geom_point(aes(y = Value),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = dodge_width)) +
  geom_errorbar(aes(ymin = len - sd, ymax = len + sd, group = interaction(Condition, Diagnosis_2)),
                data = df.summary2, width = 0.25,
                position = position_dodge(dodge_width), color = "black") +
  scale_color_manual(values = c("grey", "red")) +
  theme(legend.position = "top") + geom_hline(yintercept = 100, lty = 3) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10)) + ylab("Cell Proportion (% Control)")


# save
arrange <- ggarrange(plotlist=list(g1,NULL,g2), nrow=2, ncol=2, widths = c(2,2))
ggsave("Boxplots_iNM_colors.pdf", arrange,width = 8, height = 6)

