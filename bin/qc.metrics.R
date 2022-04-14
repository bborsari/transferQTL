.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

#************
# LIBRARIES *
#************

library(ggplot2)
library(pheatmap)
library(data.table)


#************
# FUNCTIONS *
#************

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}


function1 <- function(m, metric) {
  
  out <- data_summary(m, varname=metric, groupnames=c("Tissue"))
  out$metric <- metric
  colnames(out)[2] <- "value"
  return(out)
  
}


mcc <- function(TN, TP, FN, FP) {
  
  num = (TP*TN - FP*FN)
  denom = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  return((num/denom))
  
}


#********
# BEGIN *
#********

# 1. read palette
pal <- read.delim("/nfs/users2/rg/dgarrido/PhD/projects/GTEx/V8/sqtls/run_1/info/info.txt", 
                  h=T, sep="\t", stringsAsFactors = F)
pal$tissue <- gsub("_-_", "_", pal$tissue)
pal$tissue <- gsub("\\(", "", pal$tissue)
pal$tissue <- gsub("\\)", "", pal$tissue)

palette <- pal$tissue_color_hex
names(palette) <- pal$tissue


# 2. set wd
setwd("/users/rg/bborsari/eQTLs.model.nf/")


# 3. list of entex tissues
entex.tissues <- read.delim("EN-TEx.chromatin.tissues.txt", h=F, sep="\t", stringsAsFactors = F)
entex.tissues <- entex.tissues[, "V1"]
palette <- palette[names(palette) %in% entex.tissues]


# 4. read QC matrix for every tissue
m <- data.frame(stringsAsFactors = F)

for (tissue in entex.tissues) {
  
  tmp <- fread(paste0("run/", tissue, "/plots/quality.metrics.tsv"), header = F, data.table = F)
  colnames(tmp) <- c("Tissue",
                     "Accuracy",
                     "Kappa",
                     "AccuracyLower",
                     "AccuracyUpper",
                     "AccuracyNull",
                     "AccuracyPValue",
                     "McnemarPValue",
                     "Sensitivity",
                     "Specificity",
                     "Pos.Pred.Value",
                     "Neg.Pred.Value",
                     "Precision",
                     "Recall",
                     "F1",
                     "Prevalence",
                     "Detection.Rate",
                     "Detection.Prevalence",
                     "Balanced.Accuracy",
                     "TN", "FP", "FN", "TP",
                     "donor_tissue")
  
  print(nrow(tmp))
  
  #tmp <- tmp[, c("Tissue", "Balanced.Accuracy", "donor_tissue")]
  
  m <- rbind(m, tmp)
  
}

rm(tmp)

mcc.v <- c()
for ( i in 1:nrow(m) ){
  
  mcc.v <- c(mcc.v, mcc(TN = as.numeric(m[i, "TN"]),
                        TP = as.numeric(m[i, "TP"]),
                        FN = as.numeric(m[i, "FN"]),
                        FP = as.numeric(m[i, "FP"])))
  
}


m$MCC <- mcc.v


# 5. create vector w/ abbreviated labels
match.names <- c("ADRNLG", "AORTA", "ARTCRN", "ARTTBL", "BREAST", "CLNSGM",
                 "CLNTRN", "ESPGES", "ESPSQE", "ESPMSM", "HRTAA", "HRTLV",
                 "LIVER", "LUNG", "GASMED", "NERVET", "OVARY", "PNCREAS",
                 "PRSTTE", "SKINNS", "SKINS", "PEYERP", "SPLEEN", "STMACH",
                 "TESTIS", "THYROID", "UTERUS", "VAGINA")
names(match.names) <- sort(unique(m$Tissue))


# 6. heatmap of balanced.accuracy (donor: skin ) vs. all targets
m.h <- dcast(m, Tissue~donor_tissue, value.var = "Balanced.Accuracy")
rownames(m.h) <- match.names[m.h$Tissue]
m.h <- m.h[order(m.h$Skin_Sun_Exposed_Lower_leg, decreasing = T), , drop = F]
sorted.tissues <- m.h$Tissue
m.h$Tissue <- NULL
m.h <- m.h[!(rownames(m.h) %in% "SKINS"), c("Skin_Sun_Exposed_Lower_leg"), drop=F]
colnames(m.h) <- "SKINS"

ann_colors <- list(Tissue = palette,
                   eQTLs_Prevalence = c('#f7f7f7','#525252'))

# 7. row annotation = prevalence, tissue
m.anno <- dcast(m, Tissue~donor_tissue, value.var = "Prevalence")
rownames(m.anno) <- match.names[m.anno$Tissue]
m.anno <- m.anno[!(rownames(m.h) %in% "SKINS"), c("Skin_Sun_Exposed_Lower_leg", "Tissue"), drop=F]
colnames(m.anno)[1] <- "eQTLs_Prevalence"

# pdf("~/public_html/paper_ENTEx/model/heatmap.skins.pdf", height = 10, width = 6)
# pheatmap(m.h, 
#          cluster_rows = F, cluster_cols = F,
#          color = c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e'),
#          border_color = NA,
#          cellwidth = 10,
#          cellheight = 10,
#          annotation_row = m.anno,
#          annotation_colors = ann_colors,
#          annotation_names_row = F,
#          annotation_legend = T,
#          main = "Balanced Accuracy")
# dev.off()



# 8. for every target tissue, summarize balanced accuracy across different donor tissues
# aka main figure 8b

m.BA <- data_summary(m, varname="Balanced.Accuracy", 
                     groupnames=c("donor_tissue"))
m.BA <- m.BA[order(m.BA$Balanced.Accuracy, decreasing = T), ]
sorted.tissues <- sorted.tissues[sorted.tissues != "Skin_Sun_Exposed_Lower_leg"]
sorted.tissues <- c("Skin_Sun_Exposed_Lower_leg", sorted.tissues)
m.BA$donor_tissue <- factor(m.BA$donor_tissue, levels = sorted.tissues)

pdf("~/public_html/paper_ENTEx/model/barplot.BA.pdf", width = 3.5, height = 2.5)
ggplot(data = m.BA) + 
  geom_errorbar( aes(x=donor_tissue, 
                     ymin=Balanced.Accuracy-sd, 
                     ymax=Balanced.Accuracy+sd), 
                 colour="black",
                 width = 0.4) +
  geom_bar(aes(x=donor_tissue, 
               y=Balanced.Accuracy, 
               fill=donor_tissue), 
           stat="identity") +
  coord_cartesian(ylim=c(0,1)) +
  # scale_fill_continuous(low="#fcc5c0", high="#ae017e") +
  # scale_color_manual(values = rev(c("black", "white"))) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(panel.border = element_rect(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 13, hjust = .5)) +
  guides(fill=F, color=F) +
  # scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9)) +
  xlab("donor tissues (GTEx)") +
  ylab("Balanced Accuracy") +
  scale_x_discrete(labels = match.names[m.BA$donor_tissue])
dev.off()



