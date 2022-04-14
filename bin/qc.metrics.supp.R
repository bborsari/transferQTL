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
  data_sum <- rename(data_sum, c("mean" = varname))
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



mcc.v <- c()
for ( i in 1:nrow(m) ){
  
  mcc.v <- c(mcc.v, mcc(TN = as.numeric(m[i, "TN"]),
                        TP = as.numeric(m[i, "TP"]),
                        FN = as.numeric(m[i, "FN"]),
                        FP = as.numeric(m[i, "FP"])))
  
}


m$MCC <- mcc.v
m$Ratio <- m$Accuracy / m$AccuracyNull


# 5. create vector w/ abbreviated labels
match.names <- c("ADRNLG", "AORTA", "ARTCRN", "ARTTBL", "BREAST", "CLNSGM",
                 "CLNTRN", "ESPGES", "ESPSQE", "ESPMSM", "HRTAA", "HRTLV",
                 "LIVER", "LUNG", "GASMED", "NERVET", "OVARY", "PNCREAS",
                 "PRSTTE", "SKINNS", "SKINS", "PEYERP", "SPLEEN", "STMACH",
                 "TESTIS", "THYROID", "UTERUS", "VAGINA")
names(match.names) <- sort(unique(m$Tissue))


# 6. summary dfs for other metrics
m.BA <- function1(m = m, metric = "Balanced.Accuracy")
m.SN <- function1(m = m, metric = "Sensitivity")
m.SP <- function1(m = m, metric = "Specificity")
m.PR <- function1(m = m, metric = "Precision")
m.RT <- function1(m = m, metric = "Ratio")
m.mcc <- function1(m = m, metric = "MCC")
m.PV <- function1(m = m, metric = "Prevalence")

m.merged <- rbind(m.BA, m.SN, m.SP, m.PR, m.mcc, m.PV, m.RT)
m.merged$metric <- gsub("Balanced.Accuracy", "Balanced\nAccuracy", m.merged$metric)
m.merged$metric <- gsub("Ratio", "Accuracy\nRatio", m.merged$metric)
m.merged$metric <- factor(m.merged$metric, levels = c("Balanced\nAccuracy", 
                                                      "Sensitivity", 
                                                      "Specificity", 
                                                      "Precision", 
                                                      "MCC", 
                                                      "Prevalence", 
                                                      "Accuracy\nRatio"))
m.merged$Tissue <- factor(m.merged$Tissue, levels = m.BA[order(m.BA$value), "Tissue"])


dummy <- data.frame(value = c(1, 1.5, rep(c(0.5, 1), 6)),
                    metric = rep(c("Accuracy\nRatio", 
                               "Balanced\nAccuracy", 
                               "Sensitivity", 
                               "Specificity", 
                               "Precision", 
                               "MCC", 
                               "Prevalence"), each = 2), 
                    stringsAsFactors=FALSE)

lop <- list()


lop[[1]] <- ggplot(data = m.merged[!(m.merged$metric %in% c("Accuracy\nRatio", "Prevalence")), ]) +
  geom_errorbar( aes(x=Tissue, ymin=value-sd, ymax=value+sd), 
                 colour="black",
                 width = 0.4) +
  facet_wrap(~metric, nrow = 1) +
  geom_point(aes(x=Tissue, y=value, color = Tissue)) +
  # geom_bar(aes(x=Tissue, y=value, fill = Tissue), stat="identity") +
  coord_flip(ylim = c(0.5, 1)) +
  scale_y_continuous(breaks = c(0.60, 0.80, 1.00), labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = palette) +
  # scale_fill_manual(values = palette) +
  theme_bw() +
  theme(panel.border = element_rect(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 11),
        strip.background.x = element_blank()) +
  guides(color=F, fill=F) +
  xlab("target tissues (GTEx)") +
  scale_x_discrete(labels = match.names[m.merged$Tissue]) 


lop[[2]] <- ggplot(data = m.merged[(m.merged$metric %in% c("Accuracy\nRatio")), ]) +
  geom_errorbar( aes(x=Tissue, ymin=value-sd, ymax=value+sd), 
                 colour="black",
                 width = 0.4) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  geom_point(aes(x=Tissue, y=value, color = Tissue)) +
  # geom_bar(aes(x=Tissue, y=value, fill = Tissue), stat="identity") +
  coord_flip(ylim=c(0.75, 1.75)) +
  scale_y_continuous() +
  scale_color_manual(values = palette) +
  # scale_fill_manual(values = palette) +
  theme_bw() +
  theme(panel.border = element_rect(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.background.x = element_blank()) +
  guides(color=F, fill=F) +
  xlab("target tissues (GTEx)") +
  scale_x_discrete(labels = match.names[m.merged$Tissue]) 


lop[[3]] <- ggplot(data = m.merged[(m.merged$metric %in% c("Prevalence")), ]) +
  geom_errorbar( aes(x=Tissue, ymin=value-sd, ymax=value+sd), 
                 colour="black",
                 width = 0.4) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  geom_point(aes(x=Tissue, y=value, color = Tissue)) +
  # geom_bar(aes(x=Tissue, y=value, fill = Tissue), stat="identity") +
  coord_flip(ylim=c(0, 1)) +
  scale_y_continuous() +
  scale_color_manual(values = palette) +
  # scale_fill_manual(values = palette) +
  theme_bw() +
  theme(panel.border = element_rect(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.background.x = element_blank()) +
  guides(color=F, fill=F) +
  xlab("target tissues (GTEx)") +
  scale_x_discrete(labels = match.names[m.merged$Tissue]) 

pdf("~/public_html/paper_ENTEx/model/barplot.metrics.supp.pdf", 
    width = 10, height = 6, useDingbats = F)
plot_grid(plotlist = lop, nrow = 1, rel_widths = c(0.75, 0.15, 0.15))
dev.off()


