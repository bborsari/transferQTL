.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

#************
# LIBRARIES *
#************

library(ggplot2)
library(data.table)
library(ggpubr)


#************
# FUNCTIONS *
#************

#+++++++++++++++++++++++++
# Function to calculate the mean (or median) and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables

data_summary_median <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum <- plyr::ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("median" = varname))
  return(data_sum)
}


data_summary_mean <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum <- plyr::ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#********
# BEGIN *
#********

# 1. set wd
setwd("/users/rg/bborsari/eQTLs.model.nf/")

# 2. list of entex tissues
entex.tissues <- read.delim("EN-TEx.chromatin.tissues.txt", h=F, sep="\t", stringsAsFactors = F)
entex.tissues <- entex.tissues[, "V1"]

print(length(entex.tissues))

# 3. for every pair donor_tissue/target_tissue
# compute the mean / median nominal_pval for TN, TP, FN, FP sets

df_mean <- data.frame(stringsAsFactors = F)
df_median <- data.frame(stringsAsFactors = F)

for (tissue in entex.tissues) {
  
  tmp <- fread(paste0("/nfs/users/rg/bborsari/eQTLs.model.nf/run/", tissue, "/plots/nominal.pvals.tsv"), header = F, data.table = F)
  colnames(tmp) <- c("variant_id", "pval_nominal", "slope", "slope_se", "set", "donor_tissue", "target_tissue")
  
  print(nrow(tmp))  

  tmp2_median <- data_summary_median(tmp, varname="pval_nominal", 
                       groupnames=c("set", "target_tissue"))
  tmp2_median$donor_tissue <- tissue
  df_median <- rbind(df_median, tmp2_median)
  
  tmp2_mean <- data_summary_mean(tmp, varname="pval_nominal", 
                                   groupnames=c("set", "target_tissue"))
  tmp2_mean$donor_tissue <- tissue
  df_mean <- rbind(df_mean, tmp2_mean)
  
  
}


# 4. make plots

df_mean$set <- factor(df_mean$set, levels = c("TP", "FN", "FP", "TN"))
df_median$set <- factor(df_median$set, levels = c("TP", "FN", "FP", "TN"))

my_comparisons <- list( c("TP", "FN"), c("FP", "TN"))

# # 4.1. mean nominal pval
# pdf("~/public_html/paper_ENTEx/model/nominal.pvals.mean.pdf",
#     width = 3, height = 3)
# ggplot(df_mean, aes(x=set, y=-log10(pval_nominal))) +
#   geom_boxplot() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 15),
#         panel.border = element_rect(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black")) +
#   ylab(parse(text = "-log[10]~(nominal~p~-value)")) +
#   stat_compare_means(comparisons = my_comparisons,
#                      label.y = c(6.5, 3),
#                      size = 5) +
#   ylim(1, 7.5)
# dev.off()

# 4.2. median nominal pval
pdf("~/public_html/paper_ENTEx/model/nominal.pvals.median.pdf",
    width = 3, height = 3)
ggplot(df_median, aes(x=set, y=pval_nominal)) +
  #geom_boxplot() +
  geom_violin() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        panel.border = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylab(parse(text = "-log[10]~(GTEx~p~-value)")) +
  xlab("SNV sets") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(32, 6),
                     size = 4.5) +
  ylim(0, 34)
dev.off()

pdf("~/public_html/paper_ENTEx/model/nominal.pvals.median.FP.TN.pdf",
    width = 3, height = 3)
ggplot(df_median[df_median$set %in% c("FP", "TN"), ], aes(x=set, y=-log10(pval_nominal))) +
  #geom_boxplot() +
  geom_violin() +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 15),
        panel.border = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
  # ylab(parse(text = "-log[10]~(GTEx~p~-value)")) +
  # xlab("GTEx SNV sets")
  # stat_compare_means(comparisons = my_comparisons,
  #                    label.y = c(32, 6),
  #                    size = 4.5) +
  # ylim(0, 34)
dev.off()
