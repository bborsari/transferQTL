.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(ggplot2)



#************
# FUNCTIONS *
#************


function1 <- function(m, file, pval = F) {
  
  tmp <- read.table(file, h=T, sep="\t", stringsAsFactors = F)
  
  if (pval) {
    
    tmp <- tmp[tmp$AccuracyPValue < 0.05, ]
    
  }

  m <- rbind(m, tmp)
  
  return(m)
  
}


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-plyr::ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}


function2 <- function(m, metric) {
  
  out <- data_summary(m, varname=metric, groupnames=c("target_tissue", "type"))
  out$metric <- metric
  colnames(out)[2:3] <- c("type", "value")
  return(out)
  
}


mcc <- function(TN, TP, FN, FP) {
  
  num = (TP*TN - FP*FN)
  denom = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  return((num/denom))
  
}


function3 <- function(m) {
  
  mcc.v <- c()
  for ( i in 1:nrow(m) ){
    
    mcc.v <- c(mcc.v, mcc(TN = as.numeric(m[i, "TN"]),
                          TP = as.numeric(m[i, "TP"]),
                          FN = as.numeric(m[i, "FN"]),
                          FP = as.numeric(m[i, "FP"])))
    
  }
  
  m$MCC <- mcc.v
  return(m)
  
}




#*********
#  BEGIN *
#*********

setwd("/nfs/users2/rg/bborsari/eQTLs.model.nf")

# 1. define color palette
palette <- c("SKINS" = "#7777FF", "SKINNS" = "#0000FF", "GASMED" = "#AAAAFF", "PNCREAS" = "#995522")

#-------------------------------------------------
# 2. make plot w/o filtering for pval Accuracy
#-------------------------------------------------

# 2.1. initialize dfs
df1 <- data.frame(stringsAsFactors = F) # ext. vs. GTEx
df2 <- data.frame(stringsAsFactors = F) # ext. vs. prediction
df3 <- data.frame(stringsAsFactors = F) # ext. vs. prediction (only consistent snps)

# 2.2. read QC dfs from tissues
# 2.2.1. muscle
df1 <- function1(m = df1, file = "validation/muscle/quality.metrics.1.tsv")
df2 <- function1(m = df2, file = "validation/muscle/quality.metrics.2.tsv")
df3 <- function1(m = df3, file = "validation/muscle/quality.metrics.3.tsv")
# 2.2.2. skin not sun exposed
df1 <- function1(m = df1, file = "validation/skin/SKINNS/quality.metrics.1.tsv")
df2 <- function1(m = df2, file = "validation/skin/SKINNS/quality.metrics.2.tsv")
df3 <- function1(m = df3, file = "validation/skin/SKINNS/quality.metrics.3.tsv")
# 2.2.3. skin sun exposed
df1 <- function1(m = df1, file = "validation/skin/SKINS/quality.metrics.1.tsv")
df2 <- function1(m = df2, file = "validation/skin/SKINS/quality.metrics.2.tsv")
df3 <- function1(m = df3, file = "validation/skin/SKINS/quality.metrics.3.tsv")
# 2.2.4. pancreatic islet
df1 <- function1(m = df1, file = "validation/pancreatic.islet/quality.metrics.1.tsv")
df2 <- function1(m = df2, file = "validation/pancreatic.islet/quality.metrics.2.tsv")
df3 <- function1(m = df3, file = "validation/pancreatic.islet/quality.metrics.3.tsv")


# 2.3. add group labels
df1$type <- "ext vs. prediction"
df2$type <- "ext vs. GTEx"
df3$type <- "ext vs. prediction (c.)"

# 2.4. compute MCC
df <- rbind(df1, df2, df3)
df <- function3(m = df)

# 2.5. summarize different metrics
df.BA <- function2(m = df, metric = "Balanced.Accuracy")
df.SN <- function2(m = df, metric = "Sensitivity")
df.SP <- function2(m = df, metric = "Specificity")
df.PR <- function2(m = df, metric = "Precision")
df.mcc <- function2(m = df, metric = "MCC")

# 2.6. make plot
df.plot <- rbind(df.BA, df.SN, df.SP, df.PR, df.mcc)

match.names <- c("GASMED", "PNCREAS", "SKINNS", "SKINS")
names(match.names) <- c("Muscle_Skeletal", "Pancreas", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg")
df.plot$target_tissue_ab <- match.names[df.plot$target_tissue]
df.plot$target_tissue_ab <- factor(df.plot$target_tissue_ab, levels = c("GASMED", "SKINS", "SKINNS", "PNCREAS"))

df.plot$metric <- gsub("Balanced.Accuracy", "Balanced Accuracy", df.plot$metric)
df.plot$metric <- factor(df.plot$metric, levels = rev(c("Balanced Accuracy", 
                                                        "Sensitivity",
                                                        "Precision",
                                                        "Specificity",
                                                        "MCC")))

df.plot$target_tissue_ab <- factor(df.plot$target_tissue_ab,
                                   levels = c("PNCREAS", 
                                              "GASMED",
                                              "SKINNS",
                                              "SKINS"))

pdf("~/public_html/paper_ENTEx/model/validation.ext.datasets.pdf", width = 2, height = 3.5)
ggplot(df.plot[df.plot$type=="ext vs. prediction" &
                 df.plot$metric == "Sensitivity", ]) +
  geom_errorbar(aes(x=target_tissue_ab, 
                    ymin=(value-sd)*100, 
                    ymax=(value+sd)*100), 
                 colour="black",
                 width = 0.4) +
  geom_bar(aes(x=target_tissue_ab, y=value*100, fill=target_tissue_ab),
           stat="identity") +
  # facet_grid(type~.) +
  coord_cartesian(ylim=c(0,100)) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  guides(fill = F) +
  theme(panel.border = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 0.95, hjust=1),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        strip.background = element_blank()) +
  ylab("% of eQTLs correctly identified")
dev.off()


# #-------------------------------------------------
# # 3. make plot filtering for pval Accuracy
# #-------------------------------------------------
# 
# # 3.1. initialize dfs
# df1 <- data.frame(stringsAsFactors = F) # ext. vs. GTEx
# df2 <- data.frame(stringsAsFactors = F) # ext. vs. prediction
# df3 <- data.frame(stringsAsFactors = F) # ext. vs. prediction (only consistent snps)
# 
# # 3.2. read QC dfs from tissues
# # 3.2.1. muscle
# df1 <- function1(m = df1, file = "validation/muscle/quality.metrics.1.tsv", pval = T)
# df2 <- function1(m = df2, file = "validation/muscle/quality.metrics.2.tsv", pval = T)
# df3 <- function1(m = df3, file = "validation/muscle/quality.metrics.3.tsv", pval = T)
# # 3.2.2. skin not sun exposed
# df1 <- function1(m = df1, file = "validation/skin/SKINNS/quality.metrics.1.tsv", pval = T)
# df2 <- function1(m = df2, file = "validation/skin/SKINNS/quality.metrics.2.tsv", pval = T)
# df3 <- function1(m = df3, file = "validation/skin/SKINNS/quality.metrics.3.tsv", pval = T)
# # 3.2.3. skin sun exposed
# df1 <- function1(m = df1, file = "validation/skin/SKINS/quality.metrics.1.tsv", pval = T)
# df2 <- function1(m = df2, file = "validation/skin/SKINS/quality.metrics.2.tsv", pval = T)
# df3 <- function1(m = df3, file = "validation/skin/SKINS/quality.metrics.3.tsv", pval = T)
# # 3.2.4. pancreatic islet
# df1 <- function1(m = df1, file = "validation/pancreatic.islet/quality.metrics.1.tsv", pval = T)
# df2 <- function1(m = df2, file = "validation/pancreatic.islet/quality.metrics.2.tsv", pval = T)
# df3 <- function1(m = df3, file = "validation/pancreatic.islet/quality.metrics.3.tsv", pval = T)
# 
# 
# # 3.3. add group labels
# df1$type <- "ext vs. prediction"
# df2$type <- "ext vs. GTEx"
# df3$type <- "ext vs. prediction (c.)"
# 
# # 3.4. compute MCC
# df <- rbind(df1, df2, df3)
# df <- function3(m = df)
# 
# # 3.5. summarize different metrics
# df.BA <- function2(m = df, metric = "Balanced.Accuracy")
# df.SN <- function2(m = df, metric = "Sensitivity")
# df.SP <- function2(m = df, metric = "Specificity")
# df.PR <- function2(m = df, metric = "Precision")
# df.mcc <- function2(m = df, metric = "MCC")
# 
# # 3.6. make plot
# df.plot <- rbind(df.BA, df.SN, df.SP, df.PR, df.mcc)
# 
# match.names <- c("GASMED", "PNCREAS", "SKINNS", "SKINS")
# names(match.names) <- c("Muscle_Skeletal", "Pancreas", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg")
# df.plot$target_tissue_ab <- match.names[df.plot$target_tissue]
# df.plot$target_tissue_ab <- factor(df.plot$target_tissue_ab, levels = c("GASMED", "SKINS", "SKINNS", "PNCREAS"))
# 
# df.plot$metric <- gsub("Balanced.Accuracy", "Balanced Accuracy", df.plot$metric)
# df.plot$metric <- factor(df.plot$metric, levels = rev(c("Balanced Accuracy", 
#                                                         "Sensitivity",
#                                                         "Precision",
#                                                         "Specificity",
#                                                         "MCC")))
# 
# ggplot(df.plot, aes(x=metric, y=value, color=target_tissue_ab)) +
#   geom_point() +
#   facet_grid(type~.) +
#   coord_flip(ylim=c(0,1)) +
#   scale_color_manual(values = palette) +
#   theme_bw() +
#   theme(panel.border = element_rect(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.text = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust=1),
#         axis.title.x = element_text(size = 13),
#         axis.title.y = element_text(size = 13),
#         plot.title = element_text(size = 13, hjust = .5),
#         strip.background = element_blank(),
#         strip.text.y = element_text(size = 12, angle = 0))
# 
# 
# 
# 
# 
