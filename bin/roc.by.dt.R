.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

#************
# LIBRARIES *
#************

library(data.table)
library(caret)
library(randomForest)
library(pROC)
library(mlbench)
library(plotROC)
library(MLeval)
library(cowplot)
library(reshape2)


#********
# BEGIN *
#********

# 1. set wd
setwd("/users/rg/bborsari/eQTLs.model.nf/")

# 2. list of entex tissues
entex.tissues <- read.delim("EN-TEx.chromatin.tissues.txt", h=F, sep="\t", stringsAsFactors = F)
entex.tissues <- entex.tissues[, "V1"]

# 3. read palette
pal <- read.delim("/nfs/users2/rg/dgarrido/PhD/projects/GTEx/V8/sqtls/run_1/info/info.txt", 
                  h=T, sep="\t", stringsAsFactors = F)
pal$tissue <- gsub("_-_", "_", pal$tissue)
pal$tissue <- gsub("\\(", "", pal$tissue)
pal$tissue <- gsub("\\)", "", pal$tissue)

palette <- pal$tissue_color_hex
names(palette) <- pal$tissue


# 4. get vector of abbreviated names for tissues
match.names <- c("ADRNLG", "AORTA", "ARTCRN", "ARTTBL", "BREAST", "CLNSGM",
                 "CLNTRN", "ESPGES", "ESPSQE", "ESPMSM", "HRTAA", "HRTLV",
                 "LIVER", "LUNG", "GASMED", "NERVET", "OVARY", "PNCREAS",
                 "PRSTTE", "SKINNS", "SKINS", "PEYERP", "SPLEEN", "STMACH",
                 "TESTIS", "THYROID", "UTERUS", "VAGINA")
names(match.names) <- entex.tissues


# 5. make roc curve plots for each donor tissue
lop <- list()
i <- 1

for (tissue in entex.tissues) {
  
  # tissue is donor-tissue
  m.roc <- fread(paste0("run/", tissue, "/plots/roc.df.tsv"), header = T, data.table = F)
  
  lop[[i]] <- ggplot(m.roc, aes(m=y, d=factor(obs, levels = c("y", "n")), group=target_tissue, color=target_tissue)) +
    geom_roc(n.cuts=0) +
    coord_equal() +
    style_roc(ylab = "TP fraction", xlab = "FP fraction") +
    scale_color_manual(values = palette) +
    theme(axis.title = element_text(size = 13),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 11, angle = 90, hjust=1, vjust = .5),
          plot.title = element_text(size = 15, hjust = .5)) +
    guides(color = F) +
    labs(title = match.names[tissue])
    # +
    # annotate("text", x=0.5, y=0.25, label=paste("avg. AUC =", round(mean(auc.v), 3)))
  
  i <- i+1
  
}

pdf("~/public_html/paper_ENTEx/model/roc.dt.pdf", width = 20, height = 10)
plot_grid(plotlist = lop, nrow = 4, align = "hv")
dev.off()

