.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(data.table)
library(caret)
library(randomForest)
library(pheatmap)
library(RColorBrewer)

#********
# BEGIN *
#********

# 1. set workding directory
setwd("/users/rg/bborsari/eQTLs.model.nf/")


# 2. entex tissues
entex.tissues <- read.table("EN-TEx.chromatin.tissues.txt", h=F, sep="\t", stringsAsFactors = F)
entex.tissues <- entex.tissues$V1

# 3. get list of all tested features
# use stomach since it was characterized for all features in at least one donor
features <- read.table("run/Spleen/plots/corr.features.Spleen.Stomach.tsv", h=T, sep="\t")
features <- colnames(features)

# 3. 
df <- data.frame(stringsAsFactors = F)

for (donor_t in entex.tissues) {
  
  for (target_t in entex.tissues[entex.tissues != donor_t]) {
    
    tmp <- read.table(paste0("run/", donor_t, "/plots/corr.features.", donor_t, ".", target_t, ".tsv"), 
                      stringsAsFactors = F, h=T)
    missing_features <- setdiff(features, colnames(tmp))
    
    if (length(missing_features) > 0) {
      
      for ( f in missing_features ) {
        
        tmp[, f] <- NA
        
      }
      
    }
    
    df <- rbind(df, tmp)
    
  }
  
}


discard.tissues <- c("Lung", "Esophagus_Gastroesophageal_Junction",
                     "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg",
                     "Vagina", "Artery_Tibial", "Ovary", "Artery_Aorta", "Artery_Coronary", "Liver",
                     "Uterus", "Adrenal_Gland")

df2 <- df[!(df$target_tissue %in% discard.tissues), 1:35]
df2.anno <- df[!(df$target_tissue %in% discard.tissues), 37:38]




# read palette
pal <- read.delim("/nfs/users2/rg/dgarrido/PhD/projects/GTEx/V8/sqtls/run_1/info/info.txt", 
                  h=T, sep="\t", stringsAsFactors = F)
pal$tissue <- gsub("_-_", "_", pal$tissue)
pal$tissue <- gsub("\\(", "", pal$tissue)
pal$tissue <- gsub("\\)", "", pal$tissue)
pal <- pal[pal$tissue %in% entex.tissues, ]

palette <- pal$tissue_color_hex
names(palette) <- pal$tissue


match.names <- c("ADRNLG", "AORTA", "ARTCRN", "ARTTBL", "BREAST", "CLNSGM",
                 "CLNTRN", "ESPGES", "ESPSQE", "ESPMSM", "HRTAA", "HRTLV",
                 "LIVER", "LUNG", "GASMED", "NERVET", "OVARY", "PNCREAS",
                 "PRSTTE", "SKINNS", "SKINS", "PEYERP", "SPLEEN", "STMACH",
                 "TESTIS", "THYROID", "UTERUS", "VAGINA")
names(match.names) <- sort(unique(pal$tissue))


names(palette) <- match.names
df2.anno$donor_tissue <- match.names[df2.anno$donor_tissue]
df2.anno$target_tissue <- match.names[df2.anno$target_tissue]

ann_colors = list(donor_tissue = palette,
                  target_tissue = palette)

colnames(df2) <- gsub("cv", "tissue_specificity", colnames(df2))

pdf("~/public_html/paper_ENTEx/model/corr.features.pdf", width = 8, height = 5)
out <- pheatmap(df2, clustering_method = "ward.D2",
         clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(12),
         breaks = c(-0.5, -0.4, -0.3, -0.2, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.2, 0.3),
         show_rownames = F,
         main = "Pearson r",
         # annotation_row = df2.anno,
         # annotation_colors = ann_colors,
         border_color = NA,
         fontsize = 7,
         fontsize_row = 10,
         fontsize_col = 10)
         # color = c("#2166AC", "#4B8EC1", "#7CB5D5", "#BBD9E9", "#E0ECF2", "white", "#F19A76", "#D65C4B", "#B2182B"),
         # breaks = c(-0.5, -0.4, -0.3, -0.2, -0.1, -0.05,  0.05, 0.1, 0.2, 0.3))
dev.off()



# extract clustering rows 
cls <- sort(cutree(out$tree_row, k=2))
length(cls[cls==2])

cl2 <- df2.anno[names(cls[cls==2]), ]
cl1 <- df2.anno[names(cls[cls==1]), ]

# check number of samples per donor tissue
pal <- read.delim("/nfs/users2/rg/dgarrido/PhD/projects/GTEx/V8/sqtls/run_1/info/info.txt", 
                  h=T, sep="\t", stringsAsFactors = F)
pal$tissue <- gsub("_-_", "_", pal$tissue)
pal$tissue <- gsub("\\(", "", pal$tissue)
pal$tissue <- gsub("\\)", "", pal$tissue)

pal <- pal[pal$tissue %in% entex.tissues, ]
pal$abb <- match.names[pal$tissue]

cl1 <- merge(cl1, pal[, c("abb", "nb_samples_analyzed")], by.x = "target_tissue",
             by.y = "abb")
cl1$cluster <- "1"
cl2 <- merge(cl2, pal[, c("abb", "nb_samples_analyzed")], by.x = "target_tissue",
             by.y = "abb")
cl2$cluster <- "2"


cl <- rbind(cl1, cl2)
cl$cluster <- factor(cl$cluster, levels = c(2, 1))



pdf("~/public_html/paper_ENTEx/model/corr.features.boxplot.nsamples.pdf", width = 4, height = 4,
    useDingbats = F)
ggplot(cl, aes(x=cluster, y=nb_samples_analyzed)) +
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  theme(panel.border = element_rect(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = .5)) +
  scale_x_discrete(labels = rev(c("top", "bottom"))) +
  ylab("target-tissue sample size")
dev.off()
