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



# 4. retrieve, for each feature, the strongest (either positive or negative) correlation across models 
selected.features <- c("tss_distance", "cv", "slope", "H3K36me3", "H3K27me3", "H3K27ac", "H3K4me3", "H3K4me1", "H3K9me3",
                       "POLR2A", "POLR2AphosphoS5", "CTCF", "EP300", "ATAC", "DNase")


df <- df[, colnames(df) %in% selected.features]
df.pos <- apply(df, 2, function(x){max(x, na.rm = T)})
df.neg <- apply(df, 2, function(x){min(x, na.rm = T)})
final.cor <- c()
for(i in 1:length(selected.features)) {

  if (abs(df.pos[i]) > abs(df.neg[i])) {

    final.cor <- c(final.cor, df.pos[i])

  } else {

    final.cor <- c(final.cor, df.neg[i])

  }
}


final.cor <- as.data.frame(final.cor)
final.cor$feature <- rownames(final.cor)
final.cor$feature <- gsub("POLR2AphosphoS5", "POLR2ApS5", final.cor$feature)
final.cor$feature <- gsub("tss_distance", "TSS distance", final.cor$feature)
final.cor$feature <- gsub("cv", "tissue specificity", final.cor$feature)

pdf("~/public_html/paper_ENTEx/model/corr.features.simplified.pdf", width = 5, height = 4)
ggplot(final.cor, aes(x=reorder(feature, final.cor), y = final.cor)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust=0.95,vjust=0.5),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 13, hjust = .5)) +
  labs(title = "Correlation with eQTLs activity") +
  ylab(expression(paste("Pearson ", italic("r")))) +
  xlab("Predictive features")
dev.off()
