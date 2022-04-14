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

# 3. rule #1
df <- data.frame(stringsAsFactors = F)
for (tissue in entex.tissues) {
  
  tmp1 <- read.table(paste0("run/", tissue, "/plots/simplified.rule1.tsv"), h=T, sep="\t", stringsAsFactors = F)
  tmp2 <- read.table(paste0("run/", tissue, "/plots/simplified.rule2.tsv"), h=T, sep="\t", stringsAsFactors = F)
  
  m <- cbind(tmp1, tmp2[, c("TN_subset", "total_subset")])
  print(nrow(m))
  
  df <- rbind(df, m)
  
}

df$combo <- paste(df$donor_tissue, df$target_tissue, sep="_")
df$positive_class_p <- df$is_eQTL / df$total
df$positive_class_p_rule1 <- df$TP_subset / df[, 2]
df$positive_class_p_rule2 <- 1- (df$TN_subset / df[, 8])

sorted_combos <- df[order(df$positive_class_p), "combo"]



df.melt <- melt(df[, c("combo", "positive_class_p", "positive_class_p_rule1", "positive_class_p_rule2")])
df.melt$combo <- factor(df.melt$combo, levels = sorted_combos)



pdf("~/public_html/paper_ENTEx/model/simplified.pdf", width = 5.5, height = 4.5)
ggplot() +
  geom_line(data = df.melt, 
            aes(x=combo, y=value, color=variable, group=variable)) +
  geom_vline(xintercept = 362, linetype = "dashed") +
  theme_bw() +
  theme(panel.border = element_rect(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = .5)) +
  xlab("756 donor-target pairs") +
  ylab("fraction of target-tissue eQTLs") +
  guides(color=guide_legend(nrow=3,byrow=TRUE)) +
  scale_color_manual(values = c("positive_class_p" = "#d95f02",
                                  "positive_class_p_rule1" = "#1b9e77",
                                  "positive_class_p_rule2" = "#7570b3"), 
                       labels = c("all donor-tissue eQTLs", 
                                  "low tissue-specif. OR high chromatin activity",
                                  "high tissue-specif. AND low chromatin activity"))
dev.off()




m <- data.frame(set=c("low OR high", "high AND low"),
                value = c(67, 23))

pdf("~/public_html/paper_ENTEx/model/fig.8g.pdf", width = 5, height = 2.5)
ggplot(m, aes(x=set, y=value)) +
  geom_bar(stat = "identity", fill="white", color="black") +
  coord_flip() +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = .5)) +
  ylab("eQTLs transferred to target tissue (%)") +
  labs(title="Summarizing our results with two rules") 

dev.off()