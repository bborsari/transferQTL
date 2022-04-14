.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  
  make_option( c("-d", "--donor_tissue"), default = NULL,
               help = "Target tissue." ),
  
  make_option( c("-t", "--target_tissue"), default = NULL,
               help = "Target tissue." ),
  
  make_option( c("--out_file1"), default = NULL,
               help = "Output file for rule #1." ),
  
  make_option( c("--out_file2"), default = NULL,
               help = "Output file for rule #2." )
  
  
)


parser <- OptionParser(
  usage = "%prog [options] file",
  option_list=option_list,
  description = "\nBuilds a simplified decision tree"
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options



#************
# LIBRARIES *
#************

library(data.table)
library(caret)


#************
# FUNCTIONS *
#************

function1 <- function(feature, df, feature_imp) {
  
  mean_y <- mean(df[df$is_eQTL == "y", feature])
  mean_n <- mean(df[df$is_eQTL == "n", feature])
  
  if (mean_y < mean_n) {
    
    df[, paste0("pred_", feature)] <- ifelse(df[, feature] < mean_y, 1, 0)*feature_imp
    
  } else {
    
    df[, paste0("pred_", feature)] <- ifelse(df[, feature] < mean_n, 0, 1)*feature_imp
    
  }                 
  
  return(df)
  
}


#********
# BEGIN *
#********

# debugging opts
# opt <- list()
# opt$donor_tissue <- "Spleen"
# opt$target_tissue <- "Breast_Mammary_Tissue"

# 1. read input file for model
m <- fread(paste0("~/eQTLs.model.nf/run/", opt$donor_tissue, "/output.test/input.tables/", opt$target_tissue, ".input4model.tsv"), 
           h=T, data.table = F)

# 2. evaluate rule #1
m.pos <- m[m$cv <= 0.8 | m$sum > 3, ] #| (m$tss_distance < 10000 & m$sum > 1 & m$cv > 5), ]
m.pos$prediction <- "y"
nrow(m.pos[m.pos$prediction == m.pos$is_eQTL, ]) / nrow(m.pos)

df.pos <- data.frame(TP_subset = nrow(m.pos[m.pos$prediction == m.pos$is_eQTL, ]),
                     total_subset = nrow(m.pos),
                     is_eQTL = nrow(m[m$is_eQTL == "y", ]),
                     total = nrow(m))
df.pos$donor_tissue <- opt$donor_tissue
df.pos$target_tissue <- opt$target_tissue

write.table(df.pos, file = opt$out_file1, append = T,
            row.names = F, col.names = F, sep="\t", quote = F)


# m.neg <- m[m$cv > 5 & ((m$sum < 1 & m$tss_distance < 10000) | (m$sum < 2 & m$H3K27me3 == 1)), ]
m.neg <- m[m$cv > 5 & m$sum == 0, ]
m.neg$prediction <- "n"
nrow(m.neg[m.neg$prediction == m.neg$is_eQTL, ]) / nrow(m.neg)

df.neg <- data.frame(TN_subset = nrow(m.neg[m.neg$prediction == m.neg$is_eQTL, ]),
                     total_subset = nrow(m.neg),
                     isNot_eQTL = nrow(m[m$is_eQTL == "n", ]),
                     total = nrow(m))
df.neg$donor_tissue <- opt$donor_tissue
df.neg$target_tissue <- opt$target_tissue

write.table(df.neg, file = opt$out_file2, append = T,
            row.names = F, col.names = F, sep="\t", quote = F)


# 2. keep only tested snps
# tested.snps <- fread(paste0("~/eQTLs.model.nf/run/", opt$donor_tissue, "/input.test/tested/unique/", opt$target_tissue, ".bed"), 
#                      data.table = F, h=F)
# tested.snps <- paste(tested.snps$V1, tested.snps$V2, tested.snps$V3, sep="_")
# m <- m[m$SNP %in% tested.snps, ]
# snps <- m$SNP
# m$SNP <- NULL

# m.safe <- m

# # 3. check most important features
# # 3.1. cv
# m <- function1(feature = "cv", df = m, feature_imp = 0.8)
# 
# # 3.2 tss_distance
# m <- function1(feature = "tss_distance", df = m, feature_imp = 0.25)
# 
# # 3.3 slope
# m <- function1(feature = "slope", df = m, feature_imp = 0.4)
# 
# # 3.4 H3K27me3_k
# m <- function1(feature = "H3K27me3_k", df = m, feature_imp = 0.1)
# 
# # 3.5. H3K9me3_k
# m <- function1(feature = "H3K9me3_k", df = m, feature_imp = 0.1)
# 
# # 3.6. H3K4me1_k
# m <- function1(feature = "H3K4me1_k", df = m, feature_imp = 0.1)
# 
# # 3.7. CTCF_k
# m <- function1(feature = "CTCF_k", df = m, feature_imp = 0.1)
# 
# # 3.8. H3K27ac_k
# m <- function1(feature = "H3K27ac_k", df = m, feature_imp = 0.1)
# 
# # 3.9. POLR2A_k
# m <- function1(feature = "POLR2A_k", df = m, feature_imp = 0.1)
# 
# # 3.10. ATAC_k
# m <- function1(feature = "ATAC_k", df = m, feature_imp = 0.1)
# 
# # 3.11. DNase_k
# m <- function1(feature = "DNase_k", df = m, feature_imp = 0.1)
#
# # 3.12 H3K27ac_p
# m <- function1(feature = "H3K27ac_p", df = m, feature_imp = 0.1)
# 
# # 3.13 H3K4me3_p
# m <- function1(feature = "H3K4me3_p", df = m, feature_imp = 0.1)
# 
# # 3.14 H3K4me1_p
# m <- function1(feature = "H3K4me1_p", df = m, feature_imp = 0.1)
# 
# # 3.15 H3K27me3_p
# m <- function1(feature = "H3K27me3_p", df = m, feature_imp = 0.1)
# 
# # 3.16 H3K36me3_p
# m <- function1(feature = "H3K36me3_p", df = m, feature_imp = 0.1)
# 
# # 3.17 H3K9me3_p
# m <- function1(feature = "H3K9me3_p", df = m, feature_imp = 0.1)
# 
# # 3.18 CTCF_p
# m <- function1(feature = "CTCF_p", df = m, feature_imp = 0.1)
# 
# # 3.19 POLR2A_p
# m <- function1(feature = "POLR2A_p", df = m, feature_imp = 0.1)
# 
# # 3.20 POLR2AphosphoS5_p
# m <- function1(feature = "POLR2AphosphoS5_p", df = m, feature_imp = 0.1)
# 
# # 3.21 EP300_p
# m <- function1(feature = "EP300_p", df = m, feature_imp = 0.1)
# 
# # 3.22 ATAC_p
# m <- function1(feature = "ATAC_p", df = m, feature_imp = 0.1)
# 
# # 3.23 DNase_p
# m <- function1(feature = "DNase_p", df = m, feature_imp = 0.1)


# # 4. summarize predictive values
# m$pred <- apply(m[, c(39:48)], 1, sum)
# boxplot(m[m$is_eQTL == "y", "pred"], 
#         m[m$is_eQTL == "n", "pred"])
#                   
# m$prediction <- ifelse(m$pred > 0.5, "y", "n")
# nrow(m[m$is_eQTL == m$prediction, ]) / nrow(m)
# nrow(m[m$is_eQTL == "y" & m$prediction == "y", ]) / nrow(m[m$is_eQTL == "y", ]) # 80% sensitivity
# nrow(m[m$is_eQTL == "n" & m$prediction == "n", ]) / nrow(m[m$is_eQTL == "n", ]) # 30% specificity
# 
# confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")


