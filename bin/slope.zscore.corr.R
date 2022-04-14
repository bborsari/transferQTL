library(dplyr)
library(tidyr)


setwd("/users/rg/bborsari/eQTLs.model.nf/application.blood/input.test/tss_distance.slope")

common.eqtls <- fread("common.blood.eQTLs.gtex.ext.tsv", data.table = F, sep="\t", h=F, quote = F, stringsAsFactors = F)
common.eqtls <- common.eqtls$V1

gtex.data <- fread("Whole_Blood.GTEx.slope.tsv", data.table = F, sep="\t", h=T, quote = F, stringsAsFactors = F)
gtex.data <- gtex.data[gtex.data$SNP %in% common.eqtls, ]

ext.data <- fread("../../eQTLs.hg38.tsv", data.table = F, sep="\t", h=F, quote = F, stringsAsFactors = F)
ext.data$SNP <- paste(ext.data$V1, ext.data$V2, ext.data$V3, sep="_")
ext.data <- ext.data[ext.data$SNP %in% common.eqtls, ]
ext.data <- ext.data %>%
  separate(V7, c("gene_id", "zscore"), ";")


m <- merge(gtex.data[, c("SNP", "slope")], ext.data[, c("SNP", "zscore")], by = "SNP")
m$zscore <- as.numeric(m$zscore)

m$direction <- ifelse(m$zscore > 0, 1, -1)
m$abs_zscore <- abs(m$zscore)
m$zscore_log2 <- (log2(m$abs_zscore+1))*m$direction
m$zscore_log2_noP <- (log2(m$abs_zscore))*m$direction
m$zscore_log10 <- (log10(m$abs_zscore+1))*m$direction
m$zscore_log10_noP <- (log10(m$abs_zscore))*m$direction



cor(m$slope, m$zscore_log10, method = "pearson")
