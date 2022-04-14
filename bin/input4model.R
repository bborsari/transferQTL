#!/usr/bin/env /users/rg/bborsari/software/R-3.5.2/bin/Rscript

.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (

  make_option( c("--target_tissue"), default = NULL,
               help = "Tissue in which the model is built (target tissue)." ),
  
  make_option( c("--input_matrix"), default = NULL,
               help = "Input matrix." ),

  make_option( c("--outFile"), default = NULL,
               help = "Output file." ),
  
  make_option( c("--keep_only_tested"), default = F,
               help = "Whether to only keep snps tested in target tissue." ),
  
  make_option( c("--tested_snps_tt"), default = NULL,
               help = "Tested snps in target tissue." ),
  
  make_option( c("--one_gene_eqtls"), default = NULL,
               help = "eQTLs in donor tissue that are associated to only one gene." ),
  
  make_option( c("--entex_rnaseq_m"), default = NULL,
               help = "ENTEx matrix for RNA-seq." ),
  
  make_option( c("--slope_distance"), default = NULL,
               help = "Table w/ slope and tss_distance info for eqtls in donor tissue." ),
  
  make_option( c("--out_repeats_eqtls"), default = NULL,
               help = "eqtls in donor tissue that lie outside repeated regions." ),
  
  make_option( c("--in_cCREs_eqtls"), default = NULL,
               help = "eqtls in donor tissue that fall inside cCREs." ),
  
  make_option( c("--proximal_eqtls"), default = NULL,
               help = "eqtls in donor tissue that are proximal to TSSs." ),
  
  make_option( c("--signal_tables"), default = NULL,
               help = "Tables showing, for eqtls in donor tissue (rows), signals of chromatin features (cols)." ),
  
  make_option( c("--proportion_marked_tissues"), default = NULL,
               help = "Tables showing, for eqtls in donor tissue (rows), proportion of marked tissues for every chromatin feature (cols)." )
  
)

parser <- OptionParser(
  usage = "%prog [options] file",
  option_list=option_list,
  description = "\nPrepares a table to perform a rf model to predict eQTLs."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options



#************
# LIBRARIES *
#************

library(data.table)


#********
# BEGIN *
#********

# debugging options
# opt <- list()
# opt$tissue <- "Stomach"
# opt$donor_tissue <- "Spleen"


# 2. read binary table
# this table contains info of presence/absence of peaks
# of the different chromatin features
# for snps that are eQTLs in the donor tissue (aka spleen)
# + info on whether the snp is an eQTL in the target tissue (aka opt$tissue)
m <- fread(opt$input_matrix, data.table = F)

# 3. filtering step. 1: 
# keep only eQTLs tested in target tissue
if (opt$keep_only_tested) {
  
  tested.snps <- fread(opt$tested_snps_tt, data.table = F, h=F)
  tested.snps <- paste(tested.snps$V1, tested.snps$V2, tested.snps$V3, sep="_")
  m <- m[m$SNP %in% tested.snps, ]
  rm(tested.snps)
  
}


# 4. filtering step. 2: 
# keep eQTLs that control only one gene in donor tissue
one.gene <- fread(opt$one_gene_eqtls, data.table = F, h=F)
m <- m[m$SNP %in% one.gene$V1, ]


# 5. add total number of peaks found per eQTL
m$sum <- apply(m[, 3:ncol(m)], 1, sum)


# 6. read entex expression matrix
# mapping to reference
m.reference <- read.delim(gzfile(opt$entex_rnaseq_m), sep="\t", h=T)
m.reference$gene_id <- gsub("\\..*", "", m.reference$gene_id)


# 7. read slope and tss_distance info in donor tissue
# for selected snps
slope <- read.delim(opt$slope_distance, stringsAsFactors = F)
slope <- slope[slope$SNP %in% m$SNP, ]
slope$gene_id <- gsub("\\..*", "", slope$gene_id)

# 8. compute coefficient of variation for every gene in slope
# and keep only genes with cv != NA

m.reference <- m.reference[m.reference$gene_id %in% slope$gene_id, ]
m.reference$mean <- apply(m.reference[, 2:98], 1, mean)
m.reference$sd <- apply(m.reference[, 2:98], 1, sd)
m.reference$cv <- m.reference$sd / m.reference$mean
m.reference <- m.reference[, c("gene_id", "mean", "sd", "cv")]
m.reference <- m.reference[complete.cases(m.reference), ]


# 9. add for every gene
# info of cv
slope <- merge(slope, m.reference, by = "gene_id", all.x = T)
m <- merge(m, slope, by = "SNP") 
m <- m[complete.cases(m), ] # not all genes associated with eQTLs have cv
m$tss_distance <- abs(m$tss_distance)
m$mean <- NULL
m$sd <- NULL
m$gene_id <- NULL


# 10. add info of whether snp lies out of repeated regions
# unmarked eQTLs fall more often than marked eQTLs
# in repeated regions
repeats.snps <- read.delim(opt$out_repeats_eqtls, h=F)
colnames(repeats.snps) <- c("SNP", "is_out_repeat")
m <- merge(m, repeats.snps, by = "SNP", all.x = T)
m$is_out_repeat <- ifelse(is.na(m$is_out_repeat), 0, m$is_out_repeat)


# 11. add info of whether snp falls in cCRE
eQTLs.cCREs <- read.delim(opt$in_cCREs_eqtls, h=F)
colnames(eQTLs.cCREs) <- c("SNP", "is_cCRE")
m <- merge(m, eQTLs.cCREs, by = "SNP", all.x = T)
m$is_cCRE <- ifelse(is.na(m$is_cCRE), 0, m$is_cCRE)


# 12. add info of number of tissues the eQTL is marked in for a particular target
marking.prop <- fread(opt$proportion_marked_tissues, data.table = F)
colnames(marking.prop)[2:ncol(marking.prop)] <- 
  paste0(colnames(marking.prop)[2:ncol(marking.prop)], "_p")
m <- merge(m, marking.prop, by = "SNP", all.x = T)


# 13. add info of chromatin signal around snp
m.k <- fread(opt$signal_tables, data.table = F)
m.k <- m.k[m.k$SNP %in% m$SNP, ]
colnames(m.k)[2:ncol(m.k)] <- paste0(colnames(m.k)[2:ncol(m.k)], "_k")
colnames(m.k)[1] <- "SNP"
m <- merge(m, m.k, by = "SNP")


# 14. add info on whether the snp is close (+-2 Kb)
# to any annotated TSS
is.prox <- fread(opt$proximal_eqtls, data.table = F)
colnames(is.prox) <- c("SNP", "is_proximal")
m <- merge(m, is.prox, all.x = T, by = "SNP")
m$is_proximal <- ifelse(is.na(m$is_proximal), 0, m$is_proximal)


# 15. convert 0/1 response to factor
m$is_eQTL <- factor(m$is_eQTL, levels = c(0, 1), labels = c("n", "y"))


# 16. save table
write.table(m, file = opt$outFile,
            col.names = T, row.names = F, quote = F, sep = "\t")
