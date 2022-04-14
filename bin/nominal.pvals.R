.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  
  make_option( c("-t", "--target_tissue"), default = NULL,
               help = "Target tissue." ),
  
  make_option( c("-d", "--donor_tissue"), default = NULL,
               help = "Donor tissue." ),
  
  make_option( c("-o", "--outFile"), default = NULL,
               help = "Output file." )
  
  
)

parser <- OptionParser(
  usage = "%prog [options] file",
  option_list=option_list,
  description = "\nRetrieves nominal p-vals for TN, TP, FN and FP."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt.new <- arguments$options



#************
# LIBRARIES *
#************

library(data.table)
library(caret)


#********
# BEGIN *
#********

# 1. load image of model
load(paste0("/nfs/users/rg/bborsari/eQTLs.model.nf/run/", 
            opt.new$donor_tissue, 
            "/output.test/output.objs/",
            opt.new$target_tissue,
            ".RData"))

# 2. predict on test set
set.seed(123)
prediction.test <-  predict(model, newdata = test.set)

# 3. add the prediction's result to test-set 
test.set$prediction <- prediction.test

# 4. subset snps vector for those snps included in test set, and add variant_id
snps <- snps[as.numeric(rownames(test.set))]
test.set$id <- snps

# 5. read table with info about nominal pval and slope
info <- fread(paste0("/nfs/users/rg/bborsari/eQTLs.model.nf/tested.pv.slope/", opt.new$target_tissue, ".tested.pv.slope.tsv"), header = F, data.table = F)
colnames(info) <- c("variant_id", "pval_nominal", "slope", "slope_se")

# 6. inside info, snps may be repeated 
# b/c they may have been tested for association w/ more than 1 gene
# also, untested snps will be removed
info <- info[info$variant_id %in% snps, ]

# 7. b/c of 6., each snp is counted only once
# considering the smallest p-val
info <- info[order(info$variant_id, info$pval_nominal), ] 
info <- info[ !duplicated(info$variant_id), ] 


# 8. get sets of TN, FP, FN and TP
TN.set <- test.set[test.set$is_eQTL == "n" & test.set$prediction == "n", "id"]
FP.set <- test.set[test.set$is_eQTL == "n" & test.set$prediction == "y", "id"]
FN.set <- test.set[test.set$is_eQTL == "y" & test.set$prediction == "n", "id"]
TP.set <- test.set[test.set$is_eQTL == "y" & test.set$prediction == "y", "id"]

# 9. retrieve for these 4 different sets of snps the nominal pval
TN.pv <- info[info$variant_id %in% TN.set, ]
TN.pv$set <- "TN"
FP.pv <- info[info$variant_id %in% FP.set, ]
FP.pv$set <- "FP"
FN.pv <- info[info$variant_id %in% FN.set, ]
FN.pv$set <- "FN"
TP.pv <- info[info$variant_id %in% TP.set, ]
TP.pv$set <- "TP"

# 10. generate a unique df
df <- rbind(TN.pv, FP.pv, FN.pv, TP.pv)
df$donor_tissue <- opt.new$donor_tissue
df$target_tissue <- opt.new$target_tissue


# 11. save output
write.table(df, file = opt.new$outFile, append = T, row.names = F, col.names = F, sep="\t", quote = F)



# df$set <- factor(df$set, levels = c("TN", "FP", "FN", "TP"))
# 
# pdf("~/public_html/paper_ENTEx/model/nominal.pvals.pdf",
#     width = 4, height = 4)
# ggplot(df, aes(x=set, y=pval_nominal)) +
#   geom_boxplot() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 15),
#         panel.border = element_rect(), 
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(), 
#           axis.line = element_line(colour = "black"),)
# dev.off()
# 
