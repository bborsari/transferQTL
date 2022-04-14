.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c("-d", "--donor_tissue"), default = NULL,
               help = "Donor tissue." ),
  
  make_option( c("-t", "--target_tissue"), default = NULL,
               help = "Target tissue." ),
  
  make_option( c("-o", "--out_file"), default = NULL,
               help = "Output file." )
)

parser <- OptionParser(
  usage = "%prog [options] file",
  option_list=option_list,
  description = "\nCompute correlations between prob. of being an eQTL and features used in the model."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt.new <- arguments$options



#************
# LIBRARIES *
#************

library(data.table)
library(caret)
library(randomForest)
library(pheatmap)

#********
# BEGIN *
#********

# debugging options
# opt.new <- list()
# opt.new$donor_tissue <- "Vagina"
# opt.new$target_tissue <- "Muscle_Skeletal"


# 1. load image of model
load(paste0("/nfs/users2/rg/bborsari/eQTLs.model.nf/run/", opt.new$donor_tissue, "/output.test/output.objs/", opt.new$target_tissue, ".RData"))


# 2. predict on test set
prediction.test <-  predict(model, newdata = test.set, type = "prob")

test.set <- cbind(test.set, prediction.test)
test.set$is_eQTL <- NULL

# 3. compute correlations between probability of being 
# classified as eQTL and all other features
cor.features.response <- cor(test.set)


# 4. save result
cor.features.response <- as.data.frame(cor.features.response["y", !(colnames(cor.features.response) %in% c("n")), drop = F])
cor.features.response$donor_tissue <- opt.new$donor_tissue
cor.features.response$target_tissue <- opt.new$target_tissue

write.table(cor.features.response, file = opt.new$out_file,
            row.names = F, col.names = T, sep="\t", quote = F)

