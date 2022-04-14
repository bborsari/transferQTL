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
  
  make_option( c("--out_file1"), default = NULL,
               help = "Output file (is_eQTL_ext vs. prediction)." ),
  
  make_option( c("--out_file2"), default = NULL,
               help = "Output file (is_eQTL_ext vs. is_eQTL)." ),
  
  make_option( c("--out_file3"), default = NULL,
               help = "Output file (is_eQTL_ext vs. prediction, only consistent eQTLs)." ),
  
  make_option( c("--external_signif"), default = NULL,
               help = "Significant eQTLs from external dataset." ),
  
  make_option( c("--external_tested"), default = NULL,
               help = "Tested eQTLs from external dataset." )
  
)

parser <- OptionParser(
  usage = "%prog [options] file",
  option_list=option_list,
  description = "\nEvaluates the model on a response vector other than GTEx."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt.new <- arguments$options



#************
# LIBRARIES *
#************

library(data.table)
library(caret)
library(randomForest)


#********
# BEGIN *
#********

# debugging options
# opt.new <- list()
# opt.new$donor_tissue <- "Vagina"
# opt.new$target_tissue <- "Muscle_Skeletal"
# opt.new$external_signif <- "/users/rg/bborsari/eQTLs.model.nf/validation/muscle/FUSION_ge_muscle_naive.signif.tsv"
# opt.new$external_tested <- "/users/rg/bborsari/eQTLs.model.nf/validation/muscle/FUSION_ge_muscle_naive.tested.tsv"


# 1. load image of model
load(paste0("/nfs/users2/rg/bborsari/eQTLs.model.nf/run/", opt.new$donor_tissue, "/output.test/output.objs/", opt.new$target_tissue, ".RData"))


# 2. predict on test set
prediction.test <-  predict(model, newdata = test.set)
test.set$SNP <- snps[as.numeric(rownames(test.set))]
test.set$prediction <- prediction.test

# 3. read external tested & keep only snps tested in external dataset
external.tested <- fread(opt.new$external_tested, header = F, data.table = F)
colnames(external.tested) <- "SNP"
test.subset <- test.set[test.set$SNP %in% external.tested$SNP, ]


# 4. add external response info
external.signif <- fread(opt.new$external_signif, header = F, data.table = F)
colnames(external.signif) <- "SNP"
test.subset$is_eQTL_ext <- ifelse(test.subset$SNP %in% external.signif$SNP, "y", "n")


# 5. save results of prediction is_eQTL_ext vs. prediction
cm1 <- confusionMatrix(reference = as.factor(test.subset$is_eQTL_ext), data = as.factor(test.subset$prediction), positive = "y")

tocsv1 <- data.frame(cbind(t(cm1$overall),t(cm1$byClass), t(as.numeric(cm1$table))))
colnames(tocsv1)[19:22] <- c("TN", "FP", "FN", "TP")
tocsv1$target_tissue <- opt.new$target_tissue
tocsv1$donor_tissue <- opt.new$donor_tissue

write.table(tocsv1, file = opt.new$out_file1, append = T,
            row.names = F, col.names = F, sep="\t", quote = F)


# 6. save results of prediction is_eQTL_ext vs. is_eQTL
cm2 <- confusionMatrix(reference = as.factor(test.subset$is_eQTL_ext), data = as.factor(test.subset$is_eQTL), positive = "y")

tocsv2 <- data.frame(cbind(t(cm2$overall),t(cm2$byClass), t(as.numeric(cm2$table))))
colnames(tocsv2)[19:22] <- c("TN", "FP", "FN", "TP")
tocsv2$target_tissue <- opt.new$target_tissue
tocsv2$donor_tissue <- opt.new$donor_tissue

write.table(tocsv2, file = opt.new$out_file2, append = T,
            row.names = F, col.names = F, sep="\t", quote = F)


# 7. save results of prediction is_eQTL_ext vs. prediction, using only 
# snps consistent between GTEx and external dataset
cm3 <- confusionMatrix(reference = as.factor(test.subset[test.subset$is_eQTL_ext == test.subset$is_eQTL, "is_eQTL_ext"]), 
                       data = as.factor(test.subset[test.subset$is_eQTL_ext == test.subset$is_eQTL, "prediction"]), 
                       positive = "y")

tocsv3 <- data.frame(cbind(t(cm3$overall),t(cm3$byClass), t(as.numeric(cm3$table))))
colnames(tocsv3)[19:22] <- c("TN", "FP", "FN", "TP")
tocsv3$target_tissue <- opt.new$target_tissue
tocsv3$donor_tissue <- opt.new$donor_tissue

write.table(tocsv3, file = opt.new$out_file3, append = T,
            row.names = F, col.names = F, sep="\t", quote = F)
