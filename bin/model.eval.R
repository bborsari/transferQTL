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
  
  make_option( c("-m", "--outFile1"), default = NULL,
               help = "Output file for quality metrics." ),
  
  make_option( c("-r", "--outFile2"), default = NULL,
               help = "Output file for roc curve." ),
  
  make_option( c("-i", "--outFile3"), default = NULL,
               help = "Output file for importance features." ),
  
  make_option( c("-o", "--outFolder"), default = NULL,
               help = "Output folder." )
  
  
)

parser <- OptionParser(
  usage = "%prog [options] file",
  option_list=option_list,
  description = "\nComputes quality metrics on the predictive model."
)



arguments <- parse_args(parser, positional_arguments = TRUE)
opt.new <- arguments$options



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



#********
# BEGIN *
#********

# debugging options
# opt.new <- list()
# opt.new$target_tissue <- "Colon_Transverse"
# opt.new$outFile1 <- "plots/quality.metrics.tsv"
# opt.new$outFile2 <- "plots/roc.df.tsv"



# 1. load image of model
load(paste0("output.test/output.objs/", opt.new$target_tissue, ".RData"))

# 2. compute features importance
imp.obj <- varImp(model)
imp.df <- as.data.frame(imp.obj$importance)
imp.df$feature <- rownames(imp.df)
imp.df <- imp.df[order(imp.df$Overall, decreasing = T), ]
imp.df <- imp.df[1:20, ]
imp.df$target_tissue <- opt.new$target_tissue
imp.df$donor_tissue <- opt.new$donor_tissue


write.table(imp.df, file = paste0(opt.new$outFolder, "/", opt.new$outFile3), append = T,
            row.names = F, col.names = F, sep="\t", quote = F)


# 3. predict on training set
set.seed(123)
prediction.training <-  predict(model, newdata = training.set)
print(confusionMatrix(reference = as.factor(training.set$is_eQTL), data = prediction.training, positive = "y"))


# 4. predict on test set
# and save quality metrics of prediction
set.seed(123)
prediction.test <-  predict(model, newdata = test.set)
cm <- confusionMatrix(reference = as.factor(test.set$is_eQTL), data = prediction.test, positive = "y")

tocsv <- data.frame(cbind(t(cm$overall),t(cm$byClass), t(as.numeric(cm$table))))
colnames(tocsv)[19:22] <- c("TN", "FP", "FN", "TP")
rownames(tocsv) <- opt.new$target_tissue
tocsv$donor_tissue <- opt.new$donor_tissue

write.table(tocsv, file = paste0(opt.new$outFolder, "/", opt.new$outFile1), append = T,
            row.names = T, col.names = F, sep="\t", quote = F)


# 5. input df to plot roc curve on training set
selectedIndices <- model$pred$mtry == model$bestTune$mtry
roc.df <- model$pred[selectedIndices, ]
roc.df$target_tissue <- opt.new$target_tissue
roc.df$donor_tissue <- opt.new$donor_tissue

write.table(roc.df, file = paste0(opt.new$outFolder, "/", opt.new$outFile2), append = T,
            row.names = F, col.names = F, sep="\t", quote = F)


