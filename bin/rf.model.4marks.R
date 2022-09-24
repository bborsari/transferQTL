#!/usr/bin/env /users/rg/bborsari/software/R-3.5.2/bin/Rscript

.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c("-i", "--input"), default = NULL,
               help = "Input data matrix." ),
  
  make_option( c("-t", "--tissue"), default = NULL,
               help = "Target tissue." ),
  
  make_option( c("-p", "--partition"), type = "numeric", default = 0.7,
               help = "The proportion of input rows used for the training set." ),
  
  make_option( c("-m", "--method"), default = "rf",
               help = "The method to be applied." ),
  
  make_option( c("-o", "--out_folder"), default = NULL,
               help = "Output folder to save the model." )
  
  
  
  
)

parser <- OptionParser(
  usage = "%prog [options] file",
  option_list=option_list,
  description = "\nBuilds a model to predict eQTLs in a given tissue."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options



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
library(doParallel)



#********
# BEGIN *
#********


# # 0. request 10 cores 
# cl <- makePSOCKcluster(10)
# registerDoParallel(cl)

# 1. read input matrix
m <- fread(opt$input, header = T, data.table = F)


# 2. remove snp id
snps <- m$SNP
m$SNP <- NULL

# 3. keep only specific histone marks
sel.cols <- c("is_eQTL", 
              "H3K27ac", 
              "H3K27me3",
              "H3K36me3",        
              "H3K4me1",
              "tss_distance", 
              "slope", 
              "cv",
              "is_out_repeat",
              "is_cCRE",
              "H3K27ac_p",
              "H3K4me1_p",       
              "H3K27me3_p",
              "H3K36me3_p",
              "H3K27ac_k",
              "H3K4me1_k",
              "H3K27me3_k",
              "is_proximal")
m <- m[, colnames(m) %in% sel.cols]


# 4. partition the data into training and test
set.seed(123)
training.examples <- createDataPartition(m$is_eQTL, p = opt$partition, list = FALSE)
training.set  <- m[training.examples, ]
test.set <- m[-training.examples, ]


# 5. define the cross-validation schema
cv <- trainControl(
  method = "cv",
  number = 5,
  summaryFunction=twoClassSummary, 
  classProbs=T,
  savePredictions = T)


# 6. build the model
print("building model")
set.seed(123)
model <- train(is_eQTL ~ ., 
               data = training.set,
               method = opt$method,
               trControl = cv)


# 7. compute variable importance
print(model$results)
print(varImp(model, scale = F))


# 8. save image
save.image(paste0(opt$out_folder, "/", opt$tissue, ".RData"))
# stopCluster(cl)


