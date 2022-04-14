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
               help = "The method to be applied." )
  
  
  
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



#********
# BEGIN *
#********


# 1. read input matrix
m <- fread(opt$input, header = T, data.table = F)

# 2. remove snp id
snps <- m$SNP
m$SNP <- NULL


# 3. partition the data into training and test
set.seed(123)
training.examples <- createDataPartition(m$is_eQTL, p = opt$partition, list = FALSE)
training.set  <- m[training.examples, ]
test.set <- m[-training.examples, ]


# 4. define the cross-validation schema
cv <- trainControl(
  method = "cv",
  number = 5,
  summaryFunction=twoClassSummary, 
  classProbs=T,
  savePredictions = T)


# 5. build the model
print("building model")
set.seed(123)
model <- train(is_eQTL ~ ., 
               data = training.set,
               method = opt$method,
               trControl = cv
)


# 6. compute variable importance
print(model$results)
print(varImp(model, scale = F))


# 7. save image
save.image(paste0(opt$tissue, ".RData"))


