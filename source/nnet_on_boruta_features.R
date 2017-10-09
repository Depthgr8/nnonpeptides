# Documentation -----------------------------------------------------------

# A R script to train and test peptides data using neural networks
# Date: 14 September, 2017
# Author: Rakesh Sharma, Deepak Sharma

# Import required libraries -----------------------------------------------

library(foreign)
library(caret)
library(nnet)

# Import and data processing ----------------------------------------------

features <- read.arff("boruta_selected_features.arff")
features <- features[sample(nrow(features)),]
cat("No of missing values : ",sum(is.na(features)),"\n")

indexes = sample(1:nrow(features), size=0.4*nrow(features))
test = features[indexes,]
train = features[-indexes,]

cat("No. of variables : ",ncol(features),"\n")
cat("Training cases : ",nrow(train),"\n")
cat("Test cases : ",nrow(test),"\n") 

# Training/testing neural networks ----------------------------------------

# Train neural network 
sink("/dev/null")
train.mlln <- multinom(train$output~.,train,model = TRUE)
sink()

# Test neural network using k-fold cross validation
k <- 10
nnlog <- matrix(0,k)
for (i in 1:k)
{
  cat(i,"/",k,"iteration\n")
  indexes = sample(1:nrow(features), size=0.2*nrow(features))
  test = features[indexes,]
  train = features[-indexes,]
  sink("/dev/null")
  train.mlln <- multinom(train$output~.,train,model = TRUE)
  sink() 
  test.mlln <- predict(train.mlln,test)
  nnlog[i] <- unname(confusionMatrix(table(test$output,test.mlln))$overall[1])
}
cat("Accuracy after ",k,"cross validation is :",sum(nnlog)/k*100,"%")
