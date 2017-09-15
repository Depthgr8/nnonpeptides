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

train.mlln <- multinom(train$output~.,train,model = TRUE)
test.mlln <- predict(train.mlln,test)
print(confusionMatrix(table(test$output,test.mlln)))

