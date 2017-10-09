# A R script to train and test peptides data using neural networks
# Date: 14 September, 2017
# Author: Rakesh Sharma, Deepak Sharma

# Import libraries --------------------------------------------------------

library(foreign) # To read arff files
library(caret)
library(neuralnet)

# Import data -------------------------------------------------------------

features <- read.arff("data_filtered_211.arff")
indexes = sample(1:nrow(features), size=0.2*nrow(features))
test = features[indexes,]
dim(test)  
train = features[-indexes,]
dim(train) 


# Training ----------------------------------------------------------------

train <- as.double(as.matrix(train))
test <- as.double(as.matrix(test))
train_matrix <- matrix(train,nrow = 948,ncol = 211)
test_matrix <- matrix(test,nrow = 237,ncol = 211)
xnam <- paste("train_matrix[",",",1:209,"]", sep="")
fmla <- as.formula(paste("train_matrix[,210]+train_matrix[,211] ~ ", paste(xnam, collapse= "+")))
nnnew <- neuralnet(fmla,data=train_matrix,hidden=c(1),stepmax=1e6)

# Testing -----------------------------------------------------------------

mypredict <- compute(nnnew, test_matrix[,1:209])$net.result
maxidx <- function(arr) {
  return(which(arr == max(arr)))
}
idx <- apply(mypredict, c(1), maxidx)
prediction <- c('0', '1')[idx]
cm <- table(prediction, test_matrix[,211])
print(confusionMatrix(cm))


# Visualization -----------------------------------------------------------

plot(nnnew,show.weights = F,intercept = F)

