# Ideation challenge code for InBIx-2017
# Classfication of 3D domain swapping data using random forest
# DBT-BIF team
# Author: Rakesh Saraswat, Deepak Sharma
# Date: 16 September, 2017


# Install packages if not installed ---------------------------------------

# install.packages("randomForest")
# install.packages("foreign")
# install.packages("caret")

# Import required libraries -----------------------------------------------

library(randomForest)
library(foreign)
library(caret)

# Parameter configuration -------------------------------------------------

k <- 5                    # K-cross validation
nt <- 10                  # Number of trees
mt <- 656                 # Number of mtry
mn <- 100                 # Maximum no. of nodes

# Import and data processing ----------------------------------------------

features <- read.arff("boruta_selected_features.arff") # Read data file
features <- features[sample(nrow(features)),]   # Shuffling of data records
indexes = sample(1:nrow(features), size=1/k*nrow(features)) # Indices
test = features[indexes,]   # Test data
train = features[-indexes,] # Training data

# Data summary ------------------------------------------------------------

cat("No of missing values : ",sum(is.na(features)),"\n")
cat("No. of variables : ",ncol(features),"\n")
cat("Training cases : ",nrow(train),"\n")
cat("Test cases : ",nrow(test),"\n") 
cat("Positive cases :",sum(features[,657]==1),"\n")
cat("Negative cases :",sum(features[,657]==0),"\n")

# Random Forest -----------------------------------------------------------

# Training case accuracy --------------------------------------------------

tr_forest <- randomForest(output ~., data = train,
          ntree=nt, mtry=mt,importance=TRUE, proximity=TRUE,
          maxnodes=mn,sampsize=c(350,150),classwt=c(0.385,0.614),
          keep.forest=TRUE,oob.prox=TRUE,oob.times= 1000,
          replace=TRUE,nodesize=2, do.trace=1
          )
tra <- unname(confusionMatrix(as.table(tr_forest$confusion[,-3]))$overall[1])
cat("Training accuracy is :",tra,"\n")


# Test case accuracy ------------------------------------------------------

ts_forest <- randomForest(output ~.,
          data = train, xtest=test[,-657], ytest=test[,657],
          ntree=nt, mtry=mt,importance=TRUE, proximity=TRUE,
          maxnodes=mn,sampsize=c(350,150),classwt=c(0.385,0.614),
          keep.forest=TRUE,oob.prox=TRUE,oob.times= 1000,
          replace=TRUE,nodesize=2, do.trace=1
          )
tsa <- unname(confusionMatrix(as.table(ts_forest$confusion[,-3]))$overall[1])
cat("Testing case accuracy is :",tsa,"\n")

# Running 5-cross validation ----------------------------------------------

data <- features
data <- data[sample(nrow(data)),]
folds <- cut(seq(1,nrow(data)),breaks = k,labels = FALSE)
accur_log <- matrix(0,k)

for(i in 1:k){
  indices <- which(folds==i,arr.ind = TRUE)
  test <- data[indices,]
  train <- data[-indices,]
  cat("Running ",i,"/",k,"fold\n")
  cv_forest <- randomForest(output ~.,
              data = train, xtest=test[,-657], ytest=test[,657],
              ntree=nt, mtry=mt,importance=TRUE, proximity=TRUE,
              maxnodes=mn,sampsize=c(350,150),classwt=c(0.385,0.614),
              keep.forest=TRUE,oob.prox=TRUE,oob.times= 1000,
              replace=TRUE,nodesize=2, do.trace=1
              )
  accur_log[i] <- unname(confusionMatrix(as.table(cv_forest$confusion[,-3]))$overall[1])
  Sys.sleep(1)
}
cva <- mean(accur_log)
cat("Average accuracy after ",k,"cross validation is :",cva,"\n")

