# Ideation challenge code for InBIx-2017
# Classfication of 3D domain swapping data using random forest
# DBT-BIF team
# Author: Rakesh Saraswat, Deepak Sharma
# Date: 16 September, 2017

# Install packages if not previously installed and import required --------

# install.packages("randomForest")
# install.packages("foreign")
# install.packages("caret")
library(randomForest)
library(foreign)
library(caret)

# Parameter configuration -------------------------------------------------

k <- 5 # K-cross validation

# Import and data processing ----------------------------------------------

# features <- read.arff("boruta_selected_features.csv")
features <- read.arff("boruta_selected_features.arff")
features <- features[sample(nrow(features)),]
indexes = sample(1:nrow(features), size=1/k*nrow(features))
test = features[indexes,]
train = features[-indexes,]

cat("No of missing values : ",sum(is.na(features)),"\n")
cat("No. of variables : ",ncol(features),"\n")
cat("Training cases : ",nrow(train),"\n")
cat("Test cases : ",nrow(test),"\n") 
cat("Positive cases :",sum(features[,657]==1))
cat("Negative cases :",sum(features[,657]==0))


# Random Forest -----------------------------------------------------------

attach(features))

for (i in 1:valid_k){
  cat("Validation iteration : ",i,"\n")
  for(j in 1:10){
    for(k in 1:10){
      forest <- randomForest(output ~., data = features,
                             ntree=range_ntree[k], mtry=range_mtry[j],
                             
                             
      )
    }
    accur[j] <- unname(confusionMatrix(as.table(forest$confusion[,-3]))$overall[1])
    cat("mtry",range_mtry[junname(confusionMatrix(as.table(randomForest(output ~., data = train,xtest=test[,-657],ytest=test[,657],ntree=100, mtry=100)$confusion[,-3]))$overall[1])
],"\n")
    cat("ntree",range_ntree[k],"\n")
  }
  # print(forest)
  # print(importance(forest,type = 2))
  nnlog[i] <- unname(confusionMatrix(as.table(forest$confusion[,-3]))$overall[1])
}
cat("Accuracy after ",valid_k,"cross validation is :",sum(nnlog)/valid_k*100,"%")
print(accur)

confusionMatrix(as.table(tuneRF(y=output,x=features[-657], mtryStart=100, ntreeTry=50, stepFactor=2, improve=0.01, trace=TRUE, plot=TRUE, doBest=TRUE)$confusion[,-3]))$overall[1]
confusionMatrix(as.table(tuneRF(y=output,x=features[-657], mtryStart=50,
                                stepFactor=1.5,  trace=TRUE, plot=TRUE, doBest=TRUE)$confusion[,-3]))$overall[1]
#                                 
#                                 
#                               

attach(train)
unname(confusionMatrix(as.table(randomForest(output ~.,
  data = train, xtest=test[,-657], ytest=test[,657],
  ntree=2000, mtry=2000)$confusion[,-3]))$overall[1])
detach(train)

## Classification:
##data(iris)
##
attach(train)
forest <- randomForest(output ~.,
                       data = train, xtest=test[,-657], ytest=test[,657],
                       ntree=10, mtry=656,importance=TRUE, proximity=TRUE,
                       maxnodes=450,sampsize=c(350,150),classwt=c(0.385,0.614),
                       keep.forest=TRUE,oob.prox=TRUE,oob.times= 1000,
                       replace=TRUE,nodesize=2, do.trace=1
)
forest <- randomForest(output ~.,
           data = train, xtest=test[,-657], ytest=test[,657],
           ntree=10, mtry=656,importance=TRUE, proximity=TRUE,
           maxnodes=40,sampsize=c(350,150),classwt=c(0.385,0.614),
           keep.forest=TRUE,oob.prox=TRUE,oob.times= 1000,
           replace=TRUE,nodesize=2, do.trace=1
           )

unname(confusionMatrix(as.table(forest$confusion[,-3]))[1]$overall[1])
confusionMatrix(as.table(forest$confusion[,-3]))[1]$overall[1]
round(importance(forest), 2)

unname(table(train[,657])[2])
unname(table(train[,657])[1])




# Running 5-cross validation ----------------------------------------------

data <- features
data <- data[sample(nrow(data)),]
folds <- cut(seq(1,nrow(data)),breaks = k,labels = FALSE)
accur_log <- matrix(0,k)

for(i in 1:k){
  indices <- which(folds==i,arr.ind = TRUE)
  test <- data[indices,]
  train <- data[-indices,]
  forest <- randomForest(output ~.,
                         data = train, xtest=test[,-657], ytest=test[,657],
                         ntree=100, mtry=656,importance=TRUE, proximity=TRUE,
                         maxnodes=450,sampsize=c(350,150),classwt=c(0.385,0.614),
                         keep.forest=TRUE,oob.prox=TRUE,oob.times= 1000,
                         replace=TRUE,nodesize=2, do.trace=1
  )
  accur_log[i] <- unname(confusionMatrix(as.table(forest$confusion[,-3]))$overall[1])

  Sys.sleep(1)
}
