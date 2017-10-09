# Import libraries and data -----------------------------------------------

library(randomForest)
library(foreign)
library(caret)

# Import and data processing ----------------------------------------------

features <- read.arff("boruta_selected_features.arff")
features <- features[sample(nrow(features)),]
indexes = sample(1:nrow(features), size=0.4*nrow(features))
test = features[indexes,]
train = features[-indexes,]

cat("No of missing values : ",sum(is.na(features)),"\n")
cat("No. of variables : ",ncol(features),"\n")
cat("Training cases : ",nrow(train),"\n")
cat("Test cases : ",nrow(test),"\n") 


# Random Forest -----------------------------------------------------------

attach(features)
range_mtry <- seq(100,200,300,400)
range_ntree <- seq(100,200,300,400)
accur <- matrix(0,10)
valid_k <- 5
nnlog <- matrix(0,k)

for (i in 1:valid_k){
  cat("Validation iteration : ",i,"\n")
  for(j in 1:10){
    for(k in 1:10){
      forest <- randomForest(output ~., data = features,
      ntree=range_ntree[k], mtry=range_mtry[j],
      
      
      )
    }
    accur[j] <- unname(confusionMatrix(as.table(forest$confusion[,-3]))$overall[1])
    cat("mtry",range_mtry[j],"\n")
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
