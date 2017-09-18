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
range <- c(1:10)
out_range <- matrix(0,10)
k <- 1
nnlog <- matrix(0,k)
for (i in 1:k){
  cat("iteration",i)
  for(j in 1:4)
  {
    forest <- randomForest(output ~., data = features,ntree=1000, mtry=range[j])
    out_range[j] <- unname(confusionMatrix(as.table(forest$confusion[,-3]))$overall[1])
    cat("mtry",range[j])
  }
  # print(forest)
  # print(importance(forest,type = 2))
  nnlog[i] <- unname(confusionMatrix(as.table(forest$confusion[,-3]))$overall[1])
}
cat("Accuracy after ",k,"cross validation is :",sum(nnlog)/k*100,"%")
