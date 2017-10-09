# To read arff files
library(foreign)
library(caret) 
features <- read.arff("peptide_data (2).arff")

k <- 3
for(i in 1:k){
  # Spliting in train and test data
  indexes = sample(1:nrow(features), size=0.2*nrow(features))
  test = features[indexes,]
  dim(test)  
  train = features[-indexes,]
  dim(train) 
  
  train <- cbind(train, train$output == '0')
  train <- cbind(train, train$output == '1')
  names(train)[32] <- 'negative'
  names(train)[33] <- 'positive'
  train <- train[-31]
  
  test <- cbind(test, test$output == '0')
  test <- cbind(test, test$output == '1')
  names(test)[32] <- 'negative'
  names(test)[33] <- 'positive'
  test <- test[-31]
  
  dist <- print(as.data.frame(table(test$positive)))
  nnlog <- data.frame(matrix(0, ncol = k, nrow = 1))
  tryCatch({
    nn <- neuralnet(positive + negative ~ 
                      Tiny + Small + Aliphatic + Aromatic + NonPolar + Polar + Charged + 
                      Basic + Acidic + `aIndex(swapped)` + `boman(swapped)` + `charge(swapped)` + 
                      `Hydrophobicity index` + `Alpha and turn propensities` + 
                      `Bulky properties` + `Compositional characteristic index` + 
                      `Local flexibility` + `Electronic properties` + `hmoment(swapped)` + 
                      hydrophobicity_KyteDoolittle + KF1 + KF2 + KF3 + KF4 + KF5 + 
                      KF6 + KF7 + KF8 + KF9 + KF10,
                    data=train, 
                    hidden=c(3,2))
  },
  error = function(e){nnlog[1,k] <- -1},
  warning = function(w){nnlog[1,k] <- -1}, 
  finally={
    # plot(nn)
    mypredict <- compute(nn, test[-31:-32])$net.result
    maxidx <- function(arr) {
      return(which(arr == max(arr)))
    }
    idx <- apply(mypredict, c(1), maxidx)
    prediction <- c('TRUE', 'FALSE')[idx]
    cm <- print(table(prediction, test$positive))
  })
}

confusionMatrix(cm)
