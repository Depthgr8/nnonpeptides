# Read arff files ---------------------------------------------------------

library(foreign)
library(neuralnet)
features <- read.arff("peptide_data (2).arff")

k <- 3
for(i in 1:k){
  
  # Spliting in train and test data -----------------------------------------
  
  indexes = sample(1:nrow(features), size=0.2*nrow(features))
  test = features[indexes,]
  dim(test)  
  train = features[-indexes,]
  dim(train) 
  
  # Training NN -------------------------------------------------------------
  
  train <- cbind(train, train$output == '1')
  train <- train[-31]
  names(train)[31] <- 'out'
  
  test <- cbind(test, test$output == '1')
  test <- test[-31]
  names(test)[31] <- 'out'
  
  dist <- print(as.data.frame(table(test$out)))
  nnlog <- data.frame(matrix(0, ncol = k, nrow = 1))
  tryCatch({
    nn <- neuralnet(out ~ 
          Tiny + Small + Aliphatic + Aromatic + NonPolar + Polar + Charged + 
          Basic + Acidic + `aIndex(swapped)` + `boman(swapped)` + `charge(swapped)` + 
          `Hydrophobicity index` + `Alpha and turn propensities` + 
          `Bulky properties` + `Compositional characteristic index` + 
          `Local flexibility` + `Electronic properties` + `hmoment(swapped)` + 
          hydrophobicity_KyteDoolittle + KF1 + KF2 + KF3 + KF4 + KF5 + 
          KF6 + KF7 + KF8 + KF9 + KF10,
        data=train, 
        hidden=c(10),
        threshold=0.01,err.fct="sse",linear.output=FALSE,likelihood=TRUE,
        stepmax=1e+05,rep=1,startweights=NULL,
        learningrate.limit=list(0.1,1.5),
        learningrate.factor=list(minus=0.5,plus=1.5),
        learningrate=0.5,lifesign="minimal",lifesign.step=1000,
        algorithm="backprop",act.fct = "logistic",exclude=NULL,
        constant.weights=NULL
    )
  },
  error = function(e){nnlog[1,k] <- -1},
  warning = function(w){nnlog[1,k] <- -1}, 
  finally={
    # plot(nn)
    mypredict <- compute(nn, test[-31],rep = 1)$net.result
    maxidx <- function(arr) {
      return(which(arr == max(arr)))
    }
    idx <- apply(mypredict, c(1), maxidx)
    prediction <- c('out')[idx]
    cm <- print(table(prediction, test$out))
  })
}


# Testing NN --------------------------------------------------------------



prediction(nn)
print(nn)

table(test$positive,mypredict)
cbind(testnew$subscribed,testnew.nnbp$net.result)
print(testnew.nnbp) 



# plot(nn, information = F)
# gwplot(nn,selected.covariate="Small")
