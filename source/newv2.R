# Read arff files ---------------------------------------------------------

library(foreign)
library(neuralnet)
features <- read.arff("peptide_data (2).arff")

sum(is.na(features))
require(Amelia)
missmap(features[-32],main="Missing Data", col=c("red","grey"),legend=TRUE)
summary(features)

boxplot(features)
# hist(features[,1])

# Spliting in train and test data -----------------------------------------

indexes = sample(1:nrow(features), size=0.2*nrow(features))
test = features[indexes,]
dim(test)  
train = features[-indexes,]
dim(train) 

# Processing data -------------------------------------------------------------

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

# Multinomial log linear model --------------------------------------------

##2. Use Multinomial Log Linear models using Neural Networks
train.mlln <- multinom(positive~.,train)
##USe TEST data for testing the trained model
test.mlln <- predict(train.mlln,test)
##Misclassification or Confusion Matrix
table(test$positive,test.mlln)


# Using nnet --------------------------------------------------------------

## 1. Fit a Single Hidden Layer Neural Network using Least Squares
train.nnet<-nnet(positive ~ Tiny,train,size=3,rang=0.07,Hess=FALSE,decay=15e-4,maxit=250)
## Use TEST data for testing the trained model
test.nnet<-predict(train.nnet,test)
## MisClassification Confusion Matrix
table(test$positive,test.nnet)
## One can maximize the Accuracy by changing the "size" while training the neural network. SIZE refers to the number of nodes in the hidden layer.
which.is.max(test.nnet)  ## To Find which row break ties at random (Maximum position in vector)


# Using neuralnet ---------------------------------------------------------

k <- 5
for(i in 1:k){
  
  # Training NN -------------------------------------------------------------
  
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
    mypredict <- compute(nn, test[-31:-32],rep = 1)$net.result
    maxidx <- function(arr) {
      return(which(arr == max(arr)))
    }
    idx <- apply(mypredict, c(1), maxidx)
    prediction <- c('positive', 'negative')[idx]
    cm <- print(table(prediction, test$positive))
  })
}


# Testing NN --------------------------------------------------------------

prediction(nn)
print(nn)
table(test$positive,sum(mypredict[,1]))
cbind(test$positive,test$net.result)
print(test) 
# plot(nn, information = F)
gwplot(nn,selected.covariate="Small")
