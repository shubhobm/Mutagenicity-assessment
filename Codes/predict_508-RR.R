##
rm(list=ls())
# setwd('D:/Study/My projects/Statchem-Diudea/Codes')
source('RobustQSAR_functions.R')

library(glmnet)
library(ddalpha)
library(pls)
library(ncvreg)
library(randomForest)
library(caret)

get.ridge.preds = function(y508,X508){

  # apply robust scaling
  delta = 1e-3
  spa = spatial.median(X508, delta)
  mu = spa$mu
  ep = spa$ep
  sigma.vec = apply(X508, 2, mad)
  X508 = as.matrix(scale(X508, mu, sigma.vec))
  X508 = X508[,-which(is.na(apply(X508,2,var)))]
  df508 = data.frame(cbind(y508, X508))
  
  n = nrow(X508)
  p = ncol(X508)
  
  ## multiple train-test split and predict
  err.vec = rep(0, nsplit)
  for(split in 1:nsplit){
    set.seed(split*12012017)
    train = sample(1:n, ceiling(.8*n), replace=F)
    ntrain = length(train)
    
    ## Ridge regression
    mod508.ridge = cv.glmnet(X508[train,], y508[train], alpha=0, nfolds=10, family="binomial")
    pred.ridge = predict(mod508.ridge, X508[-train,], s="lambda.min")

    err.vec[split] = as.numeric(auc(roc(y508[-train],as.numeric(pred.ridge))))
    cat("Split",split,"done!\n")
  }
  
  err.vec
}

nsplit = 100

# combined descriptors
combined508 = read.csv("../Data/Combined-descriptors-508.csv")
y508 = as.numeric(combined508[,2])
X508 = as.matrix(combined508[,-(1:2)])
err1 = get.ridge.preds(y508, X508)

# basak descriptors
basak508 = read.csv("../Data/Basak-descriptors-508.csv")
y508 = as.numeric(basak508[-1,309])
X508 = as.matrix(basak508[-1,-c(1,309)])
err2 = get.ridge.preds(y508, X508)

## Diudea descriptors
diudea508 = read.csv("../Data/Diudea-descriptors-508.csv")
X508 = diudea508[,-(1:2)] # remove index and character columns
X508 = X508[,-which(apply(X508, 2, function(x) sum(is.na(x)))>0)] # remove NA columns
X508 = as.matrix(X508)
err3 = get.ridge.preds(y508, X508)

err.mat = cbind(err1, err2, err3)
save(err.mat, file="err508RR.Rda")

# myform = as.formula(paste0('y508~', paste0(names(df508)[-1], collapse="+")))
# mod508.nn = neuralnet(myform, data=df508, hidden=c(5,3))
