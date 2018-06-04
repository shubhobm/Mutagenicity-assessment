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

get.ridge.preds = function(y95,X95){
  # apply robust scaling
  delta = 1e-3
  spa = spatial.median(X95, delta)
  mu = spa$mu
  ep = spa$ep
  sigma.vec = apply(X95, 2, mad)
  X95 = as.matrix(scale(X95, mu, sigma.vec))
  X95 = X95[,-which(is.na(apply(X95,2,var)))]
  df95 = data.frame(cbind(y95, X95))
  
  n = nrow(X95)
  p = ncol(X95)
  
  ## multiple train-test split and predict
  err.vec = rep(0, nsplit)
  for(split in 1:nsplit){
    set.seed(split*12012017)
    train = sample(1:n, ceiling(.8*n), replace=F)
    ntrain = length(train)
    
    ## Ridge regression
    mod95.ridge = cv.glmnet(X95[train,], y95[train], alpha=0, nfolds=10)
    beta.ridge = as.numeric(coef(mod95.ridge), s="lambda.min")
    pred.ridge = cbind(1, X95[-train,]) %*% beta.ridge
    
    err.vec[split] = sum((y95[-train] - pred.ridge)^2)
    cat("Split",split,"done!\n")
  }
  
  err.vec
}

nsplit = 100

# combined descriptors
combined95 = read.csv("../Data/Combined-descriptors-95.csv")
y95 = as.numeric(combined95[,2])
X95 = as.matrix(combined95[,-(1:2)])
err1 = get.ridge.preds(y95, X95)

# basak descriptors
load('../Data/lta98.rda')
y95 = lta98$Y[-1]
X95 = as.matrix(with(lta98, cbind(ltaTS,ltaTC,lta3D, ltaQC))[-1,])
err2 = get.ridge.preds(y95, X95)

## Diudea descriptors
diudea95 = read.csv("../Data/Diudea-descriptors-95.csv")
X95 = as.matrix(diudea95)
err3 = get.ridge.preds(y95, X95)

err.mat = cbind(err1, err2, err3)
save(err.mat, file="err95RR.Rda")

# myform = as.formula(paste0('y95~', paste0(names(df95)[-1], collapse="+")))
# mod95.nn = neuralnet(myform, data=df95, hidden=c(5,3))
