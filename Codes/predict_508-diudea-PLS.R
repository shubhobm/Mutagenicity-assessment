##
rm(list=ls())
# setwd('D:/Study/My projects/Statchem-Diudea/Codes')
source('RobustQSAR_functions.R')

library(glmnet)
library(ddalpha)
library(plsRglm)
library(ncvreg)
library(randomForest)
library(caret)
library(pROC)
library(parallel)

## Diudea descriptors
basak508 = read.csv("../Data/Basak-descriptors-508.csv")
y508 = as.numeric(basak508[-1,309])
diudea508 = read.csv("../Data/Diudea-descriptors-508.csv")
X508 = diudea508[,-(1:2)] # remove index and character columns
X508 = X508[,-which(apply(X508, 2, function(x) sum(is.na(x)))>0)] # remove NA columns
X508 = as.matrix(X508)

# apply robust scaling
set.seed(11122017)
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
loopfun = function(split){
  set.seed(split*12012017)
  train = sample(1:n, ceiling(.8*n), replace=F)
  ntrain = length(train)
  
  ## PLS regression
  mod508.PLS = plsRglm(y508~., data=df508[train,], nt=10, family="binomial", modele="pls-glm-family", verbose=F)
  pred.PLS = predict(mod508.PLS, newdata=df508[-train,-1])
  cat("Split",split,"done!\n")
  
  as.numeric(auc(roc(y508[-train],pred.PLS)))
}

nsplit = 1e2
err.mat = as.numeric(lapply(1:nsplit, loopfun))
load('err508diudea.Rda')
err.mat508 = cbind(err.mat1[,1],err.mat,err.mat1[,-1])
save(err.mat508, file="err508diudea-final.Rda")

# myform = as.formula(paste0('y508~', paste0(names(df508)[-1], collapse="+")))
# mod508.nn = neuralnet(myform, data=df508, hidden=c(5,3))
