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

# combined descriptors
combined95 = read.csv("../Data/Combined-descriptors-95.csv")
y95 = as.numeric(combined95[,2])
X95 = as.matrix(combined95[,-(1:2)])

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
nsplit = 110
err.mat = matrix(0, nrow=nsplit, ncol=6)
for(split in 1:nsplit){
  set.seed(split*12012017)
  train = sample(1:n, ceiling(.8*n), replace=F)
  ntrain = length(train)
  
  ## Principal Component Regression
  Xd = X95[train,]
  depth = depth.projection(X95[train,], X95[train,])
  depth = max(depth) - depth
  for(i in 1:ntrain)
  {
    z = sqrt(sum((Xd[i,  ])^2))
    if(z > ep)
    {
      Xd[i,  ] = depth[i] * (Xd[i,  ]  )/z
    }
  }
  svd95 = svd(Xd)
  npc = min(which(cumsum(svd95$d/sum(svd95$d)) >= .95))
  Xt.train = X95[train,] %*% svd95$v[,1:npc]
  mod95.PCR = lm(X1~., data.frame(cbind(y95[train], Xt.train)))
  pred.PCR = predict(mod95.PCR, newdata=data.frame(cbind(y95[-train], X95[-train,] %*% svd95$v[,1:npc])))
  
  ## PLS regression
  mod95.PLS = plsr(y95~., data=df95[train,], validation="CV")
  pred.PLS = predict(mod95.PLS, df95[-train,])
  err.PLS = min(as.numeric(lapply(1:dim(pred.PLS)[3], function(x) sum((y95[-train] - pred.PLS[,,x])^2))))
  
  ## LS-LASSO
  mod95.lasso = cv.glmnet(X95[train,], y95[train], nfolds=10)
  beta.lasso = as.numeric(coef(mod95.lasso), s="lambda.min")
  pred.lasso = cbind(1, X95[-train,]) %*% beta.lasso
  
  ## LS-SCAD
  mod95.SCAD = cv.ncvreg(X95[train,], y95[train], family="gaussian", penalty="SCAD")
  beta.SCAD = mod95.SCAD$fit$beta[,which.min(mod95.SCAD$cve)]
  pred.SCAD = cbind(1, X95[-train,]) %*% beta.SCAD
  
  ## random forest
  mod95.rf = randomForest(y95~., data=df95[train,])
  pred.rf = predict(mod95.rf, df95[-train,])
  
  ## Gradient boosting, will train using caret
  myControl = trainControl(method="cv", number=5)
  myGrid = expand.grid(n.trees=c(100,500,1e3), interaction.depth=c(1,2), shrinkage=c(.1,.01,.001), n.minobsinnode=c(1,2,5))
  mod95.gbm = train(y95~., data=df95[train,], method="gbm", trControl=myControl, tuneGrid=myGrid,verbose=F)
  pred.gbm = predict(mod95.gbm, df95[-train,])
  
  err.mat[split,] = c(sum((y95[-train] - pred.PCR)^2),
                      err.PLS,
                      sum((y95[-train] - pred.lasso)^2),
                      sum((y95[-train] - pred.SCAD)^2),
                      sum((y95[-train] - pred.rf)^2),
                      sum((y95[-train] - pred.gbm)^2))
  cat("Split",split,"done!\n")
}

save(err.mat, file="err95.Rda")

# myform = as.formula(paste0('y95~', paste0(names(df95)[-1], collapse="+")))
# mod95.nn = neuralnet(myform, data=df95, hidden=c(5,3))
