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

# combined descriptors
basak508 = read.csv("../Data/Basak-descriptors-508.csv")
y508 = as.numeric(basak508[-1,309])
X508 = as.matrix(basak508[-1,-c(1,309)])

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
  
  ## Principal Component Regression
  Xd = X508[train,]
  depth = depth.projection(X508[train,], X508[train,])
  depth = max(depth) - depth
  for(i in 1:ntrain)
  {
    z = sqrt(sum((Xd[i,  ])^2))
    if(z > ep)
    {
      Xd[i,  ] = depth[i] * (Xd[i,  ]  )/z
    }
  }
  svd508 = svd(Xd)
  npc = min(which(cumsum(svd508$d/sum(svd508$d)) >= .95))
  Xt.train = X508[train,] %*% svd508$v[,1:npc]
  mod508.PCR = glm(X1~., family="binomial", data.frame(cbind(y508[train], Xt.train)))
  pred.PCR = predict(mod508.PCR, newdata=data.frame(cbind(y508[-train], X508[-train,] %*% svd508$v[,1:npc])),
                     type="response")
  cat("PCR done\n")
  
  # ## PLS regression
  # mod508.PLS = plsRglm(y508~., data=df508[train,], nt=10, family="binomial", modele="pls-glm-family", verbose=F)
  # pred.PLS = predict(mod508.PLS, newdata=df508[-train,-1])
  # cat("PLS done\n")
  
  ## LS-LASSO
  mod508.lasso = cv.glmnet(X508[train,], y508[train], nfolds=10, family="binomial")
  beta.lasso = as.numeric(coef(mod508.lasso), s="lambda.min")
  pred.lasso = exp(cbind(1, X508[-train,]) %*% beta.lasso)
  pred.lasso = pred.lasso/(1+pred.lasso)
  cat("Lasso done\n")
  
  ## LS-SCAD
  mod508.SCAD = cv.ncvreg(X508[train,], y508[train], family="binomial", penalty="SCAD")
  beta.SCAD = mod508.SCAD$fit$beta[,which.min(mod508.SCAD$cve)]
  pred.SCAD = cbind(1, X508[-train,]) %*% beta.SCAD
  pred.SCAD = pred.lasso/(1+pred.SCAD)
  cat("SCAD done\n")
  
  ## random forest
  mod508.rf = randomForest(as.factor(y508)~., data=df508[train,])
  pred.rf = predict(mod508.rf, df508[-train,], type="prob")[,1]
  cat("rf done\n")
  
  ## Gradient boosting, will train using caret
  myControl = trainControl(method="cv", number=5)
  myGrid = expand.grid(n.trees=c(100,500,1e3), interaction.depth=c(1,2), shrinkage=c(.1,.01,.001), n.minobsinnode=c(1,2,5))
  mod508.gbm = train(as.factor(y508)~., data=df508[train,], method="gbm", trControl=myControl, tuneGrid=myGrid,verbose=F)
  pred.gbm = predict(mod508.gbm, df508[-train,], type="prob")[,1]
  cat("gbm done\n")
  
  cat("Split",split,"done!===========\n")
  c(as.numeric(auc(roc(y508[-train],pred.PCR))),
    # as.numeric(auc(roc(y508[-train],pred.PLS))),
    as.numeric(auc(roc(y508[-train],as.numeric(pred.lasso)))),
    as.numeric(auc(roc(y508[-train],as.numeric(pred.SCAD)))),
    as.numeric(auc(roc(y508[-train],pred.rf))),
    as.numeric(auc(roc(y508[-train],pred.gbm))))
}

nsplit = 1e2
err.mat = mclapply(1:nsplit, loopfun, mc.cores=min(detectCores(),nsplit))
# err.mat = lapply(1:nsplit, loopfun)
err.mat1 = matrix(unlist(err.mat), ncol=5, byrow=T)
save(err.mat1, file="err508basak.Rda")

# myform = as.formula(paste0('y508~', paste0(names(df508)[-1], collapse="+")))
# mod508.nn = neuralnet(myform, data=df508, hidden=c(5,3))
