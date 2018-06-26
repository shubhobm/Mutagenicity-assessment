## Code to generate all model objects for 508 data
rm(list=ls())
# setwd('D:/Study/My projects/Statchem-Diudea/Codes')
source('RobustQSAR_functions.R')

library(glmnet)
library(ddalpha)
library(pls)
library(ncvreg)
library(randomForest)
library(caret)

## function to train all models
train.all.models = function(X,y){
  
  n = nrow(X)
  p = ncol(X)
  
  Xy = data.frame(cbind(y,X))
  names(Xy) = c("y",paste0("X",1:p))
  model.list = list()
  
  # apply robust scaling to X
  delta = 1e-3
  spa = spatial.median(X, delta)
  mu = spa$mu
  ep = spa$ep
  sigma.vec = apply(X, 2, mad)
  X = as.matrix(scale(X, center=mu, scale=sigma.vec))
  X = X[,-which(is.na(apply(X,2,var)))]
  
  # Principal Component Regression
  Xd = X
  depth = depth.projection(X, X)
  depth = max(depth) - depth
  for(i in 1:n)
  {
    z = sqrt(sum((Xd[i,  ])^2))
    if(z > ep)
    {
      Xd[i,  ] = depth[i] * (Xd[i,  ]  )/z
    }
  }
  svdmod = svd(Xd)
  npc = min(which(cumsum(svdmod$d/sum(svdmod$d)) >= .95))
  Xt = X %*% svdmod$v[,1:npc]
  model.list[[1]] = lm(y~., data.frame(cbind(y, Xt)))
  cat("1- PCR done!\n")
  
  ## PLS regression
  model.list[[2]] = plsr(y~., data=Xy, ncomp=75)
  cat("2- PLS done!\n")
  
  ## LS-LASSO
  model.list[[3]] = cv.glmnet(X, y, nfolds=10)
  cat("3- Lasso done!\n")
  
  ## LS-SCAD
  model.list[[4]] = cv.ncvreg(X, y, family="gaussian", penalty="SCAD")
  # beta.SCAD = mod95.SCAD$fit$beta[,which.min(mod95.SCAD$cve)]
  cat("4- SCAD done!\n")
  
  ## random forest
  model.list[[5]] = randomForest(y~., data=Xy)
  cat("5- rf done!\n")
  
  ## Gradient boosting, will train using caret
  myControl = trainControl(method="cv", number=5)
  myGrid = expand.grid(n.trees=c(100,500,1e3), interaction.depth=c(1,2), shrinkage=c(.1,.01,.001), n.minobsinnode=c(1,2,5))
  model.list[[6]] = train(y~., data=Xy, method="gbm", trControl=myControl, tuneGrid=myGrid,verbose=F)
  cat("6- gbm done!\n-----\n")
  
  names(model.list) = c("PCR","PLS","Lasso","SCAD","rf","gbm")
  model.list
}

extract.info = function(model.list){
  
  beta.lasso = as.numeric(coef(model.list[[3]]), s="lambda.min")
  beta.SCAD = model.list[[4]]$fit$beta[,which.min(model.list[[4]]$cve)]
  gbm.info = with(model.list[[6]]$finalModel,
                  c(n.trees, interaction.depth, shrinkage, n.minobsinnode))
  c(length(model.list[[1]]$coef)-1, # PCR no. of PCs
    model.list[[2]]$ncomp, # PLS no. of comps
    sum(beta.lasso[-1]!=0), # Lasso no. of non-zero coefs
    sum(beta.SCAD[-1]!=0), # SCAD no. of non-zero coefs
    gbm.info # gbm info on tuning parameters
    )
}

# basak descriptors
load('../Data/lta98.rda')
y95 = lta98$Y[-1]
X95 = as.matrix(with(lta98, cbind(ltaTS,ltaTC,lta3D, ltaQC))[-1,])
model95b = train.all.models(X95, y95)
save(model95b, file="modelb_95amine.Rda")

# combined descriptors
combined95 = read.csv("../Data/Combined-descriptors-95.csv")
y95 = as.numeric(combined95[,2])
X95 = as.matrix(combined95[,-(1:2)])
model95c = train.all.models(X95, y95)
save(model95c, file="modelc_95amine.Rda")

## Diudea descriptors
load('../Data/lta98.rda')
y95 = lta98$Y[-1]
diudea95 = read.csv("../Data/Diudea-descriptors-95.csv")
X95 = as.matrix(diudea95)
model95d = train.all.models(X95, y95)
save(model95d, file="modeld_95amine.Rda")

# extract info from models
extract.info(model95b)
extract.info(model95d)
extract.info(model95c)


