## Code to generate all model objects for 508 data
rm(list=ls())
# setwd('D:/Study/My projects/Statchem-Diudea/Codes')
source('RobustQSAR_functions.R')

library(glmnet)
library(ddalpha)
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
  
  # Principal Component Regression
  delta = 1e-3
  spa = spatial.median(X, delta)
  mu = spa$mu
  ep = spa$ep
  X = as.matrix(scale(X, center=mu, scale=F))
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
  
  # apply robust scaling to X
  sigma.vec = apply(X, 2, mad)
  X = as.matrix(scale(X, center=F, scale=sigma.vec))
  X = X[,-which(is.na(apply(X,2,var)))]
  
  ## LS-LASSO
  model.list[[3]] = cv.glmnet(X, y, nfolds=10, family="binomial")
  cat("3- Lasso done!\n")
  
  ## LS-SCAD
  model.list[[4]] = cv.ncvreg(X, y, family="binomial", penalty="SCAD")
  # beta.SCAD = mod508.SCAD$fit$beta[,which.min(mod508.SCAD$cve)]
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
    sum(beta.lasso[-1]!=0), # Lasso no. of non-zero coefs
    sum(beta.SCAD[-1]!=0), # SCAD no. of non-zero coefs
    gbm.info # gbm info on tuning parameters
    )
}

# basak descriptors
basak508 = read.csv("../Data/Basak-descriptors-508.csv")
y508 = as.numeric(basak508[-1,309])
X508 = as.matrix(basak508[-1,-c(1,309)])
model508b = train.all.models(X508, y508)
save(model508b, file="modelb_508amine.Rda")

# combined descriptors
combined508 = read.csv("../Data/Combined-descriptors-508.csv")
y508 = as.numeric(combined508[,2])
X508 = as.matrix(combined508[,-(1:2)])
model508c = train.all.models(X508, y508)
save(model508c, file="modelc_508amine.Rda")

## Diudea descriptors
diudea508 = read.csv("../Data/Diudea-descriptors-508.csv")
X508 = diudea508[,-(1:2)] # remove index and character columns
X508 = X508[,-which(apply(X508, 2, function(x) sum(is.na(x)))>0)] # remove NA columns
X508 = as.matrix(X508)
model508d = train.all.models(X508, y508)
save(model508d, file="modeld_508amine.Rda")

# extract info from models
extract.info(model508b)
extract.info(model508d)
extract.info(model508c)


