##
rm(list=ls())
setwd('D:/Study/My projects/Statchem-Diudea/Codes')
source('RobustQSAR_functions.R')

library(glmnet)
library(ddalpha)
library(pls)
library(ncvreg)
library(randomForest)
library(neuralnet)
load('../Data/lta98.rda')

set.seed(11122017)
y = lta98$Y[-1]
X = as.matrix(with(lta98, cbind(ltaTS,ltaTC,lta3D, ltaQC))[-1,])

# combined descriptors
combined95 = read.csv("../Data/Combined-descriptors-95.csv")
y95 = as.numeric(combined95[,2])
X95 = as.matrix(combined95[,-(1:2)])
X95 = X95[,-which(is.na(apply(X95,2,var)))]

# apply robust scaling
delta = 1e-3
spa = spatial.median(X95, delta)
mu = spa$mu
ep = spa$ep
sigma.vec = apply(X95, 2, mad)
X95 = as.matrix(scale(X95, mu, sigma.vec))
df95 = data.frame(cbind(y95, X95))

n = nrow(X95)
p = ncol(X95)

## Principal Component Regression
Xd = X95
# depth = depth.halfspace(X95, X95)
# depth = max(depth) - depth
# for(i in 1:n)
# {
#   z = sqrt(sum((Xd[i,  ])^2))
#   if(z > ep)
#   {
#     Xd[i,  ] = depth[i] * (Xd[i,  ]  )/z
#   }
# }
svd95 = svd(X95)
npc = min(which(cumsum(svd95$d/sum(svd95$d)) >= .95))
Xt = X95 %*% svd95$v[,1:npc]
mod95.PCR = lm(y95~., data.frame(cbind(y95, Xt)))

## PLS regression
mod95.PLS = plsr(y95~., data=df95, validation="CV")

## LS-LASSO
mod95.lasso = cv.glmnet(X95, y95, nfolds=5)
beta.lasso = as.numeric(coef(mod95.lasso), s="lambda.min")

## LS-SCAD
mod95.SCAD = cv.ncvreg(X95, y95, family="gaussian", penalty="SCAD")
beta.SCAD = mod95.SCAD$fit$beta[,which.min(mod95.SCAD$cve)]

## random forest
mod95.rf = randomForest(y95~., data=combined95[,-1])

## Gradient boosting, will train using caret
myControl = trainControl(method="cv", number=5)
myGrid = expand.grid(n.trees=c(100,500,1e3), interaction.depth=c(1,2), shrinkage=c(.1,.01,.001), n.minobsinnode=c(1,2,5))
mod95.gbm = train(y95~., data=df95, method="gbm", trControl=myControl, tuneGrid=myGrid)

myform = as.formula(paste0('y95~', paste0(names(df95)[-1], collapse="+")))
mod95.nn = neuralnet(myform, data=df95, hidden=c(5,3))
