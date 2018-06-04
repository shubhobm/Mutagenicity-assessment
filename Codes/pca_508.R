##
rm(list=ls())
# setwd('c:/Study/My projects/Statchem-Diudea/Codes')
source('RobustQSAR_functions.R')

library(ddalpha)

# combined descriptors
combined508 = read.csv("../Data/Combined-descriptors-508.csv")
y508 = as.numeric(combined508[,2])
X508 = as.matrix(combined508[,-(1:2)])

# apply robust scaling
delta = 1e-3
spa = spatial.median(X508, delta)
mu = spa$mu
ep = spa$ep
sigma.vec = apply(X508, 2, mad)
X508 = as.matrix(scale(X508, mu, sigma.vec))
# X508 = as.matrix(scale(X508, mu, scale=F))
which.na = which(is.na(apply(X508,2,var)))
# which.na = which(apply(X508,2,var) < 1e-3)
X508 = X508[,-which.na]
names508 = names(combined508)[-(1:2)][-which.na]
df508 = data.frame(cbind(y508, X508))

n = nrow(X508)
p = ncol(X508)

## Principal Component Analysis
Xd = X508
depth = depth.projection(X508, X508, seed=04172018)
depth = max(depth) - depth
for(i in 1:n)
{
  z = sqrt(sum((Xd[i,  ])^2))
  if(z > ep)
  {
    Xd[i,  ] = depth[i] * (Xd[i,  ]  )/z
  }
}
svd508 = svd(Xd)
df.list = list()
for(i in 1:10){
  V1 = names508[order(abs(svd508$v[,i]), decreasing=T)][1:10]
  V2 = svd508$v[order(abs(svd508$v[,i]), decreasing=T),i][1:10]
  idf = data.frame(cbind(V1,round(V2,2)))
  colnames(idf) = c("Descriptor","Loading")
  df.list[[i]] = idf
}
names(df.list) = paste0("PC",1:10)

## explained variance proportions in data
props = svd508$d/sum(svd508$d)
cumsum(props)

# V1 = names508[order(abs(svd508$v[,2]), decreasing=T)][1:10]
# V2 = svd508$v[order(abs(svd508$v[,2]), decreasing=T),2][1:10]
# data.frame(cbind(V1,V2))
# 
# V1 = names508[order(abs(svd508$v[,3]), decreasing=T)][1:10]
# V2 = svd508$v[order(abs(svd508$v[,3]), decreasing=T),3][1:10]
# data.frame(cbind(V1,V2))
# 
# V1 = names508[order(abs(svd508$v[,4]), decreasing=T)][1:10]
# V2 = svd508$v[order(abs(svd508$v[,4]), decreasing=T),4][1:10]
# data.frame(cbind(V1,V2))
# 
