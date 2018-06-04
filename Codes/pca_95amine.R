##
rm(list=ls())
# setwd('d:/Study/My projects/Statchem-Diudea/Codes')
source('RobustQSAR_functions.R')

library(ddalpha)
library(ICSNP)

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
# X95 = as.matrix(scale(X95, mu, scale=F))
which.na = which(is.na(apply(X95,2,var)))
# which.na = which(apply(X95,2,var) < 1e-3)
X95 = X95[,-which.na]
names95 = names(combined95)[-(1:2)][-which.na]
df95 = data.frame(cbind(y95, X95))

n = nrow(X95)
p = ncol(X95)

## Principal Component Analysis
set.seed(04172018)
Xd = X95
depth = depth.projection(X95, X95)
depth = max(depth) - depth
for(i in 1:n)
{
  z = sqrt(sum((Xd[i,  ])^2))
  if(z > ep)
  {
    Xd[i,  ] = depth[i] * (Xd[i,  ]  )/z
  }
}
svd95 = svd(Xd)
df.list = list()
for(i in 1:10){
  V1 = names95[order(abs(svd95$v[,i]), decreasing=T)][1:10]
  V2 = svd95$v[order(abs(svd95$v[,i]), decreasing=T),i][1:10]
  idf = data.frame(cbind(V1,round(V2,2)))
  colnames(idf) = c("Descriptor","Loading")
  df.list[[i]] = idf
}
names(df.list) = paste0("PC",1:10)

## explained variance proportions in data
props = svd95$d/sum(svd95$d)
cumsum(props)

# V1 = names95[order(abs(svd95$v[,1]), decreasing=T)][1:10]
# V2 = svd95$v[order(abs(svd95$v[,1]), decreasing=T),1][1:10]
# data.frame(cbind(V1,round(V2,2)))
# 
# V1 = names95[order(abs(svd95$v[,2]), decreasing=T)][1:10]
# V2 = svd95$v[order(abs(svd95$v[,2]), decreasing=T),2][1:10]
# data.frame(cbind(V1,round(V2,2)))
# 
# V1 = names95[order(abs(svd95$v[,3]), decreasing=T)][1:10]
# V2 = svd95$v[order(abs(svd95$v[,3]), decreasing=T),3][1:10]
# data.frame(cbind(V1,round(V2,2)))
# 
# V1 = names95[order(abs(svd95$v[,4]), decreasing=T)][1:10]
# V2 = svd95$v[order(abs(svd95$v[,4]), decreasing=T),4][1:10]
# data.frame(cbind(V1,round(V2,2)))
