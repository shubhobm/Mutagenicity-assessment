setwd("d:/Study/My projects/StatChem_Robust")
rm(list=ls());
source('RobustQSAR_functions.R')
library(rrcov)
library(fda.usc)
library(ldr)

## prediction function
get.preds = function(X,y,varpct=90){
  # SDR
  set.seed(03212016)
  pred1 = rep(0 ,nrow(X))
  for(i in 1:nrow(X)){
    imod = sir(X=X[-i,], y=y[-i], varpct=varpct, nslices=5)  
    pred1[i] = predict.rsir(imod, t(as.matrix(X[i,])))
  }
  
  # robust SDR
  set.seed(03212016)
  pred2 = rep(0 ,nrow(X))
  for(i in 1:nrow(X)){
    imod = rsir(X=X[-i,], y=y[-i], varpct=varpct, nslices=5)  
    pred2[i] = predict.rsir(imod, t(as.matrix(X[i,])))
  }
  
  # principal component regression
  set.seed(03212016)
  pred3 = rep(0 ,nrow(X))
  for(i in 1:nrow(X)){
    iPC = PcaClassic(X[-i,])
    eigen.cumsum = cumsum(iPC@eigenvalues)
    d = min(which(eigen.cumsum >= (varpct/100) * eigen.cumsum[length(eigen.cumsum)]))
    iS = iPC@scores[,1:d]
    imod = lm.fit(x=as.matrix(iS), y=y[-i])
    
    pred3[i] = sum((X[i,] %*% iPC@loadings)[1:d] * imod$coef)
  }
  
  # robust principal component regression
  set.seed(03212016)
  pred4 = rep(0 ,nrow(X))
  for(i in 1:nrow(X)){
    iPC = PcaRank(X[-i,])
    eigen.cumsum = cumsum(iPC@eigenvalues)
    d = min(which(eigen.cumsum >= (varpct/100) * eigen.cumsum[length(eigen.cumsum)]))
    iS = iPC@scores[,1:d]
    imod = lm.fit(x=as.matrix(iS), y=y[-i])
    
    pred4[i] = sum((X[i,] %*% iPC@loadings)[1:d] * imod$coef)
  }
  
  c(sum((pred1 - ta98y)^2),sum((pred2 - ta98y)^2),
    sum((pred3 - ta98y)^2),sum((pred4 - ta98y)^2))  
}

## compare with principal component regression


## application on 95 amine data
load('lta98.rda')
ta98X = with(lta98, cbind(ltaTS,ltaTC,lta3D,ltaQC))[-1,]
#ta98X = logt(ta98X)
#ta98X = scale(ta98X)
ta98X = log(-floor(min(ta98X)) + ta98X)
ta98y = as.numeric(lta98$Y[-1])

outmat = rbind(get.preds(logt(lta98$ltaTS[-1,]), ta98y),
               get.preds(logt(with(lta98, cbind(ltaTS, ltaTC))[-1,]), ta98y),
               get.preds(logt(with(lta98, cbind(ltaTS, ltaTC, lta3D))[-1,]), ta98y),
               get.preds(ta98X, ta98y),
               get.preds(logt(lta98$ltaTS[-1,]), ta98y),
               get.preds(logt(lta98$ltaTC[-1,]), ta98y),
               get.preds(logt(lta98$lta3D[-1,]), ta98y),
               get.preds(logt(lta98$ltaQC[-1,]), ta98y))
round(outmat,digits=0)

##### 95 amine data
n = nrow(ta98X)
p = ncol(ta98X)
minus = c(1,194,195,196)
varnames95 = names(read.delim('varnames95.txt'))[-c(minus,p+5)]

delta = 1e-3
data = ta98X
y = data
spa = spatial.median(data, delta)
mu = spa$mu
ep = spa$ep
tt = matrix(mu, n, p, byrow=TRUE)
data = data-tt
depth = mdepth.HS(data, data)$dep
depth = max(depth) - depth
for(i in 1:n)
{
  z = sqrt(sum((data[i,  ])^2))
  y[i,  ] = 0 * data[i,  ]
  if(z > ep)
  {
    y[i,  ] = depth[i] * (data[i,  ]  )/z
  }
}
svd95 = svd(y)
loading95 = svd95$v
  
PCi.vals = matrix(0, nrow=p, ncol=min(n,p))
PCi.names = matrix("0", nrow=p, ncol=min(n,p))

for(i in 1:min(n,ncol(loading95))){
  PCi.ind = order(abs(loading95[,i]), decreasing=T)
  PCi.vals[,i] = loading95[PCi.ind,i]
  PCi.names[,i] = varnames95[PCi.ind]
}

df95 = data.frame(cbind(PC1vals = as.numeric(PCi.vals[1:10,1]),
                        PC1names = PCi.names[1:10,1],
                        PC2vals = PCi.vals[1:10,2],
                        PC2names = PCi.names[1:10,2],
                        PC3vals = PCi.vals[1:10,3],
                        PC3names = PCi.names[1:10,3],
                        PC4vals = PCi.vals[1:10,4],
                        PC4names = PCi.names[1:10,4],
                        PC5vals = PCi.vals[1:10,5],
                        PC5names = PCi.names[1:10,5]))

write.csv(df95, 'top10for95.csv')

## application on 95 amine data: binary classification
# ta98y01 = ifelse(ta98y > 0, 1, 0)
# 
# set.seed(03212016)
# pred = rep(0 ,nrow(ta98X))
# for(i in 1:nrow(ta98X)){
#   imod = rsdr(X=ta98X[-i,], y=ta98y01[-i], d=6, nslices=5)
#   
#   pred[i] = predict.rsir(imod, t(as.matrix(ta98X[i,])))
# }
# 
# plot.roc(ta98y01, pred)
# pred01 = ifelse(pred > mean(ta98y01), 1, 0)
# sum(pred01 == ta98y01) / 95
# sum((pred01==1 & ta98y01==1)) / sum(ta98y01==1)
# sum((pred01==0 & ta98y01==0)) / sum(ta98y01==0)

