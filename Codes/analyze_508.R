setwd("C:/Study/My projects/StatChem_Robust")
rm(list=ls());
source('misc_functions.R')
library(rrcov)
library(fda.usc)
library(ldr)

# predictor functions
get.preds2 = function(X,y,varpct=90){
  set.seed(03212016)
  pred1 = rep(0 ,nrow(X))
  pred2 = pred1
  for(i in 1:nrow(X)){
    imod = sdr(X=X[-i,], y=y[-i], varpct=varpct, nslices=5)
    pred1[i] = predict.rsir(imod, t(as.matrix(X[i,])))
    
    imod = rsdr(X=X[-i,], y=y[-i], varpct=varpct, nslices=5)
    pred2[i] = predict.rsir(imod, t(as.matrix(X[i,])))
  }
  
  #plot.roc(data508y, pred)
  pred101 = ifelse(pred1 > mean(y), 1, 0)
  pred201 = ifelse(pred2 > mean(y), 1, 0)
  100*c(sum(pred101 == y) / 508,
    sum((pred101==0 & y==0)) / sum(y==0),
    sum((pred101==1 & y==1)) / sum(y==1),
    sum(pred201 == y) / 508,
    sum((pred201==0 & y==0)) / sum(y==0),
    sum((pred201==1 & y==1)) / sum(y==1))
}

## application on 95 amine data
data508 = read.table('data508.txt', header=T)
#data508atom = read.table('data508atom.txt')

vartype = as.numeric(data508[1,2:308])
data508X1 = as.matrix(data508[2:509,2:308])
data508X1 = logt(data508X1)
# data508X2 = as.matrix(data508atom[-1,-2518])
# data508X = cbind(data508X1, data508X2)
data508X1 = log(-floor(min(data508X1)) + data508X1)
#data508X1 = scale(data508X1)
data508y = as.numeric(data508[2:509,309])

outmat = rbind(get.preds2(data508X1[,which(vartype==1)], data508y),
               get.preds2(data508X1[,which(vartype==1 | vartype==2)], data508y),
               get.preds2(data508X1[,which(vartype==1 | vartype==2 | vartype==3)], data508y),
               get.preds2(data508X1, data508y),
               get.preds2(data508X1[,which(vartype==1)], data508y),
               get.preds2(data508X1[,which(vartype==2)], data508y),
               get.preds2(data508X1[,which(vartype==3)], data508y))
round(outmat, digits=2)
               
# principal component regression
get.preds.logistic = function(X,y,varpct=90){
  set.seed(03212016)
  pred3 = rep(0 ,nrow(X))
  for(i in 1:nrow(X)){
    iPC = PcaClassic(X[-i,])
    eigen.cumsum = cumsum(iPC@eigenvalues)
    d = min(which(eigen.cumsum >= (varpct/100) * eigen.cumsum[length(eigen.cumsum)]))
    iS = iPC@scores[,1:d]
    imod = glm.fit(x=as.matrix(iS), y=y[-i])
    
    X.beta = sum((X[i,] %*% iPC@loadings)[1:d] * imod$coef)
    pred3[i] = exp(X.beta) / (1+exp(X.beta))
  }
  
  pred301 = ifelse(pred3 > mean(y), 1, 0)
  100*c(sum(pred301 == y) / 508,
        sum((pred301==0 & y==0)) / sum(y==0),
        sum((pred301==1 & y==1)) / sum(y==1))
}

outmat = rbind(get.preds.logistic(data508X1[,which(vartype==1)], data508y),
               get.preds.logistic(data508X1[,which(vartype==1 | vartype==2)], data508y),
               get.preds.logistic(data508X1[,which(vartype==1 | vartype==2 | vartype==3)], data508y),
               get.preds.logistic(data508X1, data508y),
               get.preds.logistic(data508X1[,which(vartype==1)], data508y),
               get.preds.logistic(data508X1[,which(vartype==2)], data508y),
               get.preds.logistic(data508X1[,which(vartype==3)], data508y))
round(outmat, digits=2)

## 508 PCA
# convert to rank
delta = 1e-3
data = data508X1
y = data
n = nrow(y)
p = ncol(y)
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
svd508 = svd(y)
loading508 = svd508$v

PCi.vals = matrix(0, nrow=p, ncol=min(n,p))
PCi.names = matrix("0", nrow=p, ncol=min(n,p))
varnames508 = names(read.delim('varnames508.txt'))[-1]

for(i in 1:min(n,p)){
  PCi.ind = order(abs(loading508[,i]), decreasing=T)
  PCi.vals[,i] = loading508[PCi.ind,i]
  PCi.names[,i] = varnames508[PCi.ind]
}

df508 = data.frame(cbind(PC1vals = as.numeric(PCi.vals[1:10,1]),
                        PC1names = PCi.names[1:10,1],
                        PC2vals = PCi.vals[1:10,2],
                        PC2names = PCi.names[1:10,2],
                        PC3vals = PCi.vals[1:10,3],
                        PC3names = PCi.names[1:10,3],
                        PC4vals = PCi.vals[1:10,4],
                        PC4names = PCi.names[1:10,4],
                        PC5vals = PCi.vals[1:10,5],
                        PC5names = PCi.names[1:10,5]))

write.csv(df508, 'top10for508.csv')
