## analyze outputs of 95 and 508 compound analysis
rm(list=ls())
setwd("D:/Study/My projects/Statchem-Diudea/Codes")

## 508 compounds
load("err508.Rda")
mean.full = round(apply(err.mat1, 2, mean),2)
sd.full = round(apply(err.mat1, 2, sd),3)
rm(err.mat1)

load("err508basak-final.Rda")
mean.b = round(apply(err.mat508, 2, mean),2)
sd.b = round(apply(err.mat508, 2, sd),3)
rm(err.mat508)

load("err508diudea-final.Rda")
mean.d = round(apply(err.mat508, 2, mean),2)
sd.d = round(apply(err.mat508, 2, sd),3)
rm(err.mat508)

z = data.frame(cbind(mean.full, sd.full, mean.b, sd.b, mean.d, sd.d),
               row.names = c("PCR","PLS","LASSO","SCAD","RF","GBM"))
z

## 95 compounds
load("err95.Rda")
mean.full = round(apply(err.mat, 2, median),2)
sd.full = round(apply(err.mat, 2, mad),3)
rm(err.mat)

load("err95basak.Rda")
mean.b = round(apply(err.mat, 2, mean),2)
sd.b = round(apply(err.mat, 2, sd),3)
rm(err.mat)

load("err95diudea.Rda")
mean.d = round(apply(err.mat, 2, median),2)
sd.d = round(apply(err.mat, 2, mad),3)
rm(err.mat)

z = data.frame(cbind(mean.full, sd.full, mean.b, sd.b, mean.d, sd.d),
               row.names = c("PCR","PLS","LASSO","SCAD","RF","GBM"))
z
