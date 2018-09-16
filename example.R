#####  Gaussian Process with Errors in Variables (GPEV) #########

rm(list = ls())

library(fields)
library(mvtnorm)
library(matrixStats) # for logSumExp
library(stats)
library(FSA)
library(ggplot2)

source("appGPEV.R")

sig = 0.2       # regression error sd
ntrain = 100    # sample size
J = 40          # Truncated number
delta = 0.01    # measurement error variance 
niter = 1500    # Number of iterations 

# regression function 
myfun <- function(x) sin(pi*x/2)/(1+2*x^2*(sign(x)+1))

# generate data 
locs=runif(ntrain,min=-3,max=3)
locs.err = locs + rnorm(ntrain,0,sqrt(delta))
ztrue = myfun(locs)
z = ztrue + rnorm(ntrain,0,sig)
train_data = rbind(locs, locs.err, z, ztrue)

# test data 
ntest = 100   
xtest = seq(-3,3,length.out=ntest)
ytest = myfun(xtest) 
test_data = rbind(xtest,ytest)


## Running GPEV method 
r = appGPEV(train_data, test_data, delta, J, niter, theta = c(5,1,10))

# out of sample prediction on test grid
ypred = r[[3]]
# out of sample 95% point-wise credible intervals
CIpred = r[[7]]

data = data.frame(xtest, ypred)
p <- ggplot(data, aes(x=xtest, y=ytest), col = variable)+
     geom_line(aes(x=xtest, y=ytest), colour="red")+
     geom_line(aes(x=xtest, y=ypred), colour="black")+
     geom_ribbon(aes(ymin=CIpred[1,], ymax=CIpred[2,]), alpha=0.2)

plabs<- p+labs(x = "x", y = "f(x)")
plabs + theme( 
  axis.title.x = element_text(size=8, face="bold"),
  axis.title.y = element_text(size=8, face="bold"),
  legend.text = element_text(colour="blue", size=10, face="bold"))
     

