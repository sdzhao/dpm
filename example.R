## example.R

rm(list=ls());
dyn.load("dpm.so");
source("dpm.R");
library(MASS);

## **************************************************************
## generate data
## **************************************************************
n0 <- 100;
n1 <- 100;
p <- 10;

O1 <- matrix(0.5,nrow=p,ncol=p)+diag(0.5,p);
D0 <- diag(c(1,1,rep(0,p-2)));
O0 <- O1+D0;

X1.t <- mvrnorm(n1,rep(0,p),solve(O1));
X0.t <- mvrnorm(n0,rep(0,p),solve(O0));

## **************************************************************
## calculate difference in precision matrices
## **************************************************************
fit.aic <- dpm(X1.t,X0.t,nlambda=10,tuning="aic");
fit.cv <- dpm(X1.t,X0.t,nlambda=10,tuning="cv",folds=3);
