## estimation ###

rm(list=ls())

#setwd
setwd("C:/Users/u0105757/Desktop/sequential_design.sas/workdirec")

#library
library(xlsx)
library(mlogit)
library(dplyr)
library(nnet)
library(rootSolve)
library(maxLik)
library(mnlogit)

#load design 
design<-read.xlsx("end.design3.1.xlsx",1)

#global 
n_alts<-length(unique(design$alt))
beta_true<-c(0.5, 0.5, 0.8, -0.5, 0.9, 0, 0, 1.2) 
beta_start1<-c(0, 0, 0, 0, 0, 0, 0, 0) #all zero
beta_start2<-c(0.8, 0.2, 1.1, 0.0, -0.7, 0, 0.1, -1.0) #random

#Functions:----------------------------------------------
#generate repsonses (Y)
respond<-function (beta, des) {
  
  p<-t(t(des) * beta)
  p<-.rowSums(p,m= nrow(des),n=length(beta))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,rep(seq(1,nrow(des)/n_alts,1), each = n_alts)), each=n_alts)
  
  Y<- as.numeric(p>0.5)
  
  return(Y)
}

#log likelihood
logLik<-function(beta, des, Y){
  
  p<-t(t(des) * beta)
  p<-.rowSums(p,m= nrow(des),n=length(beta))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,rep(seq(1,nrow(des)/n_alts,1), each = n_alts)), each=n_alts)
  
  logL<-sum(Y*log(p))
  
  return(logL)
  
}

#-------------------------------------------------------

###generate responses
des<-select(design,X1:X8)
Y<-respond(beta = beta_true, des)

######## ESTIMATION OF BETAS ###################

### mlogit {mlogit}
data<-as.data.frame(cbind(design, Y))
mlogit.design<-mlogit.data(data, choice = "Y", shape = "long", alt.var = "alt")

mlogit(Y ~ X1 + X2 + X3 +X4 + X5 + X6 + X7 + X8 -1, data = mlogit.design, method= "nr")
mlogit(Y ~ X1 + X2 + X3 +X4 + X5 + X6 + X7 + X8 -1, data = mlogit.design, start = beta_start1)
mlogit(Y ~ X1 + X2 + X3 +X4 + X5 + X6 + X7 + X8 -1, data = mlogit.design, start = beta_start2)

#mnlogit {mnlogit}
mnlogit(Y ~ X1 + X2 + X3 +X4 + X5 + X6 + X7 + X8 -1, data = mlogit.design)

#multinomial {nnet}
multinom(Y~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 -1, data = mlogit.design)
multinom(Y~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 -1, data = mlogit.design,start = beta_start1)
multinom(Y~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 -1, data = mlogit.design,start = beta_start2)

#optim {stats}
optim(par=beta_start1, logLik, gr=NULL, des=des, Y=Y)
optim(par=beta_start2, logLik, gr=NULL, des=des, Y=Y)
optim(par=beta_true, logLik, gr=NULL, des=des, Y=Y)

optim(par=beta_start1, logLik, gr=NULL, des=des, Y=Y, method = "BFGS") #quasi newton
optim(par=beta_start1, logLik, gr=NULL, des=des, Y=Y, method = "CG")
optim(par=beta_start1, logLik, gr=NULL, des=des, Y=Y, method = "L-BFGS-B", lower= -2, upper= 2)

#nlm {stats} 
nlm(logLik, p= beta_start1, des, Y=Y)

#maxLik
maxLik(logLik, start= beta_start1, des=des, Y=Y)
maxLik(logLik, start= beta_start2, des=des, Y=Y)
maxLik(logLik, start= beta_true, des=des, Y=Y)

#maxNR {maxLik}
maxNR(logLik, start= beta_start1, des=des, Y=Y)
maxNR(logLik, start= beta_start2, des=des, Y=Y)
maxNR(logLik, start= beta_true, des=des, Y=Y)


###compare responses estimated and true beta
beta_est<-c(  9.978675, 7.695254,  8.799287, -1.974818, 17.412430, -9.030795,  4.966539, 13.110765)
beta_est2<-c( 7.6559,  6.4852,   5.3865,  -2.1931,  15.1105,  -7.0944,  -1.1323,  13.3017 )
respond(beta = beta_est, des = des)
Y

logLik(beta = beta_est, des, Y) 
logLik(beta = beta_est2, des, Y) 
logLik(beta = beta_true, des, Y)







