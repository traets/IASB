rm(list=ls())

###compute the mode and the hessian matrix H of the posterior distribution (newton raphson) 
library(xlsx)
library(dplyr)
library(mlogit)
library(mvtnorm)
library(optimbase)

#load a design 
setwd("C:/Users/u0105757/Desktop/SHINY.7")
design<-read.xlsx("end.design3.1.xlsx",1)
design<-select(design,X1:X8 )
design.resp<-cbind(design, Y)

#set
n_alts<-2
beta_true<-c(0.5, 0.5, 0.8, -0.5, 0.9, 0, 0, 1.2)
beta_prior<-c(0, 0, 0, 0, 0, 0, 0, 0)
var<-c(1,1,1,1,1,1,1,1)
covar<-diag(var)
df<-8

#FUNCTIONS#########################################################################
#loglikelihood function
logLik<-function(beta, des, n_alts, Y){
  
  p<-t(t(des) * beta)
  p<-.rowSums(p,m= nrow(des),n=length(beta))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,rep(seq(1,nrow(des)/n_alts,1), each = n_alts)), each=n_alts)

  logL<-sum(Y*log(p))
  
  return(logL)
  
}

beta<-beta_true
design

#likelihood function
Lik<-function(beta, des, n_alts, Y){
  
  p<-t(t(des) * beta)
  p<-.rowSums(p,m= nrow(des),n=length(beta))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,rep(seq(1,nrow(des)/n_alts,1), each = n_alts)), each=n_alts)
  
  L<-prod(p[Y])
  
  return(L)
  
}

#Function to generate repsonses (Y)
respond<-function (beta, des) {
  p<-t(t(des) * beta)
  p<-.rowSums(p,m= nrow(des),n=length(beta))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,rep(seq(1,nrow(des)/n_alts,1), each = n_alts)), each=n_alts)
  
  Y<- as.numeric(p>0.5)
  
  return(Y)
}

#Function to generate lattice points
lattice.points <- function(K, b = 2, m =9){
  base <- function(num){
    a1 <- c1 <- rep(0, m)
    a1[m] <- c1[m] <- num %% b
    for(j in seq(m - 1, 1, by = -1)){
      tem <- (num - c1[j + 1]) / (b^(m - j))
      a1[j] <- tem %% b
      c1[j] <- c1[j + 1] + a1[j] * (b^(m - j))
    }
    return(a1)
  }
  
  a <- 1571
  N <- b^m # number of lattice points
  u <- runif(K)
  av <- rep(1, K)
  
  for(i in 2:K) {av[i] <- (a * av[i - 1]) %% N}
  
  e <- matrix(0, N, K)
  seq <- rep(0, m)
  kk <- -m
  
  for(k in 1:m) {seq[k] <- b^kk; kk <- kk + 1;}
  
  for(i in 1:N){
    ei <- crossprod(seq, base(i - 1)) * av + u
    e[i, ] <- ei - floor(ei)
  }
  
  latt <- matrix(0, N, K)
  for(i in 1:K){
    for(j in 1:N){latt[j, i] <- 1 - abs(2 * e[j, i] - 1)}
  }
  
  latt <- qnorm(latt)
  
  return(latt)
}

#Function to compute hessian
hessian<-function (beta, des, covar){
  
  des<-as.matrix(des)
  p<- des %*% diag(beta)
  p<-.rowSums(p,m= nrow(des),n=length(beta))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,rep(seq(1,nrow(des)/n_alts,1), each = n_alts)), each=n_alts)
  
  info<-crossprod(des*p,des)-crossprod(rowsum(des*p,rep(seq(1,nrow(des)/n_alts,1),each=n_alts)  ))
  
  hess<-(-info-solve(covar))
  
  return(hess)
  
}

#Function to shift lattice points to mode posterior (draws from importance density)
grid.MVT<- function (mode , cvar, lattice){
  
  mean<-t(beta_est)
  parm<-length(beta_est)
  X<-matrix(NA,nrow(grid),parm)
  A<-chol(cov)
  
  
  for(i in 1:nrow(grid))
  {
    
    inv<-rgamma(1,df/2)
    invw<-inv/(df/2)
    W<-1/invw
    
    Z<-ranz[i,]  
    
    r<-mean + sqrt(W) * (Z %*% t(A))
    X[i, ]<-t(r)
  }
  
  return (X)
}

#Function to calculate density of draw of importance density 
dens.MVT<-function (beta_sample, beta_est, cvar){
  
  #maybe do t() first 
  n<-nrow(beta_sample)
  dif<-beta_est-tbeta_sample
  invcov<-solve(cvar)  
  differ<-as.numeric(t(dif)%*% invcov %*% dif)
  
  iMVSTd=1/(det(cvar)^(0.5))*(1+((1/df)*differ))^(-(df+length(beta_sample))/2)
}
###################################################################################

Y<-respond(beta_true, des=design )

### (1) Compute the mode B and the hessian matrix 
# mode B
mlogit_des<-as.data.frame(cbind(design, Y, alt= rep(seq(1:n_alts), nrow(design)/n_alts)))
mlogit_design<-mlogit.data(mlogit_des, choice = "Y", shape = "long", alt.var = "alt")
est<-mlogit(Y ~ X1 + X2 + X3 +X4 + X5 + X6 + X7 + X8 -1, data = mlogit_design, method= "nr")
beta_est<-as.numeric(est$coefficients)

#hessian
H<-hessian(des=design, beta = beta_est, covar = covar)
covar_2<--solve(H)
covar_2

### (2) Take R draws from the importance density 
#imp dens = multivar t with mean= B, covar= -HESS^-1, v degrees of freedom (8)

#draws
g_draws<-grid.MVT(beta_est, covar_2, grid)

### BAYESIAN 
den_prior<-numeric(ncol(g_draws))
LK<-numeric(ncol(g_draws))

for (r in 1:ncol(g_draws)){

#prior
prior1=1/((2*3.1415926)^(length(beta_est/2)*det(covar)^(0.5)))
den_prior[r]<-prior1*exp(-0.5* (g_draws[,r]-beta_prior) %*% solve(covar) %*% transpose(g_draws[,r]-beta_prior))

#likelihood
LK[r]<-Lik(beta =g_draws[,r], des=design, n_alts,  Y)


}




















