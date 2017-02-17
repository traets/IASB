rm(list=ls())

#libraries
library(xlsx)
library(dplyr)
library(mvtnorm)
library(maxLik)
library(sas7bdat)

#load a design 
#setwd("C:/Users/u0105757/Desktop/SHINY.7")
setwd("C:/Users/u0105757/Desktop/sequential design/R.code/updating")
#design<-read.xlsx("end.design3.1.xlsx",1)
design<-read.xlsx("init_des.xlsx",1, header=F)


#set
n_alts<-2
beta_true<-read.sas7bdat("C:/Users/u0105757/Desktop/sequential design/commented sas/IASB with comments/betatrue.sas7bdat")$COL1
beta_prior<-c(-0.5,0,-0.5,0)
var_prior<-c(1,1,1,1)
covar_prior<-diag(var_prior)
beta_start<-rep(0,length(beta_prior))
df<-8
logprior1<-df/2*log(2*3.1415926)-0.5*log(det(covar_prior))


#FUNCTIONS#########################################################################
#loglikelihood function
logLik<-function(beta, des, Y){
  
  p<-t(t(des) * beta)
  p<-.rowSums(p,m= nrow(des),n=length(beta))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,rep(seq(1,nrow(des)/n_alts,1), each = n_alts)), each=n_alts)
  
  logL<-sum(Y*log(p))
  
  return(logL)
  
}

#likelihood function
Lik<-function(beta, des, Y){
  
  p<-t(t(des) * beta)
  p<-.rowSums(p,m= nrow(des),n=length(beta))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,rep(seq(1,nrow(des)/n_alts,1), each = n_alts)), each=n_alts)
  
  L<-prod(p^Y)
  
  return(L)
  
}

#log Posterior
logPost<-function(beta, des, Y){
  
  p<-t(t(des) * beta)
  p<-.rowSums(p,m= nrow(des),n=length(beta))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,rep(seq(1,nrow(des)/n_alts,1), each = n_alts)), each=n_alts)
  
  logL<-sum(Y*log(p))
  logprior2=-0.5*(beta-t(beta_prior))%*%solve(covar_prior)%*%(t(beta)-beta_prior)
  
  logpost<-logprior1 + logprior2 + logL
  
  return(logpost)
  
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
lattice <- function(K, b = 2, m =9){
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
  A<-chol(cvar)
  
  
  for(i in 1:nrow(grid))
  {
    
    inv<-rgamma(1,df/2)
    invw<-inv/(df/2)
    W<-1/invw
    
    Z<-lattice[i,]  
    
    r<-mean + sqrt(W) * (Z %*% t(A))
    X[i, ]<-t(r)
  }
  
  return (X)
}

#Function to calculate density of draw of importance density 
dens.MVT<-function (beta_sample, beta_est, cvar){
  
  n<-length(beta_sample)
  dif<-beta_est-beta_sample
  invcov<-solve(cvar)  
  differ<-as.numeric(t(dif)%*% invcov %*% dif)
  
  iMVSTd=1/(det(cvar)^(0.5))*(1+((1/df)*differ))^(-(df+length(beta_sample))/2)
}

###################################################################################

#Create responses
#Y<-respond(beta_true, des=design )
Y<-c(0,1,0,1,0,1,0,1,1,0)

### (1) Compute the mode B and the hessian matrix 
beta_est<-maxNR(logPost, start= beta_start, des=design, Y=Y)$estimate

#constant in calculation dens of MV normal prior
prior1<-(2*3.1415926)^(-length(beta_est)/2)*(det(covar_prior))^(-0.5) 


#hessian
H<-hessian(des=design, beta = beta_est, covar = covar_prior)
covar_mvt<--solve(H)

### (2) Take R draws from the importance density 
#imp dens = multivar t with mean= B, covar= -HESS^-1, v degrees of freedom (8)

#draws
#grid<-lattice(length(beta_est))
grid<-read.sas7bdat("C:/Users/u0105757/Desktop/sequential design/R.code/updating/latfull.sas7bdat", debug=FALSE)
grid<-grid[1:16,]
g_draws<-grid.MVT(beta_est, covar_mvt, grid)


### BAYESIAN ##################################################################################???
dens_prior<-numeric(nrow(g_draws))
LK<-numeric(nrow(g_draws))
dens_g<-numeric(nrow(g_draws))
weights<-numeric(nrow(g_draws))

for (r in 1:nrow(g_draws)){
  
  #prior
  dens_prior[r]<-prior1*exp(-0.5* (g_draws[r,]-beta_prior) %*% solve(covar_prior) %*% transpose(g_draws[r,]-beta_prior))
  
  #likelihood
  LK[r]<-Lik(beta =g_draws[r,], des=design,  Y)
  
  #density of g 
  dens_g[r]<-dens.MVT(beta_sample = g_draws[r,], beta_est = beta_est, cvar = covar_post)
  
}



### (3) compute the weights of samples
weights<-LK*dens_prior/dens_g 
weights

grid_post<-weights/sum(weights)

plot(g_draws[,1], grid_post)

post_dens<-data.frame(cbind(g_draws, posterior))

### (4) 

