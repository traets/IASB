rm(list=ls())

#libraries
library(xlsx)
library(dplyr)
library(maxLik)

#load a design
#setwd("C:/Users/u0105757/Desktop/SHINY.7")
setwd("C:/Users/u0105757/Desktop/sequential design/R.code/updating")
#design<-read.xlsx("end.design3.1.xlsx",1)
design<-read.xlsx("init_des.xlsx",1, header=F)

#set
n_alts<-2
beta_true<-c(-0.46,-0.21, 0.19, 0.22)
beta_true
beta_prior<-c(-0.5,0,-0.5,0)
var_prior<-c(1,1,1,1)
covar_prior<-diag(var_prior)
inv_covar_prior<-solve(covar_prior)
beta_start<-rep(0,length(beta_prior))
df<-8
logprior1<-df/2*log(2*pi)-0.5*log(det(covar_prior))


#FUNCTIONS#########################################################################
#generate start design

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
logPost<-function(b, des, Y){

  b<-as.vector(b)

  p<-t(t(des) * b)
  p<-.rowSums(p,m= nrow(des),n=length(b))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,rep(seq(1,nrow(des)/n_alts,1), each = n_alts)), each=n_alts)

  logL<-sum(Y*log(p))
  logprior2=-0.5*(b -t(beta_prior))%*%solve(covar_prior)%*%(as.matrix(b)-beta_prior)

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
grid.MVT<- function (mode, cvar, lattice){

  mean<-t(mode)
  parm<-length(mode)
  X<-matrix(NA,nrow(lattice),parm)
  A<-chol(cvar)
  t(chol(cvar))

  for(i in 1:nrow(lattice))
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
Y<-c(0,1,1,0,0,1,0,1,0,1)

### (1) Compute the mode B and the hessian matrix
beta_est<-maxNR(logPost, start= beta_start, des=design, Y=Y)$estimate

#constant in calculation dens of MV normal prior
prior1<-(2*pi)^(-length(beta_est)/2)*(det(covar_prior))^(-0.5)

#hessian
H<-hessian(des=design, beta = beta_est, covar = covar_prior)
covar_mvt<--solve(H)

### (2) Take R draws from the importance density
#imp dens = multivar t with mean= B, covar= -HESS^-1, v degrees of freedom (8)

#draws
grid<-lattice(length(beta_est))
#grid<-read.sas7bdat("C:/Users/u0105757/Desktop/sequential design/R.code/updating/latfull.sas7bdat", debug=FALSE)
grid<-as.matrix(grid)
g_draws<-grid.MVT(beta_est, covar_mvt, grid)


### IMPORTANCE SAMPLING ##################################################################################
prior<-numeric(nrow(g_draws))
LK<-numeric(nrow(g_draws))
dens_g<-numeric(nrow(g_draws))
weights<-numeric(nrow(g_draws))

for (r in 1:nrow(g_draws)){

  #prior
  prior[r]<-prior1*exp(-0.5* (g_draws[r,]-beta_prior) %*% solve(covar_prior) %*% as.matrix(g_draws[r,]-beta_prior))

  #likelihood
  LK[r]<-Lik(beta =g_draws[r,], des=design,  Y)

  #density of g
  dens_g[r]<-dens.MVT(beta_sample = g_draws[r,], beta_est = beta_est, cvar = covar_mvt)

}



### (3) compute the weights of samples
w<-LK*prior/dens_g
w<-w/sum(w)


#plot posterior
h<-hist(g_draws[,1], breaks=20)
group_fac<-factor(rep(1:(length(h$breaks)-1), h$counts))

weights_groups<-split(w, group_fac)
mean_weight<-unlist(lapply(weights_groups, mean))

dens<-mean_weight*h$counts

plot(h$mids, dens, col= 'red',type = 'l')

### generate new choice sets ###
# this program takes as input an initial design, beta samples and weights for the samples.
# as output it produces an efficient choice set.

#FUNCTIONS

#all profile combinations
full.profiles<-function (levels_attributes,contrasts) {

  #as list
  levels_attributes<-as.list(levels_attributes)

  #Make list with all atribute levels
  fun<-function (x){return(1:x)}
  levels_attributes<-lapply(X = levels_attributes, fun)

  #generate every combination
  full_design<-as.data.frame(expand.grid(levels_attributes))

  #create factors:
  col_names <- names(full_design)
  full_design[,col_names] <- lapply(full_design[,col_names] , factor)

  #contrast list
  con <- list()
  for(i in 1:length(levels_attributes)){
    name <- paste('Var',i,sep='')
    con[name] <- contrasts
  }

  #coding
  full_D<-as.data.frame(model.matrix(~ ., full_design, contrasts = con))
  #delete intercept
  full_D<-full_D[,-1]

  return(full_D)

}

#all possible choice sets
full.choice.sets<- function(full_design, n_alts){


  fun<-function (x){return(1:x)}
  N.row.design<-nrow(full_design)

  sets.1<-rep(N.row.design,n_alts)
  sets.2<-as.list(sets.1)
  sets.3<-lapply(X = sets.2, fun)
  set.design<-as.data.frame(expand.grid(sets.3))

  #to global env
  assign("set.design", set.design, envir=globalenv())

  #return
  return (set.design)
}

#information of choice set
info.set<-function (beta, c_set, n_alts){

  p<-c_set %*% diag(beta)
  p<-.rowSums(p,m= n_alts, n=length(beta))
  p<-exp(p)/sum(exp(p))

  info_set<-t(c_set)%*% (diag(p)-p%*%t(p)) %*% c_set

  return (info_set)
}

#information of design
info.design<-function (beta, des, n_alts){

  des<-as.matrix(des)
  group<-rep(seq(1,nrow(des)/n_alts,1),each=n_alts)

  p<- des %*% diag(beta)
  p<-.rowSums(p,m= nrow(des),n=length(beta))
  p<-exp(p)/rep(rowsum(exp(p),group), each=2)

  info_design<-crossprod(des*p, des)-crossprod(rowsum( des*p,group))
}

################################################################

#all possible profiles
candit_set<-full.profiles(c(3,3), 'contr.sum')

#all possible combinations choice sets
full_combinations<-full.choice.sets(candit_set,2)

#start DB-error
DB_best<-1000

#for each potential set:
for (p in 1:nrow(candit_set)){

  set<-as.matrix(candit_set[as.numeric(full_combinations[8,]), ])

  #for each beta draw calculate d error:
    DB_error<-0

    for ( s in 1: nrow(g_draws)){

      info_set<-info.set(g_draws[s,], set, n_alts)
      info_des<-info.design(g_draws[s,], design, n_alts)

      DB_error=DB_error+(det(info_des+info_set+inv_covar_prior)*(-1/ncol(g_draws))*w[s])

    }

    #select new set if smaller Db error
    if (DB_error< DB_best){ final_set<-set; DB_best<-DB_error}
}



final_set


#info init design

#info






