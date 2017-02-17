rm(list=ls())


# negative KL values ?
# adjust full_sets function in rest of shiny scripts


#libraries
library(aprof)
library(choice)

#initials needed
n_alts=2
levels<-c(3,2,3,2,2,2)
p_mod<- c( 1, 2, 1, -3, 1, 5, 1, -3)
var<-1
p_cov<- diag(rep(var,length(p_mod)))
des<-matrix(data = c(-1, -1,	1, 1, 0, -1, 1, -1, 0, 1, -1,	-1,	-1,	1, -1, 1), nrow=2	)
y_bin<-c(0,1)

#samples
samples<-imp_sampling(prior_mode = p_mod, prior_covar = p_cov, design = des, n_alts = n_alts, Y=y_bin)
sam<-samples[[1]]
w<-samples[[2]]

#list
 fcomb <- full_sets(lvls = levels, n_alts=n_alts, mindiff = 4)
 fcomb <-as.matrix(fcomb)
fp<-profiles(lvls=levels)

 #####################################################################################################


KL_original<-function (lvls, par_samples, weights, n_alts) #34 sec
{

  minKL <- 0
  for (s in 1:nrow(fcomb)) {

    set <- as.matrix(fp[as.numeric(fcomb[s, ]), ])
    num <- set %*% t(par_samples)

    mmat <- t(apply(num, 2, max))
    nummax <- exp(sweep(num, MARGIN = 2, mmat, FUN = "-"))
    denom <- colSums(nummax)
    probs <- sweep(nummax, MARGIN = 2, denom, FUN = "/")

    # expnum<-exp(num)
    # total<-colSums(expnum)
    # probs<-sweep(expnum, MARGIN = 2, total, FUN = "/")
    #

    logprobs <- log(probs)
    wprob <- sweep(probs, MARGIN = 2, weights, FUN = "*")
    totwprob <- rowSums(wprob)
    logwprob <- sweep(logprobs, MARGIN = 2, weights, FUN = "*")
    totlogwprob <- rowSums(logwprob)
    klinfo <- 0
    for (a in 1:n_alts) {
      klinfo <- klinfo + (totwprob[a] * (log(totwprob[a]) -
                                           totlogwprob[a]))
    }

    if (klinfo > minKL) {
      minKL <- klinfo
      new_set <- set
    }
  }
  return(new_set)
}


KL_1<-function (fcomb, fp, par_samples, weights, n_alts) #5.7 sec
{

  minKL <- 0
  new_set<-numeric()

  for (s in 1:nrow(fcomb)) {


    set <- fp[as.numeric(fcomb[s,]), ]
    num <- tcrossprod(par_samples, set)

    expnum<-exp(num)
    total<-rowSums(expnum)
    probs<-expnum / total

    totwprob <- t(probs) %*% w

    logprobs <- log(probs)
    totlogwprob <- t(logprobs) %*% w


    klinfo <- 0

    for (a in 1:n_alts) {
      klinfo <- klinfo + (totwprob[a] * (log(totwprob[a]) -
                                           totlogwprob[a]))
    }

    if (klinfo > minKL) {
      minKL <- klinfo
      new_set <- set
    }
  }
  return(new_set)
}


##### profiling ########
## save function to a source file and reload
dump("KL_1",file="KL_ori.R")
source("KL_ori.R")

## create file to save profiler output
tmp<-tempfile()

## Profile the function
Rprof(tmp,line.profiling=TRUE)
KL_1(lvls = levels, par_samples = sam, weights = w, n_alts = n_alts)
Rprof(append=FALSE)

## Create a aprof object
fooaprof<-aprof("KL_ori.R",tmp)
## display basic information, summarize and plot the object
fooaprof
summary(fooaprof)
plot(fooaprof)
profileplot(fooaprof)
######################



KL_1(fcomb = fcomb, fp= fp, par_samples = sam, weights = w)
lvls = levels; par_samples = sam; weights = w; n_alts = n_alts


 library(devtools)
 install_github("traets/choice/choice")





