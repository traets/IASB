#' log Posterior
#'
#' Calculates the logposterior with a normal prior density
#' @param par parameter value(s) for wich the logposterior is calculated
#' @param design The Design matrix
#' @param Y the response vector
#' @param n_alts The number of alternatives in each choice set
#' @param prior_mode vector containing the prior mode
#' @param prior_covar matrix containing the prior covariance
#' @return the logposterior probability
logPost<-function(par, design, Y, n_alts){

  des<-design
  p<-t(t(des) * par)
  p<-.rowSums(p, m= nrow(des), n=length(par))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp, rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)

  log_L<-sum(Y*log(p))

  logprior2=-0.5*(par -t(prior_mode))%*%solve(prior_covar)%*%(as.matrix(par) - prior_mode)

  logpost<-(length(prior_mode)/2*log(2*pi)-0.5*log(det(prior_covar))) + logprior2 + log_L

  return(logpost)

}

#' Hessian
#'
#' @param par parameter value(s)
#' @param design The design matrix
#' @param covar the covariance matrix
#' @param n_alts The number of alternatives in each choice set
#' @return the hessian matrix
hessian<-function (par, design, covar, n_alts){

  des<-as.matrix(design)
  p<- des %*% diag(par)
  p<-.rowSums(p, m= nrow(des), n=length(par))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp, rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)

  info<-crossprod(des*p, des) - crossprod(rowsum(des*p, rep(seq(1, nrow(des)/n_alts, 1), each=n_alts)))

  hess<-(-info - solve(covar))

  return(hess)

}

#' Likelihood function
#'
#' @param par parameter value(s)
#' @param design The design matrix
#' @param n_alts The number of alternatives in each choice set
#' @param Y the response vector
#' @return the likelihood
Lik<-function(par, design, n_alts, Y){

  des<-as.matrix(design)

  p<-t(t(des) * par)
  p<-.rowSums(p, m= nrow(des), n=length(par))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp, rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)

  L<-prod(p^Y)

  return(L)

}

#' Importance sampling
#'
#' This functions samples from an imortance density (multivariate t-distribution),
#' and gives weightes to the samples according to the posterior distribution. The prior is
#' a normal distribution.
#' @param prior_mode Vector containing the mode of the multivariate normal distribution which is the prior.
#' @param prior_covar The covariance matrix of the multivariate normal prior distribution.
#' @param design A design matrix in which each row is a profile.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param Y A binary response vector.
#' @return A list containing samples and their associated weights, which together represent the posterior distribution.
#' @export
imp_sampling <- function (prior_mode, prior_covar, design,  n_alts, Y, ...){

  #cte
  prior1<-(2*pi)^(-length(mode)/2)*(det(prior_covar))^(-0.5)

  #estimate importance mode
  imp_mode<-maxLik::maxNR(logPost, start= prior_mode, design=des, Y=Y, n_alts=n_alts)$estimate

  #draws from importance density
  H<-hessian(par = imp_mode, design = des, covar = prior_covar, n_alts = n_alts)
  g_covar<--solve(H)
  g_draws<-lattice_mvt(mode=imp_mode, cvar = g_covar, df=length(imp_mode), ...)

  #importance density
  g_dens<-function (par, g_mode, g_covar, df=length(imp_mode)){

    n<-length(par)
    dif<-g_mode-par
    invcov<-solve(g_covar)
    differ<-as.numeric(t(dif)%*% invcov %*% dif)

    iMVSTd=1/(det(g_covar)^(0.5))*(1+((1/df)*differ))^(-(df+length(par))/2)
  }

  #vectors
  prior=LK=dens_g=weights<- numeric(nrow(g_draws))

  for (r in 1:nrow(g_draws)){

    #prior
    prior[r]<-prior1*exp(-0.5* (g_draws[r, ]-prior_mode) %*% solve(prior_covar) %*% as.matrix(g_draws[r, ]- prior_mode))
    #likelihood
    LK[r]<-Lik(par = g_draws[r, ], design=des, Y=Y, n_alts = n_alts)
    #density of g
    dens_g[r]<-g_dens(par = g_draws[r, ], g_mode = imp_mode, g_covar = g_covar)

  }

  #compute the weights of samples
  w<-LK*prior/dens_g
  w<-w/sum(w)

  return(list(g_draws,w))
}

#' modes
#'
#' This functions returns a vector with the posterior modes,
#' given samples from the posterior and their weights (importance sampling).
#' @param samples Matrix with samples coming from the posterior distribution. Each row is a sample.
#' @param weights Vector containing the weights of the samples.
#' @return The mode of the posterior distribution.
#' @export
modes<-function(samples, weights, s){

  mode<-numeric(ncol(samples))

  for (i in 1:ncol(samples)){

    ### posterior probability
    #mean weight intervals
    h<-hist(sam[ ,i], breaks = nrow(samples)/s)
    group_fac<-factor(rep(1:(length(h$breaks)-1), h$counts))
    w_groups<-split(w, group_fac)
    mean_weight<-unlist(lapply(w_groups, mean))

    #posterior= mean weights * frequency
    freq<-h$counts[h$counts != 0]
    int_mids<-h$mids[h$counts != 0]
    pos<-freq*mean_weight/sum(freq*mean_weight)
    d<-as.data.frame(cbind(int_mids,pos))

    d<-d[with(d, order(-pos)), ]
    mode[i]<-d$int_mids[1]

  }

  return(mode)

}




