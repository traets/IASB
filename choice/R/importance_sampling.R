
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
#' @export
imp_sampling <- function (prior_mode, prior_covar, design,  n_alts, Y){

  #cte
  prior1<-(2*pi)^(-length(mode)/2)*(det(prior_covar))^(-0.5)
  lp1<-length(prior_mode)/2*log(2*pi)-0.5*log(det(prior_covar))

  #log Posterior
  logPost<-function(par, design, Y, n_alts, logprior1=lp1){

    des<-design
    p<-t(t(des) * par)
    p<-.rowSums(p, m= nrow(des), n=length(par))
    expp<-exp(p)
    p<-expp/rep(rowsum(expp,rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)

    log_L<-sum(Y*log(p))

    logprior2=-0.5*(par -t(prior_mode))%*%solve(prior_covar)%*%(as.matrix(par)-prior_mode)

    logpost<-logprior1 + logprior2 + log_L

    return(logpost)

  }

  #Hessian
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

  #estimate importance mode
  imp_mode<-maxLik::maxNR(logPost, start= prior_mode, design=des, Y=Y, n_alts=n_alts)$estimate

  #draws from importance density
  H<-hessian(par = imp_mode, design = des, covar = prior_covar, n_alts = 2)
  g_covar<--solve(H)
  g_draws<-lattice_mvt(mode=imp_mode, cvar = g_covar, df=length(imp_mode))

  #likelihood function
  Lik<-function(par, design, n_alts, Y){

    des<-as.matrix(design)

    p<-t(t(des) * par)
    p<-.rowSums(p, m= nrow(des), n=length(par))
    expp<-exp(p)
    p<-expp/rep(rowsum(expp, rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)

    L<-prod(p^Y)

    return(L)

  }

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




