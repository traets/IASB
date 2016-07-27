
#' D-error
#'
#' Function to calculate d error given a design, and parameter values.
#' @param par Vector containing parameter values.
#' @param design Design matrix.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @return D-error.
#' @export
d_err<-function (par, design, n_alts){

  des<-as.matrix(design)

  n_sets<-nrow(des)/n_alts
  group<-rep(seq(1, nrow(des)/n_alts, 1), each=n_alts)

  #probability
  p<-des %*% diag(par)
  p<-.rowSums(p,m= n_sets*n_alts,n=length(par))
  p<-exp(p)/rep(rowsum(exp(p),group), each=n_alts)

  infoM_des<- crossprod(des*p, des) - crossprod(rowsum(des*p, group))
  derr<- det(infoM_des)^(-1/length(par))

  return (derr)
}


#' Fisher Information of choice set
#'
#' Return the fisher information matrix of a particular choice set,
#' given parameter values.
#' @param par A vector containing the parameter values
#' @param c_set A matrix representing the choice set. Each row is a profile.
#' @return Fisher Information matrix
info_set<-function (par, c_set){

  p<-c_set %*% diag(par)
  p<-.rowSums(p, m = nrow(c_set), n = length(par))
  p<-exp(p)/sum(exp(p))

  i_set<-t(c_set) %*% (diag(p) - p %*% t(p)) %*% c_set

  return (i_set)
}

#' Fisher Information of design
#'
#' Returns the Fisher Information of a design, given parameter values
#' @param par A vector containing the parameter values
#' @param design A design matrix. Each row is a profile, successive rows form choice sets.
#' @return Fisher Information matrix.
info_design<-function (par, design, n_alts){

  des<-as.matrix(design)
  group<-rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)

  p<- des %*% diag(par)
  p<-.rowSums(p, m= nrow(des), n=length(par))
  p<-exp(p)/rep(rowsum(exp(p), group), each=n_alts)

  info_design<- crossprod(des*p, des) - crossprod(rowsum( des*p, group))
}
