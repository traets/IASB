
#' Profile generation
#'
#' Function to generate all possible combinations of attribute levels.
#' @param lvls  A vector which contains for each attribute, the number of levels.
#' @param intercept Logical argument indicating whether an intercept should be included. The Default is False.
#' @param contr See model.matrix for more information. The default is effect coded.
#' @return A matrix containig all possible profiles.
#' @export
profiles<-function (lvls, intercept=FALSE, contr='contr.sum'){


  #Make list with all atribute levels
  lvls <-lapply(X = as.list(lvls), function (x) (1:x))

  #generate every combination
  full_design<-as.data.frame(expand.grid(lvls))

  #create factors:
  colnames <- names(full_design)
  full_design[ ,colnames] <- lapply(full_design[ ,colnames], factor)

  #contrast list
  con <- list()

  for(i in 1:length(lvls)){
   name <- paste('Var',i,sep='')
   con[name] <- contr
  }

  #coding
  full_D <- as.data.frame(model.matrix(~ ., full_design, contrasts = con))

  #delete intercept
  if (intercept==F) {full_D <- full_D[,-1]}

  return(as.matrix(full_D))

}

#' Design generation
#'
#' Function to generate a Design matrix given the number of choice sets,
#' the number of alternatives per choice set and a matrix containting possible profiles.
#' @param n_sets Numeric value indicating the number of choide sets.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param profiles Matrix containing (all) possible profiles.
#' @return A design matrix
#' @export
design.gen<-function (lvls, n_sets, n_alts, intercept=FALSE, contr='contr.sum'){


  profs<-profiles(lvls = lvls, contr = contr, intercept = intercept)

  R<-round(runif((n_alts*n_sets), 1, nrow(profs)))
  design<-profs[R,]
  design<-as.matrix(design)

  #return design
  return(design)
}

#' Response generation
#'
#' Function to generate responses given parameter values and a design matrix. Chooses the alternative with max probability.
#' @param par Vector containing parameter values.
#' @param design Design matrix.
#' @return Vector with binary responses.
#' @export
respond<-function (par, design, n_alts){

  d<-design

  p<-t(t(d) * par)
  p<-.rowSums(p, m= nrow(d), n=length(par))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp, rep(seq(1, nrow(d)/n_alts,1), each = n_alts)), each=n_alts)

  if(nrow(d) > n_alts){

    Y<-numeric()
    set<-seq(1,length(p)+n_alts, n_alts)

      for(i in 1:(nrow(d)/n_alts)){

      temp<-rep(0, n_alts)
      temp[which.max( p[set[i]:(set[i+1]-1)] )]<-1

      Y<-c(Y,temp)
    }
  }else{
    Y<-rep(0, length(p))
    Y[which.max(p)]<-1
  }

  return(Y)
}

#' Lattice standard normal generation.
#'
#' Generates a grid of points (N=b^m), coming from a multivariate standard normal distribution with dim = K.
#' @param K Numeric value indicating the dimensionality of the grid.
#' @param b Numeric value indicating the base.
#' @param m Numeric value. N=b^m.
#' @return Matrix of lattice points drawn from a multivariate standard normal distribution.
lattice <- function (K, b = 2, m=10){

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

#' Lattice t-distribution generation.
#'
#' Generates a grid of points (N=b^m), coming from a multivariate t-distribution.
#' @param K Numeric value indicating the dimensionality of the distribution.
#' @param b Numeric value indicating the base.
#' @param m Numeric value. N=b^m.
#' @return Matrix of lattice points drawn from a multivariate t-distribution.
#' @export
lattice_mvt<- function (mode, cvar, df, ...){

  dim<-length(mode)

  #gen lattice from standard normal
  lattice<- lattice(K = dim, ...)

  mean <- t(mode)
  X <- matrix(NA, nrow(lattice), dim)
  A <- chol(cvar)
  t(chol(cvar))

  for(i in 1:nrow(lattice))
  {

    inv<-rgamma(1, df/2)
    invw<-inv/(df/2)
    W<-1/invw

    Z<-lattice[i,]

    r<-mean + sqrt(W) * (Z %*% t(A))
    X[i, ]<-t(r)
  }

  return (X)
}

#' All choice sets
#'
#' Generates all possible combinations of choice sets, given all candidate profiles
#' and the number of alternatives.
#' @param lvls A vector which contains for each attribute, the number of levels.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param mindiff The minimal ammount of atrribute levels that needs to be different.
#' @return Matrix with all possible combinations of profiles, taking into account the mindiff constraint.
#' @export
full_sets<- function(lvls, n_alts, mindiff, intercept= F, contr='contr.sum'){

  cand<-profiles(lvls = lvls)
  fun<-function(x){return(1:x)}
  N<- nrow(cand)

  s1<-rep(N, n_alts)
  s2<-as.list(s1)
  s3<-lapply(X = s2, fun)
  fcomb<-as.data.frame(expand.grid(s3))

  parplace<-lvls-1
  par2<-cumsum(parplace)
  par1<-c(1,par2+1)
  par1<-par1[-length(par1)]
  diff<-numeric()
  newfcomb<-numeric()


  for (s in 1:nrow(fcomb)){

    set <- fp[as.numeric(fcomb[s,]), ]

    for (pp in 1:length(par1)){
    ifelse((sum(abs(diff(set[ , par1[pp] :par2[pp]],1))) !=0), dif<-1, dif<-0)
    diff<-c(diff,dif)
    }

    if (sum(diff) > mindiff){newfcomb<-rbind(newfcomb, fcomb[s, ])}
    diff<-numeric(0)
    }

    return(newfcomb)
}


#roxygen2::roxygenise()
