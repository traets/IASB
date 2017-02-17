##################################
# updating the prior             #
# based on initial choice sets   #
##################################


#generate lattice points from standard normal 
lattice2 <- function(K, b, m){
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
  for(i in 2:K) av[i] <- (a * av[i - 1]) %% N
  e <- matrix(0, N, K)
  seq <- rep(0, m)
  kk <- -m
  for(k in 1:m){ seq[k] <- b^kk; kk <- kk + 1;}
  for(i in 1:N){
    ei <- tcrossprod(as.vector(crossprod(seq, base(i - 1))), av) + u
    e[i, ] <- ei - floor(ei)
  }
  latt <- matrix(0, N, K)
  for(i in 1:K){
    for(j in 1:N){
      latt[j, i] <- abs(2 * e[j, i] - 1) # This is a different "bakers" transformation, differs
      # between Hickernell et al. and Yu her thesis
    }
  }
  latt <- qnorm(latt)
  return(latt)
}

#compute the mode and the hessian matrix H of the posterior distribution (newton raphson) 


#importance distribution = multivariate t distribution with mean = mode post


#use lattice points as draws from importance distribution


#compute the weight for each draw 







