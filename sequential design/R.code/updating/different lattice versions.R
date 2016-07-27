
rm(list=ls())

# Function to create lattice points
lattice <- function(K, b = 2, m = 9){
  
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
      latt[j, i] <- 1 - abs(2 * e[j, i] - 1)
    }
  }
  latt <- qnorm(latt)
  return(latt)
}

lattice2 <- function(K, b = 2, m = 9){
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

lattice3 <- function(K, m = 6){
  b <- 2 # base, we only work with base = 2 here
  # K is the dimension of the sample
  # R, the size of the sample, is determined by b^m, default = 2^6 = 64
  
  # This routine is based on Hickernell, Hong, L'Ecuyer & Lemieux, "Extensible lattice
  # sequences for quasi-Monte Carlo quadrature, 2000, Society for Industrial and Applied Mathematics,
  # vol 22, no 3, pp 1117-1138
  # This routine constructs an infinite rank-1 lattice sequence in base b with a specific
  # generating vector h and a shift delta, the shift delta randomizes the deterministic sample
  
  # Obtain the shift vector of dimension K
  delta <- runif(K)
  
  # Obtain the generating vector
  eta <- 1571
  # this eta is based on the spectral test (see 4.2 on page 1127) and works well for
  # b = 2, m in 6:12 and dimension K up to 33 with the following specific generating vector h
  # see table 4.1, on page 1127
  h <- eta^(0:(K - 1))
  
  # Obtain phi_b(i) (see page 1121)
  phi <- function(i, m){
    binar <- as.integer(intToBits(i)) # transform i into a bit representation
    binar <- binar[1:m] # number of samples is 2^m, so we require no larger representation
    sum(binar * 2^((-1):(-m)))
  }
  
  # Construct sample
  P <- matrix(nrow = 2^m, ncol = K)
  for(i in 0:(2^m - 1)) P[i + 1, ] <- phi(i, m) * h + delta
  P <- P - floor(P)
  # P <- P %% 1
  
  # Periodizing transformation
  P <- abs(2 * P - 1)
  
  # Transform to a standard normal sample
  P <- qnorm(P)
  
  # Return
  P
}

K <- 5 # dimension
m <- 15 # sample size = 2^m
set.seed(666)
test1 <- lattice(K, m = m)
set.seed(666)
test2 <- lattice2(K, m = m)
set.seed(666)
test3 <- lattice3(K, m = m)

windows()
plot(test1, pch = 19, col = "red", cex = 1, xlab = "", ylab = "")
points(test2, pch = 19, col = "green", cex = 0.75)
points(test3, pch = 19, col = "blue", cex = 0.5)

identical(test1[order(test1[, 1]), ], test2[order(test2[, 1]), ]); identical(test1[order(test1[, 1]), ], test3[order(test3[, 1]), ]); identical(test2, test3)
all.equal(test2, test3)

# The last version is clear in what it does but is not numerically stable.
# The first version is Yies implementation, which uses a weird 'baker' transformation, although
# it doesn't seem to matter with respect to the outcome. The order of the samples is different
# but in the end you average over the samples so the order is not important
# The second version is Yies implementation but with the 'baker' transformation according to the paper

