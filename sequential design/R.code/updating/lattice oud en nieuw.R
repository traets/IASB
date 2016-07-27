lattice_old <- function(K, b = 2, m = 9){
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
    ei <- tcrossprod(crossprod(seq, base(i - 1)), av) + u
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

lattice_new <- function(K, b = 2, m = 9){
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
    ei <- crossprod(seq, base(i - 1)) * av + u
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