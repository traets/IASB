KL_original <-
function (lvls, par_samples, weights, n_alts) 
{
  
  fp <- profiles(lvls)
  fcomb <- full_sets(cand = fp, n_alts = n_alts)
  rows <- apply(fcomb[, 1:ncol(fcomb)], 1, function(i) length(unique(i)) > 
                  ncol(fcomb) - 1)
  fcomb <- fcomb[rows, ]
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
