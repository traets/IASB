KL_1 <-
function (lvls, par_samples, weights, n_alts ) #5.7 sec
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
