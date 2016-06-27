

#' DB modified federov algorithm
#'
#' Modified federov algorithm. Swipes every profile of an initial design with candidate profiles.
#' By doing this it tries to minimize the DB error. The DB error is te mean D-error over samples
#' from the prior parameter distribution.
#' @param lvls A vector which contains for each attribute, the number of levels.
#' @param n_sets Numeric value indicating the number of choice sets.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param par_samples A matrix in which each row is a sample from the multivariate prior parameter distribution.
#' @param max_iter A numeric value indicating the maximum number allowed iterations.
#' @return An efficient design, and the associated DB error.
#' @export
DB_mod_fed <- function(lvls, n_sets, n_alts, par_samples, max_iter=Inf, ...){

  #Generation
  profiles<- profiles(lvls, ...)
  des<- design(profiles= profiles, n_sets, n_alts)

  #check
  if(ncol(des)!= ncol(par_samples)){
    stop('Number of paramters does not match the design matrix')
  }

  if((nrow(des)/n_alts)%%1!=0){
    stop('number of alternatives per choice set does not match the design')
  }


  #start value
  DB_start <- 99999
  converge <- FALSE
  it <- 1

  while (!converge & it<=max_iter){

    it <- it+1
    #Progress bar
    pb <- txtProgressBar(min = 0, max = nrow(des), style = 3)

    #set round design
    round_des<-des

    for (trow in 1:nrow(des)){

      for (cand in 1:nrow(profiles)){

        #Change alternative
        des[trow, ]<- profiles[cand, ]

        #calculate for all samples the D error.
        D_errors<-apply(par_samples, 1, d_err, design = des, n_alts = n_alts)

        #DB error
        DB<-mean(D_errors)
        if (is.nan(DB)){stop('NaN produced in DB error')}

        #if better --> keep
        if (DB < DB_start){
          best_row<- as.numeric(des[trow, ])
          DB_start<-DB
          change<-T
        }

        Sys.sleep(0.1)
        # update progress bar
        setTxtProgressBar(pb, trow)

      }
      # Changes after loop through one row
      if (change){
        des[trow, ] <- best_row
      }else{
        des[trow, ] <- round_des[trow, ]
      }

      change<-F
    }

    #Changes after loop through whole design
    close(pb)
    converge<-isTRUE(all.equal(des, round_des))
  }

  return(list(des, DB))
}


#' DB sequential set
#'
#' Adds the most DB efficient choice set to a design, given (updated) parameter values.
#' @param design A design matrix.
#' @param lvls A vector which contains for each attribute, the number of levels.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param par_samples A matrix in which each row is a sample from the multivariate prior parameter distribution.
#' @param weights A vector containing the weights of all the samples.
#' @param prior_covar covariance matrix of the prior distribution.
#' @return The most DB efficient choice set
#' @export
DB_seq_fed <-function(design, lvls, n_alts, par_samples, weights, prior_covar, ... ){

  #start values
  DB_best<-9999
  candis<-profiles(lvls, ...)
  full_comb<- full_sets(cand = candis, n_alts = n_alts)


  #for each potential set:
  for (p in 1:nrow(candis)){

    set<-as.matrix(candis[as.numeric(full_comb[p, ]), ])

    #for each par draw calculate d-error:
    DB_error<-0

    for (s in 1: nrow(par_samples)){

      info_s<-info_set(par = par_samples[s, ], c_set= set)
      info_d<-info_design(par= par_samples[s, ], design = design, n_alts = n_alts)
      inv_cov_prior<-solve(prior_covar)

      DB_error<- DB_error+(det(info_d + info_s + inv_cov_prior)*(-1/ncol(par_samples))*weights[s])

    }

    #select new set if smaller DB error
    if (DB_error< DB_best){final_set<-set; DB_best<- DB_error}
  }

  return(final_set)

}
