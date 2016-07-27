rm(list = ls())


##### specify #######################################################
prior_mode<- c(1,-2,3,0.4)
prior_cov<- diag(c(3,3,3,3))
true_par<-c(5,-5, 2, -2)

lvs<-c(2,3,2)
n_alts<-2

initial_sets<-10
new_sets<-2
#####################################################################


#prior samples
p_samples<-MASS::mvrnorm(n = 500, mu= prior_mode, Sigma=prior_cov)

#generate start design
design<-DB_mod_fed(lvls = lvs , n_sets= initial_sets, n_alts = n_alts, par_samples = p_samples, max_iter = 3 )
design
des<-design[[1]]

#generate responses
Y<-respond(par= true_par, design = des, n_alts = 2)


for (i in 1:new_sets){

  #update parameters
  draws<-imp_sampling(prior_mode = prior_mode, prior_covar = prior_cov, design = des, n_alts = n_alts, Y=Y )
  w<-draws[[2]]
  sam<-draws[[1]]

  #new set based on updated parameters
  new_set<-choice::DB_seq_fed(design = des, lvls = lvs, n_alts = n_alts, par_samples = sam,
                              weights = w, prior_covar = prior_cov)

  #respond
  Y_set<-respond(par= true_par, design = new_set, n_alts = n_alts)

  #update
  des <- rbind(des,new_set)
  Y<- c(Y, Y_set)
  prior_mode<-modes(samples = sam, weights = w, s=10)

 }









