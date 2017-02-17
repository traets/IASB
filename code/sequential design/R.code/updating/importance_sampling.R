#create design using KL-crit. 
rm(list=ls())

#PRIOR
p_mod<- c( -1, -1, -1, -1, -1, 1)
var<-1
p_cov<- diag(rep(var,6))

P_obs<-1
#samples
N<-1024
sam<-MASS::mvrnorm(n=N, mu=p_mod, Sigma= p_cov)
s<-sam[,P_obs]
w<-rep(1,N)/N


plot(density(x=s, weights = w), col="blue")

#------------
choice_set<-KL_info(lvls=c(3,3,3), par_samples = sam, weights = w, n_alts=2)
set<-choice_set

r<-respond(par= c(-2,-5,-2,-5,-8,2), design = set, n_alts=2)
y<-r

#posterior
impie<-imp_sampling(prior_mode = p_mod, prior_covar = p_cov, design = set, n_alts = 2, Y= y )


###
post_samples<-impie[[1]]
s_p<-post_samples[ ,P_obs]

post_weights<-impie[[2]]

lines (density(x=s_p, weights = post_weights), col="red")

for (i in 1:100)
{

##

choice_set<-KL_info(lvls=c(3,3,3), par_samples = post_samples, weights = post_weights, n_alts=2)
set<-rbind(set, choice_set)

##Respond
r<-respond(par= c(-20,-5,-2,-5,-8,2), design = choice_set, n_alts=2)
y<-cbind(y,r)

#posterior
impie<-imp_sampling(prior_mode = p_mod, prior_covar = p_cov, design = set, n_alts = 2, Y= y )

post_samples<-impie[[1]]
s_p<-post_samples[ ,1]

post_weights<-impie[[2]]

if((i/10)%%1==0) lines(density(x=s_p, weights = post_weights), col="red")

}

imp_sampling


imp_mode <- maxLik::maxNR(logPost, start = p_mod, prior_mode = p_mod, 
                          prior_covar = p_cov, design = set, Y = y, n_alts = 2)$estimate
imp_mode
respond

imp_sampling
likeli