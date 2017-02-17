# R-Code to estimate the MDFT-version used in Berkowitsch, Scheibehenne, & Rieskamp (2013), JEP: G
# Author: Nicolas Berkowitsch: nicolas.berkowitsch@unibas.ch
# If you are using this code please cite it:

#Berkowitsch, N. A. J., Scheibehenne, B., & Rieskamp, R. (2013). 
#Rigorously Testing Multialternative Decision Field Theory Against Random Utility Models. 
#Journal of Experimental Psychology: General. Advance online publication. doi: 10.1037/a0035159

# Suggestions on how to speed-up the code (e.g., running time) are welcome: nicolas.berkowitsch@unibas.ch

###########################
### Clear Working Space ###
###########################

rm(list=ls(all=T))
graphics.off()

############################
## set working directory ###
############################

workingDir="C:/Users/u0105757/Desktop/DFT/Berkowitsch code/MDFT_Code" # change!
setwd(workingDir)

######################
### set parameters ###
######################

# the next three lines indicate how many attraction, compromisen and similarity choice triplets 
# should be generated
# must be EVEN NUMBER to balance the added optionS
nAttraction  =  6 
nCompromise  =  6  
nSimilarity  =  6 
nTripletsEff =  c(nAttraction,nCompromise,nSimilarity) 
nTrials      =  sum(nTripletsEff)  # indicateS the size of the choice sets (e.g., number of trials) 
nAlt         =  3 # indicate the number of alternatives
nAttr        =  2 # indicate the number of attributes describing the options
sig1         =  1 # set to 0 obtain independent Probit model; 
                  # additionally set phi2 to 1 and phi1 to >1000 (see below)
boundVal     =  5 # indicate the upper (1Eboundval) and lower bound (1-1EboundVal) 
                  # of the estimated parameter (see below)

###########################
### load (example) data ###
###########################

source("nTrialFun.r")  # function to create random attraction, compromise, and similarity choice sets

output = nTrialFun(nTripletsEff = nTripletsEff,  # nTriplets per effect
                           distAttraction = 0.1, # distance of the dominated option to the core options 
                           distComp = 1,         # distance of the extreme option to the core options
                           distSimilarity = 0.1, # distance of the similar option to the core options
                           w1 = 0.5)             # assume some weight for attribute 1

yourDesignData = output[[1]] # all attibute values in a vector
yourObsChoiceData = output[[2]] # assuming "Target" is always chosen
                                # 1 = A, 2 = B, 3 = C)
tripletEffect = output[[3]] # indicates the order of the context effects
                            # (1 = Attraction, 2 = Compromise, 3 = Similarity)


# Or directly replace your design data here:
# Note: to adequately fit the model to the data, more than this six decisions are required
# yourDesignData = c(
#              # Attraction (A = Target)
#              0.5, 0.7, 0.45, 
#              0.5, 0.3, 0.45,
#              # Attraction (B = Target)
#              0.5, 0.7, 0.65,
#              0.5, 0.3, 0.25,
#              # compromise (A = Target)
#              0.5, 0.7, 0.3,
#              0.5, 0.3, 0.7,
#              # Compromise (B = Target)
#              0.5, 0.7, 0.9,
#              0.5, 0.3, 0.1,
#              # Similarity (B = Target)
#              0.5, 0.7, 0.45,
#              0.5, 0.3, 0.55,
#              # Similarity (A = Target)
#              0.5, 0.7, 0.75,
#              0.5, 0.3, 0.25
#              )

# yourObsChoiceData = c(1,2,1,2,2,1) # assuming "Target" is always chosen
#                                    # 1 = A; 2 = B; 3 = C
             
# transforms your data-vector into an array
# M contains the attribute values of the options
M = array(data = yourDesignData, dim = c(nAlt,nAttr,nTrials),
               dimnames = list(
                 paste("Alt_",1:nAlt, sep = ""),
                 paste("Attr_", 1:nAttr, sep = ""),
                 paste("Trial_", 1:nTrials, sep = "")))

# visualization of the (generated) choice triplets
plotM_fun = function(){
  labels = LETTERS[1:nAlt]
  par(mfrow = ceiling(c(sqrt(nTrials),sqrt(nTrials))))
  for (i in 1:dim(M)[3]) {
    plot(M[,,i], type = "n", xlim = c(0,1), ylim = c(0,1))
    text(M[,,i], labels, col = c("black","red","dark green"), cex = 1)
    title(paste("Trial_", i,sep = ""))
  }
}
windows()
plotM_fun()

#####################
## load libraries ###
#####################

library(mnormt) # for calculation multinomial density functions
library(pso) # optimization algorythm

################################
### load necessary functions ###
################################

# function for t-->Inf; calculates Xi(Inf) and Omega(Inf)
source("omegaNoLoop.r") 

#load the generalized distance function from Berkowitsch, Scheibehenne, Rieskamp, & Matthäus (2013)
source("GenDistFun.r")

#create the L-array and the C-Matrix (see Roe, Busemeyer, & Townsend, 2001)
source("createCandL.r") 

######################
### MDFT-functions ###
######################

# choice function
choice_fun <- function (L, eta, omega) { # start function: choice
  
  Leta = L%*%eta # mean difference in preference
  Lomega = L%*%omega%*%t(L) # varance-covariance matrix of the mean difference in preference
  p = pmnorm(x=Leta,varcov=Lomega)
  return(p)
} # end function: choice

pMDFT_fun <- function (M, S, W, sig2) {
  
  # checks whether NAs or negative values were assigned to the parameter. If so, assign a bad fit
  if(NA %in% c(W,sig2)==TRUE || TRUE %in% cbind(W<0,sig2<0)==TRUE) {
    p = rep(1/1000,nAlt) # assign a bad fit
  }
  
  else { # the core of MDFT
    EV=C%*%M%*%W   # expected value for each option
    eta = EV
    psi = diag(W, nrow=length(W)) - W%*%t(W) # variance-covariance matrix for the primary weights
    nM = M / norm(M); # prevents that R assign in phi very small values 0
    nC = C / norm(C); # prevents that R assign in phi very small values 0
    # variance-covariance matrix of eta
    # multiply with norm to obtain same values
    phi= (sig1*nC%*%nM%*%psi%*%t(nM)%*%t(nC))*norm(M)^2*norm(C)^2 + sig2*diag(nAlt) 
    
    omega = phi + 0*diag(nAlt)
    
    # t goes to infinity
    omega = omegaNoLoop(S,phi,nAlt)
    if(all(is.na(omega)) == T) {
      p = rep(1/1E+10,nAlt) #assign a bad fit
    } else {
      eta = try(solve((diag(nAlt)-S),EV),silent = T)        
      p = apply(L,c(nAlt),choice_fun, eta=eta,omega=omega)          
    }
  }
  return(p)
} #end function: pMDFT_fun
?apply
pMDFT_choiceSet_fun = function(W,sig2,phi2,phi1,nTrials,M,wgt) {
  
  pMDFT = matrix(NA,nTrials,nAlt) # empty matrix to speed up
  dimnames(pMDFT) = list(paste("trial_",1:nTrials,sep=""),paste("alt_",1:nAlt,sep = ""))
  
  # Calculate the feedback array S and its eigenvalues 
  S = array(apply(M,3,distfunct, phi1 = phi1, phi2 = phi2, W = W, wgt = wgt), dim = c(nAlt,nAlt,nTrials))
  eigenValS = matrix(unlist(apply(S,3,eigen, only.values = T)),nTrials,nAlt, byrow = T)

  # function assigns a bad fit if eigenvalues of S are not between 0 and 1 or if there are any NAs in S
  eigenValCheck_fun = function(S,eigenValS,nTrials,nAlt){

    if(TRUE %in% cbind(eigenValS<=0,eigenValS>=1)==TRUE || TRUE %in% is.na(S)==TRUE) {
      eigenValCheck = FALSE
    } else {
        eigenValCheck = TRUE
      }
    return(eigenValCheck)
  }
  
  if(eigenValCheck_fun(S,eigenValS,nTrials,nAlt) == TRUE){
    pMDFT = t(sapply(1:nTrials,function(x)pMDFT_fun(M[,,x],S[,,x],W,sig2))) # calculated the choice probabilities
  } else {
    pMDFT = matrix(rep(1/1E+10,nTrials*nAlt),nTrials,nAlt) # assign a bad fit
  }
  
  return(pMDFT)
} # end function: pMDFT_choiceSet_fun

optMDFT_fun = function(par){ # start optimization function
  
  # assure the attribute weights sum to 1
  Wtemp = sort(par[1:(nAttr-1)])
  Wtemp2 = c(0,Wtemp,1)
  W = diff(Wtemp2)  
  sig2 = par[nAttr]
  phi2 = par[nAttr+1]
  phi1 = par[nAttr+2]
  wgt = par[nAttr+3]
  
  # check for negative or missing values
  if(NA %in% c(W,sig2,phi2,phi1)==TRUE || TRUE %in% cbind(W<0,sig2<0,phi2<0,phi1<0)==TRUE) {
    LL_MDFT_sum = 1000 # assign a bad fit
  }
  
  pMDFT_allAlt = pMDFT_choiceSet_fun(W,sig2,phi2,phi1,nTrials,M,wgt)  
  pMDFT_chosenAlt = pMDFT_allAlt[cbind(1:nTrials,yourObsChoiceData)]  
  # prevent Inf values when taking the log of pMDFT_chosenAlt 
  pMDFT_chosenAlt[pMDFT_chosenAlt == 0]   =.0000000001
  pMDFT_chosenAlt[pMDFT_chosenAlt == +Inf]=.999
  # calculate the summed Log-Likelihood (LL)
  # multiplicate with -1 to minimize positive values
  LL_MDFT_sum = try(-1*sum(log(pMDFT_chosenAlt)), silent=TRUE) 
  
  return(LL_MDFT_sum)  
} # end function: optMDFT_fun 

#######################################
### call the optimization procedure ###
#######################################

# generate random starting values
initParW = runif(nAttr-1,.0001,.9999)
initParSig2 = runif(1,.0001,999)  
initParPhi2Phi1 = c(.05,runif(1,10,100)) # seems to be good starting values
initParWgt = runif(1,1,12) # seems to be good starting values
initPar = c(initParW,initParSig2,initParPhi2Phi1,initParWgt)

# use different optimization algorhythms:
# 1. nlminb
valuesMDFT = nlminb(start=initPar, objective = optMDFT_fun, 
                       lower = c(rep(0.1^boundVal,length(initPar)-1),1), 
                       upper = c(rep(1-0.1^boundVal,nAttr-1),0.1^-boundVal,1-0.1^boundVal,rep(0.1^-boundVal,2)),
                       control = list(trace = T, rel.tol = 1e-15))
# 2. psoptim (use best values from nlminb as starting parameters)
valuesMDFT = psoptim(par=valuesMDFT$par, fn = optMDFT_fun, 
                         lower = c(rep(0.1^boundVal,length(initPar)-1),1), 
                         upper = c(rep(1-0.1^boundVal,nAttr-1),0.1^-boundVal,1-0.1^boundVal,rep(0.1^-boundVal,2)),
                        control = list(trace = 6, reltol = 1e-15, maxit=400, REPORT = 20))
# 3. nlminb (use best values from psoptim as starting parameters)
valuesMDFT = nlminb(start=valuesMDFT$par, objective = optMDFT_fun, 
                        lower = c(rep(0.1^boundVal,length(initPar)-1),1), 
                        upper = c(rep(1-0.1^boundVal,nAttr-1),0.1^-boundVal,1-0.1^boundVal,rep(0.1^-boundVal,2)),
                        control = list(trace = T, rel.tol = 1e-15))

# valuesMDFT = nlminb(start= initPar, objective = optMDFT_fun, control = list(trace = T))
# valuesMDFT = psoptim(par=valuesMDFT$par, fn = optMDFT_fun, control = list(trace = T, maxit = 600, REPORT = 50))

bestPar = valuesMDFT$par

Wtemp = sort(bestPar[1:(nAttr-1)])
Wtemp2 = c(0,Wtemp,1)
W = diff(Wtemp2)  
sig2 = bestPar[nAttr]
phi2 = bestPar[nAttr+1]
phi1 = bestPar[nAttr+2]
wgt = bestPar[nAttr+3]

# fit of MDFT 
LL_MDFT = -1*valuesMDFT$objective
# fit of random model
LL_Random = log(1/nAlt)*nTrials
LL_Models = c(LL_MDFT,LL_Random)
names(LL_Models) = c("LL_MDFT", "LL_Random")


# calculate probabilities based on estimated parameters
pPred = pMDFT_choiceSet_fun(W,sig2,phi2,phi1,nTrials,M,wgt)  
pPred
# probabilities of chosen options:
pObsPred = pPred[cbind(1:nTrials,yourObsChoiceData)]
cbind(pObsPred,tripletEffect)

# show parameter values
bestPar
# show model fit
LL_Models


