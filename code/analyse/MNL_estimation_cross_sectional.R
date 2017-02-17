source("function_input.R")

#####################################################################
# DATA INPUT                                                        #
#####################################################################

# define name of data file - note this needs to include a variable  
# called ID, with respondent IDs (with all observations for the same respondent grouped together)                          

datafilename="danish_data.dat"

# read data file

data=read.table(datafilename,header=T)

# OPTIONAL: apply exclusions (use subset of data)
# example: data=subset(data,data$sex==1)

# turn into data table (for efficient estimation)

data=data.table(data)

# determine number of observations in the data, treating it as cross-sectional

N=nrow(data)

#####################################################################
# DEFINE MODEL PARAMETERS                                           #
#####################################################################

# set starting values

beta=c(0,0,0)

# set names for parameters

names(beta)=c("asc1","tc","tt")

# OPTIONAL: set fixed parameters (leave empty if none, i.e. fixedparams=c())

fixedparams=c("asc1")

#####################################################################
# DEFINE MODEL STRUCTURE AND LIKELIHOOD FUNCTION                    #
#####################################################################

loglike=function(beta)
{
   # needed to be able to refer to parameters by name

   beta1=as.list(beta)
   attach(beta1)

   # define utility functions
   
   data[,U1:=asc1+tc*tc1+tt*tt1]
   data[,U2:=tc*tc2+tt*tt2]
   
   # avoid numerical issues
   
   data[U1>700,U1:=700]
   data[U2>700,U2:=700]
   data[U1< -700,U1:=-700]
   data[U2< -700,U2:=-700]
   
   # exponentiate utilities
   
   data[,U1:=exp(U1)]
   data[,U2:=exp(U2)]
   
   # calculate probability for chosen alternative for each observation
   
   data[,P:=((choice==1)*U1+(choice==2)*U2)/(U1+U2)]
   
   # take product across choices for the same person (likelihood)
   # then take the log for log-likelihood 
   
   LL <- data[, .(out_LL = log(P))][["out_LL"]]

   # remove beta names from memory so as to avoid double attachment of names in next iteration
  
   detach(beta1)
  
   # return the log-likelihood with current estimates as output
  
   return(LL)
}

#####################################################################
# RUN THE MODEL                                                     #
#####################################################################

model=maxLik(loglike,start=beta,fixed=fixedparams,method="BFGS",print.level=3)

# formatted output

modeloutput(model)

# example calculation of VTT

est["tt"]/est["tc"]*60
