#####################################################################
# MODEL APPLICATION (FORECASTING)                                   #
#####################################################################

number_of_alts=2

beta_for_forecast=est # if using the final estimates from the previous model run

# DEFINE MODEL STRUCTURE (NEEDS TO BE THE SAME AS IN ESTIMATION - BUT SAVES EACH CHOICE TASK SEPARATELY)

forecast=function(beta)
{
  # vector to save probabilities in

  probs_out<-data.table(matrix(0,nrow(data),number_of_alts+3))
  
  setnames(probs_out,c("ID","choice","P1","P2","P_choice"))
  
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
  
  # calculate probability for all alternatives
    
  probs_out[,ID:=data[, ID]]
  probs_out[,choice:=data[, choice]]
  probs_out[,P1:=data[, .(U1/(U1+U2))]]
  probs_out[,P2:=data[, .(U2/(U1+U2))]]
  probs_out[,P_choice:=(choice==1)*P1+(choice==2)*P2]
  
  
  # required to avoid double attachment of names
  
  detach(beta1)
 
  # return the calculated probabilities
  
  return(probs_out)
}

#####################################################################
# RUN THE PREDICTION                                                #
#####################################################################

predictions=forecast(beta_for_forecast)

# calculate lowest half percentile for probabilities for actual choice

limit=quantile(predictions[,P_choice],0.005)

# print out any falling below that

outliers=predictions[P_choice<limit]

setorder(outliers,P_choice)

outliers

