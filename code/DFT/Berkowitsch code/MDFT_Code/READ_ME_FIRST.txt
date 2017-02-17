Berkowitsch, N. A. J., Scheibehenne, B., & Rieskamp, R. (2013). 
Rigorously Testing Multialternative Decision Field Theory Against Random Utility Models. 
Journal of Experimental Psychology: General. Advance online publication. doi: 10.1037/a0035159


###############################
# Read me file for the R-Code #
###############################

Main File: "R-Code_to_Berkowitsch_et_al._(2013).r"

Put the following files in the same folder as the main file:
- "createCandL.r"
- "GenDistFun.r"
- "nTrialFun.r"
- "omegaNoLoop.r"

#########################################
# R-Code_to_Berkowitsch_et_al._(2013).r #
#########################################

- Contains the main code of MDFT to calculate
  - the choice probabilities
  - the model fit of MDFT and a random model
  - the model parameters
- Set here the number of 
  - alternatives
  - attributes
  - ...
- Calls the optimization procedure

#################
# createCandL.r #
#################

- Function which creates the contrast matrix C and the L-matrix (for details see Roe et al. 2001)

################
# GenDistFun.r #
################

- Generalized Distance Function (for details see Appendix B in Berkowitsch, Scheibehenne, & Rieskamp, 2013 and Berkowitsch, Scheibehenne, Rieskamp, Matthäus, 2014)

###############
# nTrialFun.r #
###############

- Function which creates random attraction, compromise, and similarity choice triplets
- for the paper the choice sets were generated using first a matching task (for details see method section of Study 2 in Berkowitsch, Scheibehenne, & Rieskamp, 2013)

#################
# omegaNoLoop.r #
#################

- Function for t-->Inf; calculates Xi(Inf) and Omega(Inf)
- Avoids loopings up to t