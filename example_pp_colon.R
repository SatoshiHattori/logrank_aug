#
# example_pp_colon.R
#---------------------------------------------------------------------------
# Illustration of the predicted power calculation
# 
# The obserational group of COLON dataset in survival package of R is used
#  as the reference data to calculate the predicted power
#
#---------------------------------------------------------------------------
#  28Dec2021 Satoshi Hattori
#---------------------------------------------------------------------------


rm( list = ls() )
setwd( "F:\\logrank_adj_R1\\web_supplementary" )
getwd()

# Calling the R function for calculating the predicted power
source( "pp_aug_logrank.R" )


# Inputting a reference data
library( survival )
colon2         <- colon[ colon$rx == "Obs", ]

# Defining covariates in the augmentation term
attach( colon2 )
covvec         <- cbind( sex, node4, age )
detach( colon2 )

time           <- colon2$time
status         <- colon2$status

result         <- pp_aug_logrank( 0.8, 0.5, 640, 1490, time, status, covvec )

result



