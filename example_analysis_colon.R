#
# example_analysis_colon.R
#---------------------------------------------------------------------------
# Illustration of the augmented logrank test
# 
# The adjuvant chemotherapy groups of COLON dataset in survival package of R
# is used as the target study
#
#---------------------------------------------------------------------------
#  28Dec2021 Sho Komukai, Satoshi Hattori
#---------------------------------------------------------------------------


rm( list = ls() )
setwd( "F:\\logrank_adj_R1\\web_supplementary" )
getwd()


source( "aug_logrank.R" )


# Inputting and creating the target study dataset
library( survival )

colon2         <- colon[ colon$rx != "Obs", ]

colon2$group[ colon2$rx == "Lev" ]      <-0
colon2$group[ colon2$rx == "Lev+5FU" ]  <-1

#-----analysis of covariate-adjusted logrank test------#

attach( colon2 )

Covariate      <- data.frame( sex, node4, age )
q.f            <- function( x ){ data.frame( 1, x[ , 1 ], x[ , 2 ],  x[ , 3 ] ) } 

adj_cox1       <- covariate.adjusted2( time, status, group, Covariate, q.f ) 

detach( colon2 )

hr_adj         <- exp( adj_cox1[[ 2 ]][[ 1 ]][ 1 ] )
ci_low_adj     <- exp( adj_cox1[[ 2 ]][[ 1 ]][ 1 ] - 1.96 * adj_cox1[[ 2 ]][[ 1 ]][ 3 ] )
ci_upp_adj     <- exp( adj_cox1[[ 2 ]][[ 1 ]][ 1 ] + 1.96 * adj_cox1[[ 2 ]][[ 1 ]][ 3 ] )
ci_adj         <- cbind( hr_adj, ci_low_adj, ci_upp_adj )
pvalue         <- adj_cox1[[ 2 ]][[ 1 ]][ 13 ]

pvalue
ci_adj





