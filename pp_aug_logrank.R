#
# pp_aug_logrank.R
#---------------------------------------------------------------------------
# Calculate the predicted power of randomized two groups 
#   with the method by Hattori, Komukai and Friede
#
#---------------------------------------------------------------------------
# hr   : hazard ratio to detect
# p    : allocation ratio (=0.5 if 1:1 randomization)
# L    : the number of events
# n    : the number of subjects
# dat  : dataset
# time : survival time
# delta: survival indicatior (1: event, 0: censored)
# cov  : covariates
#
# [Remark] 
#  n should be determined accounting for accrual and follow-up durations.
#  See Section 3.1 in Hattori, Komukai and Friede.
#--------------------------------------------------------------------------
#  27Dec2021 Satoshi Hattori
#---------------------------------------------------------------------------

library( survival )

pp_aug_logrank       <- function( hr, p, L, n, time, delta, cov ){

	beta                 <- log( hr )

	cox_null             <- coxph( formula = Surv( time, delta ) ~ 1 )

	mcov                 <- cox_null$residuals * cov
	mmcov                <- apply( mcov, 2, mean )

	op                   <- function( x ){ outer( x, x ) }
	qq                   <- matrix( apply( apply( cov, 1, op ), 1, mean ), nrow = ncol( cov ), ncol = ncol( cov ) )
	e2                   <- t( as.matrix( mmcov ) ) %*% solve( qq ) %*% as.matrix( mmcov ) * p * ( 1 - p )

	LP                   <- L * p * ( 1 - p )
	v2                   <- 1 / LP * ( LP - n * e2 ) * 1 / LP
	v2_logrank           <- 1 / LP

	power_LT             <- pnorm( - 1.96 - beta / sqrt( v2 ) ) + 1 - pnorm( 1.96 - beta / sqrt( v2 ) )
	power_logrank        <- pnorm( - 1.96 - beta / sqrt( v2_logrank ) ) + 1 - pnorm( 1.96 - beta / sqrt( v2_logrank ) )

	output               <- list( hr, L, n, e2, power_LT, power_logrank )
	names( output )      <- c( "HR", "L", "n", "e2", "Power(augmented logrank)", "Power(logrank)" )

	return( output )
}




