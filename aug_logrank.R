#
# aug_logrank.R
#---------------------------------------------------------------------------
# Augmented logrank test of randomized two groups 
#   with the method by Hattori, Komukai and Friede
#
#---------------------------------------------------------------------------
# covariate.adjusted2( time, delta, Z, Covariate, q.f): main function
# time     : survival time
# delta    : survival indicatior (1: event, 0: censored)
# Z        : treatment indicator
# Covariate: covariates
# q.f      : function "q" of equation (3) in the main paper
#
#--------------------------------------------------------------------------
#  28Dec2021 Sho Komukai, Satoshi Hattori
#---------------------------------------------------------------------------


############################################################################################
beta_PH.function   <- function( n, beta.size, tau.size, delta.d, mat.Z.d, mat.Z.d.cro, int.ind, count.ind, risk.ind, t.risk.ind ){

	cro          <- function(x)x%x%x

	beta0        <- rep( 0, beta.size )
	bz           <- as.numeric( mat.Z.d %*% beta0 )
	exp.bz       <- exp( bz )

	s0tZ         <- as.data.frame( t.risk.ind %*% as.matrix( exp.bz ) )
	s1tZ         <- as.data.frame( t.risk.ind %*% ( mat.Z.d * exp.bz ) )
	s2tZ         <- as.data.frame( t.risk.ind %*% as.matrix( mat.Z.d.cro * exp.bz ) )
	logs0tZ      <- as.data.frame( log( s0tZ ) )
	s1tZ_s0tZ    <- as.data.frame( s1tZ/as.matrix( s0tZ ) )
	s1tZ_s0tZ2   <- as.data.frame( if( beta.size == 1 ){ s1tZ_s0tZ^2 }else{ t( apply( s1tZ_s0tZ, 1, cro ) ) } )
	s2tZ_s0tZ    <- as.data.frame( s2tZ/as.matrix( s0tZ ) )

	logL.NULL    <- apply( delta.d * ( exp.bz - int.ind %*% as.matrix( logs0tZ ) ), 2, sum )
	Score.NULL   <- apply( delta.d * ( mat.Z.d - int.ind %*% as.matrix( s1tZ_s0tZ ) ), 2, sum )
	Fisher.NULL  <- matrix( apply( delta.d * as.data.frame( int.ind %*% ( as.matrix( s2tZ_s0tZ - s1tZ_s0tZ2 ) ) ), 2, mean ), beta.size )

	wi_integrand <- 	aperm( 
				 array( rep( mat.Z.d, tau.size ), dim = c( n, beta.size, tau.size ) ), perm = c( 1, 3, 2 ) 
				) - array( 
				 rep( as.matrix( s1tZ_s0tZ ), rep( n, tau.size * beta.size ) ), dim = c( n, tau.size, beta.size ) 
				)
	wi_d         <- 	array( 
				 rep( count.ind - 
				  t( 
				   as.data.frame( 
				    t( risk.ind * exp.bz ) 
				   )/as.matrix( s0tZ/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
				  ), beta.size 
				 ), dim = c( n, tau.size, beta.size ) 
				)
	wi           <- apply( wi_integrand * wi_d, c( 1, 3 ), sum )

	Ewi2.NULL    <- if( beta.size == 1 ){ 
				 matrix( mean( wi^2 ) ) 
				}else{ 
				 (
				  matrix( 
				   apply( 
				    t( apply( wi, 1, cro ) ), 2, mean 
				   ), ncol = beta.size 
				  ) 
				 )
				}

	beta.new     <- as.numeric( beta0 + solve( Fisher.NULL * n ) %*% Score.NULL )
	for(j in 1:100){
		beta.old     <- beta.new

		bz           <- as.numeric( mat.Z.d %*% beta.old )
		exp.bz       <- exp( bz )

		s0tZ         <- as.data.frame( t.risk.ind %*% as.matrix( exp.bz ) )
		s1tZ         <- as.data.frame( t.risk.ind %*% ( mat.Z.d * exp.bz ) )
		s2tZ         <- as.data.frame( t.risk.ind %*% as.matrix( mat.Z.d.cro * exp.bz ) )
		logs0tZ      <- as.data.frame( log( s0tZ ) )
		s1tZ_s0tZ    <- as.data.frame( s1tZ/as.matrix( s0tZ ) )
		s1tZ_s0tZ2   <- as.data.frame( if( beta.size == 1 ){ s1tZ_s0tZ^2 }else{ t( apply( s1tZ_s0tZ, 1, cro ) ) } )
		s2tZ_s0tZ    <- as.data.frame( s2tZ/as.matrix( s0tZ ) )

		Score        <- apply( delta.d * ( mat.Z.d - int.ind %*% as.matrix( s1tZ_s0tZ ) ), 2, sum )
		Fisher       <- matrix( apply( delta.d * as.data.frame( int.ind %*% ( as.matrix( s2tZ_s0tZ - s1tZ_s0tZ2 ) ) ), 2, sum ), beta.size )/n

		beta.new     <- as.numeric( beta.old + solve( Fisher*n ) %*% Score )
		if( sum( abs( c( beta.new - beta.old ) ) ) < 0.00001 )break
	}

	beta.Hat     <- beta.new
	bz           <- as.numeric( mat.Z.d %*% beta.Hat )
	exp.bz       <- exp( bz )

	s0tZ         <- as.data.frame( t.risk.ind %*% as.matrix( exp.bz ) )
	s1tZ         <- as.data.frame( t.risk.ind %*% ( mat.Z.d * exp.bz ) )
	s2tZ         <- as.data.frame( t.risk.ind %*% as.matrix( mat.Z.d.cro * exp.bz ) )
	logs0tZ      <- as.data.frame( log( s0tZ ) )
	s1tZ_s0tZ    <- as.data.frame( s1tZ/as.matrix( s0tZ ) )
	s1tZ_s0tZ2   <- as.data.frame( if( beta.size == 1 ){ s1tZ_s0tZ^2 }else{ t( apply( s1tZ_s0tZ, 1, cro ) ) } )
	s2tZ_s0tZ    <- as.data.frame( s2tZ/as.matrix( s0tZ ) )

	logL.Hat     <- apply( delta.d * ( exp.bz - int.ind %*% as.matrix( logs0tZ ) ), 2, sum )
	Score.Hat    <- apply( delta.d * ( mat.Z.d - int.ind %*% as.matrix( s1tZ_s0tZ ) ), 2, sum )
	Fisher.Hat   <- matrix( apply( delta.d * as.data.frame( int.ind %*% ( as.matrix( s2tZ_s0tZ - s1tZ_s0tZ2 ) ) ), 2, mean ), beta.size )

	wi_integrand <- 	aperm( 
				 array( rep( mat.Z.d, tau.size ), dim = c( n, beta.size, tau.size ) ), perm = c( 1, 3, 2 ) 
				) - array( 
				 rep( as.matrix( s1tZ_s0tZ ), rep( n, tau.size * beta.size ) ), dim = c( n, tau.size, beta.size ) 
				)
	wi_d         <- 	array( 
				 rep( count.ind - 
				  t( 
				   as.data.frame( 
				    t( risk.ind * exp.bz ) 
				   )/as.matrix( s0tZ/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
				  ), beta.size 
				 ), dim = c( n, tau.size, beta.size ) 
				)
	wi           <- apply( wi_integrand * wi_d, c( 1, 3 ), sum )

	Ewi2         <- if( beta.size == 1 ){ 
				 matrix( mean( wi^2 ) ) 
				}else{ 
				 (
				  matrix( 
				   apply( 
				    t( apply( wi, 1, cro ) ), 2, mean 
				   ), ncol = beta.size 
				  ) 
				 )
				}

	v.cov.beta   <- solve( Fisher.Hat*n )                                   ## estimator
	v.cov2.beta  <- solve( Fisher.Hat ) %*% Ewi2 %*% solve( Fisher.Hat )/n  ## Robust estimator

	output_NULL  <- list( beta0, logL.NULL, Score.NULL, Fisher.NULL, Ewi2.NULL ) ; names( output_NULL ) <- c( "beta", "logL", "Score", "Fisher", "Ewi2" )
	output_Hat   <- list( 
				 beta.Hat, bz, exp.bz, 
				 s0tZ, s1tZ, s2tZ, logs0tZ, s1tZ_s0tZ, s1tZ_s0tZ2, s2tZ_s0tZ, 
				 logL.Hat, Score.Hat, Fisher.Hat, 
				 wi, Ewi2, v.cov.beta, v.cov2.beta
				)
	names( output_Hat ) <- c( 
					 "beta", "bz", "exp.bz", 
					 "s0tZ", "s1tZ", "s2tZ", "logs0tZ", "s1tZ_s0tZ", "s1tZ_s0tZ2", "s2tZ_s0tZ", 
					 "logL", "Score", "Fisher", 
					 "wi", "Ewi2", "Var", "RobustVar"
					)
	output       <- list( output_NULL, output_Hat ) ; names( output ) <- c( "outputNULL", "outputHat" )
	return( output )
}
############################################################################################
covariate.adjusted2 <- function( time, delta, Z, Covariate, q.f){

	#-------------------#
	#--- Preliminary ---#
	#-------------------#
	cro          <- function(x)x%x%x

	time.ind     <- order( time )
	time.d       <- time[ time.ind ]
	delta.d      <- delta[ time.ind ]
	Z.d          <- as.matrix( Z )[ time.ind, ]
	Covariate.d  <- Covariate[ time.ind, ]
	n            <- length( time )

	mat.Z.d      <- as.matrix( Z.d )
	beta.size    <- length( mat.Z.d[ 1, ] )
	mat.Z.d.cro  <- if( beta.size == 1 ){ mat.Z.d^2 }else{ t( sapply( 1:n, function( x )mat.Z.d[ x, ] %x% mat.Z.d[ x, ] ) ) }

	df.Z.d       <- as.data.frame( Z.d )
	df.Z.d.cro   <- as.data.frame( mat.Z.d.cro )

	tau.d        <- unique( time.d[ which( delta.d == 1 ) ] )
	tau          <- tau.d[ order( tau.d ) ]

	tau.size     <- length( tau )

	tau.dd       <- c( tau, max( time.d )+10 )
	int.ind      <- sapply( 1:tau.size, function( x ) ifelse( time.d >= tau.dd[ x ], ifelse( time.d < tau.dd[ x+1 ], 1, 0 ), 0 ) )
	risk.ind     <- sapply( 1:tau.size, function( x ) ifelse( time.d >= tau[ x ], 1, 0 ) )
	count.ind    <- sapply( 1:tau.size, function( x ) ifelse( time.d == tau[ x ], 1, 0 ) )
	one.matrix   <- matrix( rep( 1, n * tau.size ), ncol = tau.size )
	t.int.ind    <- t( int.ind )
	t.risk.ind   <- t( risk.ind )
	t.count.ind  <- t( count.ind )
	t.one.matrix <- t( one.matrix )
	cum.ind      <- matrix( rep( 1, tau.size^2 ), ncol = tau.size )
	cum.ind[ upper.tri( cum.ind, diag = FALSE ) ] <- 0

	#---------------#
	#--- beta_PH ---#
	#---------------#
	beta0             <- rep( 0, beta.size )
	beta_PH.res       <- beta_PH.function( n, beta.size, tau.size, delta.d, mat.Z.d, mat.Z.d.cro, int.ind, count.ind, risk.ind, t.risk.ind )

	logL_PH.NULL      <- beta_PH.res$outputNULL$logL
	Score_PH.NULL     <- beta_PH.res$outputNULL$Score
	Fisher_PH.NULL    <- beta_PH.res$outputNULL$Fisher
	Ewi2_PH.NULL      <- beta_PH.res$outputNULL$Ewi2

	beta_PH.Hat       <- beta_PH.res$outputHat$beta
	bz_PH             <- beta_PH.res$outputHat$bz
	exp.bz_PH         <- beta_PH.res$outputHat$exp.bz

	s0tZ_PH           <- beta_PH.res$outputHat$s0tZ
	s1tZ_PH           <- beta_PH.res$outputHat$s1tZ
	s2tZ_PH           <- beta_PH.res$outputHat$s2tZ
	logs0tZ_PH        <- beta_PH.res$outputHat$logs0tZ
	s1tZ_s0tZ_PH      <- beta_PH.res$outputHat$s1tZ_s0tZ
	s1tZ_s0tZ2_PH     <- beta_PH.res$outputHat$s1tZ_s0tZ2
	s2tZ_s0tZ_PH      <- beta_PH.res$outputHat$s2tZ_s0tZ

	logL_PH.Hat       <- beta_PH.res$outputHat$logL
	Score_PH.Hat      <- beta_PH.res$outputHat$Score
	Fisher_PH.Hat     <- beta_PH.res$outputHat$Fisher

	wi_PH             <- beta_PH.res$outputHat$wi
	Ewi2_PH           <- beta_PH.res$outputHat$Ewi2
	Var.beta_PH       <- beta_PH.res$outputHat$Var
	RobustVar.beta_PH <- beta_PH.res$outputHat$RobustVar

	hazard_PH         <- as.data.frame( apply( delta.d * int.ind, 2, sum )/s0tZ_PH )
	Hazard_PH         <- as.data.frame( cum.ind %*% as.matrix( hazard_PH ) )
	survival_PH       <- exp( -Hazard_PH )

	#--- for baseline cumulative hazard ---#

	Var.bch_PH        <- cum.ind %*% as.matrix( 
				 apply( count.ind, 2, mean )/( s0tZ_PH/n )^2 + diag( 
				  as.matrix( 
				   -s1tZ_s0tZ_PH/as.matrix( s0tZ_PH/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
				  ) %*% solve( Fisher_PH.Hat ) %*% t( 
				   as.matrix( 
				    -s1tZ_s0tZ_PH/as.matrix( s0tZ_PH/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
				   ) 
				  ) 
				 ) 
				)/n


	qid1_PH           <- t( 
				 ( 
				  cum.ind %*% as.matrix( 
				   - s1tZ_s0tZ_PH/as.matrix( s0tZ_PH/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
				  ) 
				 ) %*% t( wi_PH ) 
				)
	qid2_PH           <- t( 
					 cum.ind %*% as.matrix( 
					  as.data.frame( 
					   t( 
					    count.ind - t( 
					     as.data.frame( t( risk.ind * exp.bz_PH ) )/as.matrix( s0tZ_PH/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
					    ) 
					   ) 
					  )/as.matrix( s0tZ_PH/n ) 
					 ) 
					)

	qi_PH             <- qid1_PH + qid2_PH

	RobustVar.bch_PH  <- apply( qi_PH^2, 2, mean )/n

	#-----------------------------#
	#--- Goodness of fit tests ---#
	#-----------------------------#
		#--- likelihood ratio test ---#
		logL_PH           <- logL_PH.Hat
		logL_PH.s         <- - 2 * ( logL_PH.NULL - logL_PH )
		logL_PH.p         <- 1-round( pchisq( logL_PH.s, beta.size), 4 )
		LRStat_PH         <- data.frame( logL_PH.NULL, logL_PH ) ; names( LRStat_PH ) <- c( "NULL", "Optimal" ) ; rownames( LRStat_PH ) <- "LR"
		LRTest_PH         <- data.frame( logL_PH.s, logL_PH.p )  ; names( LRTest_PH ) <- c( "Stat", "p-value" ) ; rownames( LRTest_PH ) <- "LR"

		#--- Score test ---#
		Score_PH          <- Score_PH.Hat
		Score_PH.s        <- Score_PH.NULL %*% solve( Fisher_PH.NULL*n ) %*% matrix( Score_PH.NULL )
		Score_PH.p        <- 1-round( pchisq( Score_PH.s, beta.size), 4 )
		ScStat_PH         <- data.frame( Score_PH.NULL, Score_PH ) ; names( ScStat_PH ) <- c( "NULL", "Optimal" )
		ScTest_PH         <- data.frame( Score_PH.s, Score_PH.p )  ; names( ScTest_PH ) <- c( "Stat", "p-value" ) ; rownames( ScTest_PH ) <- "Score"

		#--- Robust Score test ---#
		Score_PH          <- Score_PH.Hat
		ScoreRobust_PH.s  <- Score_PH.NULL %*% solve( Ewi2_PH*n ) %*% matrix( Score_PH.NULL )
		ScoreRobust_PH.p  <- 1-round( pchisq( ScoreRobust_PH.s, beta.size ), 4 )
		ScRobustTest_PH   <- data.frame( ScoreRobust_PH.s, ScoreRobust_PH.p )  ; names( ScRobustTest_PH ) <- c( "Stat", "p-value" ) ; rownames( ScRobustTest_PH ) <- "ScRobust"

		#--- Wald test ---#
		Wald_PH           <- beta_PH.Hat
		Wald_PH.s         <- beta_PH.Hat %*% solve( Var.beta_PH ) %*% matrix( beta_PH.Hat )
		Wald_PH.p         <- 1-round( pchisq( Wald_PH.s,  beta.size), 4 )
		WdStat_PH         <- data.frame( beta0, Wald_PH ) ; names( WdStat_PH ) <- c( "NULL", "Optimal" )
		WdTest_PH         <- data.frame( Wald_PH.s, Wald_PH.p )  ; names( WdTest_PH ) <- c( "Stat", "p-value" ) ; rownames( WdTest_PH ) <- "Wald"

		#--- Robust Wald test ---#
		WaldRobust_PH.s   <- beta_PH.Hat %*% solve( RobustVar.beta_PH ) %*% matrix( beta_PH.Hat )
		WaldRobust_PH.p   <- 1-round( pchisq(WaldRobust_PH.s, beta.size ), 4 )
		WdRobustTest_PH   <- data.frame( WaldRobust_PH.s, WaldRobust_PH.p )  ; names( WdRobustTest_PH ) <- c( "Stat", "p-value" ) ; rownames( WdRobustTest_PH ) <- "WdRobust"

		beta_PH.s.p       <- cbind( 
					  logL_PH.s       , logL_PH.p       , 
					  Score_PH.s      , Score_PH.p      , 
					  ScoreRobust_PH.s, ScoreRobust_PH.p, 
					  Wald_PH.s       , Wald_PH.p       , 
					  WaldRobust_PH.s , WaldRobust_PH.p
					)
		names( beta_PH.s.p )    <- c(
						"Stat.logL"       , "p-value.logL", 
						"Stat.Score"      , "p-value.Score", 
						"Stat.ScoreRobust", "p-value.ScoreRobust",
						"Stat.Wald"       , "p-value.Wald", 
						"Stat.WaldRobust" , "p-value.WaldRobust"
						)

	GofTest_PH        <- rbind( LRTest_PH, ScTest_PH, ScRobustTest_PH, WdTest_PH, WdRobustTest_PH )


	#--------------------------#
	#--- covariate adjusted ---#
	#--------------------------#

	#-----------#
	#--- pai ---#
	#-----------#
	pai          <- mean( mat.Z.d )

	#--------------#
	#--- S2aHat ---#
	#--------------#
		#--- mHatDiZ ---#
		mHatDiZ_int  <- matrix( 
					 rep( mat.Z.d, tau.size ), ncol = tau.size 
					) - matrix( 
					 rep( as.matrix( s1tZ_s0tZ_PH ), rep( n, tau.size ) ), ncol = tau.size 
					)
		mHatDiZ_d    <- count.ind - risk.ind * matrix( 
					 rep( exp.bz_PH, tau.size ), ncol = tau.size 
					) * matrix( 
					 rep( as.matrix( hazard_PH ), rep( n, tau.size ) ), ncol = tau.size 
					)

		mHatDiZ      <- apply( mHatDiZ_int * mHatDiZ_d, 1, sum )

		#--- qVi ---#
		qVi          <- as.matrix( q.f( Covariate.d ) )

		#--- aHat ---#
		qVi.size     <- length( qVi[ 1, ] )
		aHat_d       <- solve( 
					 pai * ( 1-pai ) * matrix( 
					  apply( 
					   if( qVi.size == 1 ){ qVi^2 }else{ t( sapply( 1:n, function( x )qVi[ x, ] %x% qVi[ x, ] ) ) }, 2, sum 
					  ), ncol = qVi.size
					 )
					)
		aHat         <- aHat_d %*% apply(
					 as.data.frame( qVi ) * ( mat.Z.d - pai ) * mHatDiZ, 2, sum 
					)

		#--- f0ViaHat ---#
		f0ViaHat     <- qVi %*% aHat

	S2aHat       <- ( mat.Z.d - pai ) * f0ViaHat

	#--------------------------#
	#--- Estimation of beta ---#
	#--------------------------#
	beta0           <- rep( 0, beta.size )
	bz              <- as.numeric( mat.Z.d %*% beta0 )
	exp.bz          <- exp( bz )

	s0tZ            <- as.data.frame( t.risk.ind %*% as.matrix( exp.bz ) )
	s1tZ            <- as.data.frame( t.risk.ind %*% ( mat.Z.d * exp.bz ) )
	s2tZ            <- as.data.frame( t.risk.ind %*% as.matrix( mat.Z.d.cro * exp.bz ) )
	logs0tZ         <- as.data.frame( log( s0tZ ) )
	s1tZ_s0tZ       <- as.data.frame( s1tZ/as.matrix( s0tZ ) )
	s1tZ_s0tZ2      <- as.data.frame( if( beta.size == 1 ){ s1tZ_s0tZ^2 }else{ t( apply( s1tZ_s0tZ, 1, cro ) ) } )
	s2tZ_s0tZ       <- as.data.frame( s2tZ/as.matrix( s0tZ ) )

	logL.NULL       <- apply( delta.d * ( exp.bz - int.ind %*% as.matrix( logs0tZ ) ), 2, sum )
	ScoreTilde.NULL <- apply( delta.d * ( mat.Z.d - int.ind %*% as.matrix( s1tZ_s0tZ ) ) - S2aHat , 2, sum )
	Fisher.NULL     <- matrix( apply( delta.d * as.data.frame( int.ind %*% ( as.matrix( s2tZ_s0tZ - s1tZ_s0tZ2 ) ) ), 2, mean ), beta.size )

	wi_integrand    <- aperm( 
				 array( rep( mat.Z.d, tau.size ), dim = c( n, beta.size, tau.size ) ), perm = c( 1, 3, 2 ) 
				) - array( 
				 rep( as.matrix( s1tZ_s0tZ ), rep( n, tau.size * beta.size ) ), dim = c( n, tau.size, beta.size ) 
				)
	wi_d            <- array( 
				 rep( count.ind - 
				  t( 
				   as.data.frame( 
				    t( risk.ind * exp.bz ) 
				   )/as.matrix( s0tZ/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
				  ), beta.size 
				 ), dim = c( n, tau.size, beta.size ) 
				)
	wi              <- apply( wi_integrand * wi_d, c( 1, 3 ), sum ) - S2aHat


	Ewi2.NULL       <- if( beta.size == 1 ){ 
				 matrix( mean( wi^2 ) ) 
				}else{ 
				 (
				  matrix( 
				   apply( 
				    t( apply( wi, 1, cro ) ), 2, mean 
				   ), ncol = beta.size 
				  ) 
				 )
				}

	beta.new        <- as.numeric( beta0 + solve( Fisher.NULL * n ) %*% ScoreTilde.NULL )
	for(j in 1:100){
		beta.old     <- beta.new

		bz           <- as.numeric( mat.Z.d %*% beta.old )
		exp.bz       <- exp( bz )

		s0tZ         <- as.data.frame( t.risk.ind %*% as.matrix( exp.bz ) )
		s1tZ         <- as.data.frame( t.risk.ind %*% ( mat.Z.d * exp.bz ) )
		s2tZ         <- as.data.frame( t.risk.ind %*% as.matrix( mat.Z.d.cro * exp.bz ) )
		logs0tZ      <- as.data.frame( log( s0tZ ) )
		s1tZ_s0tZ    <- as.data.frame( s1tZ/as.matrix( s0tZ ) )
		s1tZ_s0tZ2   <- as.data.frame( if( beta.size == 1 ){ s1tZ_s0tZ^2 }else{ t( apply( s1tZ_s0tZ, 1, cro ) ) } )
		s2tZ_s0tZ    <- as.data.frame( s2tZ/as.matrix( s0tZ ) )

		#ScoreTilde   <- apply( delta.d * ( mat.Z.d - int.ind %*% as.matrix( s1tZ_s0tZ ) ) - S2aHat + S3bHat, 2, sum ) 
		ScoreTilde   <- apply( delta.d * ( mat.Z.d - int.ind %*% as.matrix( s1tZ_s0tZ ) ) - S2aHat , 2, sum )    
		Fisher       <- matrix( apply( delta.d * as.data.frame( int.ind %*% ( as.matrix( s2tZ_s0tZ - s1tZ_s0tZ2 ) ) ), 2, sum ), beta.size )/n

		beta.new     <- as.numeric( beta.old + solve( Fisher*n ) %*% ScoreTilde )
		if( sum( abs( c( beta.new - beta.old ) ) ) < 0.00001 )break
	}

	beta.Hat     <- beta.new
	bz           <- as.numeric( mat.Z.d %*% beta.Hat )
	exp.bz       <- exp( bz )

	s0tZ         <- as.data.frame( t.risk.ind %*% as.matrix( exp.bz ) )
	s1tZ         <- as.data.frame( t.risk.ind %*% ( mat.Z.d * exp.bz ) )
	s2tZ         <- as.data.frame( t.risk.ind %*% as.matrix( mat.Z.d.cro * exp.bz ) )
	logs0tZ      <- as.data.frame( log( s0tZ ) )
	s1tZ_s0tZ    <- as.data.frame( s1tZ/as.matrix( s0tZ ) )
	s1tZ_s0tZ2   <- as.data.frame( if( beta.size == 1 ){ s1tZ_s0tZ^2 }else{ t( apply( s1tZ_s0tZ, 1, cro ) ) } )
	s2tZ_s0tZ    <- as.data.frame( s2tZ/as.matrix( s0tZ ) )

	logL.Hat     <- apply( delta.d * ( exp.bz - int.ind %*% as.matrix( logs0tZ ) ), 2, sum )
	Score.Hat    <- apply( delta.d * ( mat.Z.d - int.ind %*% as.matrix( s1tZ_s0tZ ) ) - S2aHat , 2, sum )
	Fisher.Hat   <- matrix( apply( delta.d * as.data.frame( int.ind %*% ( as.matrix( s2tZ_s0tZ - s1tZ_s0tZ2 ) ) ), 2, mean ), beta.size )

	wi_integrand <- 	aperm( 
				 array( rep( mat.Z.d, tau.size ), dim = c( n, beta.size, tau.size ) ), perm = c( 1, 3, 2 ) 
				) - array( 
				 rep( as.matrix( s1tZ_s0tZ ), rep( n, tau.size * beta.size ) ), dim = c( n, tau.size, beta.size ) 
				)
	wi_d         <- 	array( 
				 rep( count.ind - 
				  t( 
				   as.data.frame( 
				    t( risk.ind * exp.bz ) 
				   )/as.matrix( s0tZ/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
				  ), beta.size 
				 ), dim = c( n, tau.size, beta.size ) 
				)
	wi           <- apply( wi_integrand * wi_d, c( 1, 3 ), sum ) - S2aHat 


	Ewi2         <- if( beta.size == 1 ){ 
				 matrix( mean( wi^2 ) ) 
				}else{ 
				 (
				  matrix( 
				   apply( 
				    t( apply( wi, 1, cro ) ), 2, mean 
				   ), ncol = beta.size 
				  ) 
				 )
				}

	Var.beta          <- solve( Fisher.Hat*n )                                   ## estimator
	RobustVar.beta    <- solve( Fisher.Hat ) %*% Ewi2 %*% solve( Fisher.Hat )/n  ## Robust estimator


	hazard            <- as.data.frame( apply( delta.d * int.ind, 2, sum )/s0tZ )
	Hazard            <- as.data.frame( cum.ind %*% as.matrix( hazard ) )
	survival          <- exp( -Hazard )

	#--- for baseline cumulative hazard ---#

	Var.bch           <- cum.ind %*% as.matrix( 
				 apply( count.ind, 2, mean )/( s0tZ/n )^2 + diag( 
				  as.matrix( 
				   -s1tZ_s0tZ/as.matrix( s0tZ/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
				  ) %*% solve( Fisher.Hat ) %*% t( 
				   as.matrix( 
				    -s1tZ_s0tZ/as.matrix( s0tZ/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
				   ) 
				  ) 
				 ) 
				)/n


	qid1              <- t( 
				 ( 
				  cum.ind %*% as.matrix( 
				   - s1tZ_s0tZ/as.matrix( s0tZ/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
				  ) 
				 ) %*% t( wi ) 
				)
	qid2              <- t( 
					 cum.ind %*% as.matrix( 
					  as.data.frame( 
					   t( 
					    count.ind - t( 
					     as.data.frame( t( risk.ind * exp.bz ) )/as.matrix( s0tZ/n ) * as.matrix( apply( count.ind, 2, mean ) ) 
					    ) 
					   ) 
					  )/as.matrix( s0tZ/n ) 
					 ) 
					)

	qi                <- qid1 + qid2

	RobustVar.bch     <- apply( qi^2, 2, mean )/n

	#-----------------------------#
	#--- Goodness of fit tests ---#
	#-----------------------------#
		#--- likelihood ratio test ---#
		logL              <- logL.Hat
		logL.s            <- - 2 * ( logL.NULL - logL )
		logL.p            <- 1-round( pchisq( logL.s, beta.size ), 4 )
		LRStat            <- data.frame( logL.NULL, logL ) ; names( LRStat ) <- c( "NULL", "Optimal" ) ; rownames( LRStat ) <- "LR"
		LRTest            <- data.frame( logL.s, logL.p )  ; names( LRTest ) <- c( "Stat", "p-value" ) ; rownames( LRTest ) <- "LR"

		#--- Score test ---#
		Score             <- Score.Hat
		Score.s           <- ScoreTilde.NULL %*% solve( Fisher.NULL*n ) %*% matrix( ScoreTilde.NULL )
		Score.p           <- 1-round( pchisq(Score.s, beta.size ), 4 )
		ScStat            <- data.frame( Score_PH.NULL, Score ) ; names( ScStat ) <- c( "NULL", "Optimal" )
		ScTest            <- data.frame( Score.s, Score.p )  ; names( ScTest ) <- c( "Stat", "p-value" ) ; rownames( ScTest ) <- "Score"

		#--- Robust Score test ---#
		Score             <- Score.Hat
		ScoreRobust.s     <- ScoreTilde.NULL %*% solve( Ewi2*n ) %*% matrix( ScoreTilde.NULL )
		ScoreRobust.p     <- 1-round( pchisq( ScoreRobust.s, beta.size ), 4 )
		ScRobustTest      <- data.frame( ScoreRobust.s, ScoreRobust.p )  ; names( ScRobustTest ) <- c( "Stat", "p-value" ) ; rownames( ScRobustTest ) <- "ScRobust"

		#--- Wald test ---#
		Wald              <- beta.Hat
		Wald.s            <- beta.Hat %*% solve( Var.beta ) %*% matrix( beta.Hat )
		Wald.p            <- 1-round( pchisq( Wald.s,  beta.size), 4 )
		WdStat            <- data.frame( beta0, Wald ) ; names( WdStat ) <- c( "NULL", "Optimal" )
		WdTest            <- data.frame( Wald.s, Wald.p )  ; names( WdTest ) <- c( "Stat", "p-value" ) ; rownames( WdTest ) <- "Wald"

		#--- Robust Wald test ---#
		WaldRobust.s      <- beta.Hat %*% solve( RobustVar.beta ) %*% matrix( beta.Hat )
		WaldRobust.p      <- 1-round( pchisq( WaldRobust.s, beta.size ), 4 )
		WdRobustTest      <- data.frame( WaldRobust.s, WaldRobust.p )  ; names( WdRobustTest ) <- c( "Stat", "p-value" ) ; rownames( WdRobustTest ) <- "WdRobust"

		beta.s.p          <- cbind( 
					  logL.s       , logL.p       , 
					  Score.s      , Score.p      , 
					  ScoreRobust.s, ScoreRobust.p, 
					  Wald.s       , Wald.p       , 
					  WaldRobust.s , WaldRobust.p
					)
		names( beta.s.p ) <- c(
						"Stat.logL"       , "p-value.logL", 
						"Stat.Score"      , "p-value.Score", 
						"Stat.ScoreRobust", "p-value.ScoreRobust",
						"Stat.Wald"       , "p-value.Wald", 
						"Stat.WaldRobust" , "p-value.WaldRobust"
						)

	GofTest           <- rbind( LRTest, ScTest, ScRobustTest, WdTest, WdRobustTest )


	#---------------#
	#--- Results ---#
	#---------------#

	#--- Unajusted ---#
		Results.beta                 <- cbind( beta_PH.Hat, sqrt( diag( Var.beta_PH ) ), sqrt( diag( RobustVar.beta_PH ) ), beta_PH.s.p )
		rownames( Results.beta )     <- names( Z )
		names( Results.beta )        <- c( "Est", "SE", "RobustSE", names( beta_PH.s.p ) )

		Results.survival             <- cbind( tau, survival_PH, sqrt( survival_PH^2 * Var.bch_PH ), sqrt( survival_PH^2 * RobustVar.bch_PH ) )
		names( Results.survival )    <- c( "Time", "Est", "SE", "RobustSE" )

		Results.Hazard               <- cbind( tau, Hazard_PH, sqrt( Var.bch_PH ), sqrt( RobustVar.bch_PH ) )
		names( Results.Hazard )      <- c( "Time", "Est", "SE", "RobustSE" )

		Results.Test                 <- GofTest_PH

		OtherResults                 <- list( Fisher_PH.Hat, wi_PH, qi_PH )
		names( OtherResults )        <- c( "FisherInformationMatrix", "wi", "qi" )

	Results.unajusted                  <- list( Results.beta,   Results.survival,   Results.Hazard,   Results.Test,   OtherResults )
	names( Results.unajusted )         <- c(   "Results.beta", "Results.survival", "Results.Hazard", "Results.Test", "OtherResults" )

	#--- Covariate-ajusted ---#
		Results.beta                 <- cbind( beta.Hat, sqrt( diag( Var.beta ) ), sqrt( diag( RobustVar.beta ) ), beta.s.p )
		rownames( Results.beta )     <- names( Z )
		names( Results.beta )        <- c( "Est", "SE", "RobustSE", names( beta.s.p ) )

		Results.survival             <- cbind( tau, survival, sqrt( survival^2 * Var.bch ), sqrt( survival^2 * RobustVar.bch ) )
		names( Results.survival )    <- c( "Time", "Est", "SE", "RobustSE" )

		Results.Hazard               <- cbind( tau, Hazard, sqrt( Var.bch ), sqrt( RobustVar.bch ) )
		names( Results.Hazard )      <- c( "Time", "Est", "SE", "RobustSE" )

		Results.Test                 <- GofTest

		OtherResults                 <- list( Fisher.Hat, wi, qi )
		names( OtherResults )        <- c( "FisherInformationMatrix", "wi", "qi" )

	Results.CovariateAjusted           <- list( Results.beta,   Results.survival,   Results.Hazard,   Results.Test,   OtherResults )
	names( Results.CovariateAjusted )  <- c(   "Results.beta", "Results.survival", "Results.Hazard", "Results.Test", "OtherResults" )

	Output                             <- list( Results.unajusted ,  Results.CovariateAjusted  )
	names( Output )                    <- c(   "Results.unajusted", "Results.CovariateAjusted" )

	return( Output )
}
############################################################################################

