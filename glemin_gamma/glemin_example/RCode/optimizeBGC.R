# Estimation of GC-biased gene conversion intensity from frequency spectrum data
# Sylvain Glemin
# sylvain.glemin@univ-montp2.fr



# source("basic_functions.R")
# source("likelihood_Poisson.R")
# source("likelihood_Error_Poisson.R")




###########################################
# Optimization of likelihood functions ####
###########################################

# Function taking SFSs and GC content in arguments and returning log-likelihood of the different models and parameters estimates
#
# Return
#		list of list:optimizeBGC$modelX$varY
#		for each model three kind of values are returned:
#		the log-likelihood of the model, the convergence statut, a list of parameter estimates (more or less according to the verbose option)
#
# Options
# 		BoundInf and BoundSup: lists of the boundaries to group SNP frequencies:
# 		  For i to length(BoundInf), frequencies between BoundInf[i] and BoundSup[i] are grouped together.
#			  Note that the user must correctly define the boundaries:
#			  ex BoundInf=c(1,2,5) and BoundInf=c(1,4,10) leads to three groups of frequencies: [1-1],[2-4], and [5-10]
#			  By default BoundInf=BoundSup=c(1:length(SFS)) so that every frequency are used separately
#			  If only BoundInf is given then BoundSup=BoundInf: frequencies are not grouped but some frequencies can be excluded.
#			  Example to exclude singleton: for 10 categories use BoundInf=c(2,9)
# 		models (boolean): to optimize model 2a, 2b. By default only model 0 and model 1 are optimized
#       If model 2b is chosen, model 2a is also optimized
# 		MAXIT: maximum number of iterations for the optim function.
#			  If the convergence message is not "OK", MAXIT can be increased.
#		  FACTR: 	precision factor. By default FACTR=10^7.
#				10^6 can be used if you think that likelihood is too low.
#				(for instance by comparison with the likelihood of the full model that can be obtained using multilnL)
#     BMAX  Constant to define the limits for searching B: ]-BMAX,+BMAX[ ==> to avoid out of limit problem especially for model2b and model3
#     All intinial parameters, except the ri factors, can be set up manually
#		  verbose: 0: basic output (default), 1: intermediate output, 2 (or any other value): full output. 
#		  PrintMessage: to print warning messages
#     PrintInit: to print initial values for optimization (except the Ri). Stored in modelX$init
#
# Cautions
#  	  The model 2b has two equivalent optimum by exchanging B0 and B1 and f and 1-f. To avoid the problem we fixed 0<f<1/2.
#     Because the search for B0 and B1 is limited by BMAX the likelihood of model 2b could be lower than for model 2a
#     If so, it means that the B0 or the B1 values are too extreme and likely not relevant biologically


optimizeBGC <- function(dataNeutral,dataWS,dataSW,dataGC,
                        BoundInf=c(1:length(dataNeutral)),BoundSup=BoundInf,models=c(FALSE,FALSE),
                        MAXIT=1000,FACTR=10^7,LMM=20,BMAX=100,
                        verbose=0,PrintMessage=TRUE,PrintInit=FALSE,
                        theta_init0=NULL,thetaWS_init0=NULL,bias_init0=NULL, #initial parameters for model0
                        theta_init1=NULL,thetaWS_init1=NULL,bias_init1=NULL,B_init1=NULL, #initial parameters for model1
                        theta_init2a=NULL,thetaWS_init2a=NULL,bias_init2a=NULL,B_init2a=NULL,f_init2a=NULL, #initial parameters of model2a
                        theta_init2b=NULL,thetaWS_init2b=NULL,bias_init2b=NULL,B0_init2b=NULL,B1_init2b=NULL,f_init2b=NULL #initial parameters of model2b
                        ) {

	n <- length(dataNeutral)
	ncat <- length(BoundInf)
  
  # Check the length of the vectors of data given to the function
  nws <- length(dataWS)
  nsw <- length(dataSW)
  if(n != nws || n != nsw) return( print("Failure: The frequency spectrum do not have compatible lengths!") ) 
  
  # Optimization of model 0
	if(is.null(theta_init0)) theta_init0 <- sum(dataNeutral)/HN(n+1)	#Watterson's estimate of theta is used as the starting value
	if(is.null(thetaWS_init0)) thetaWS_init0 <- sum(dataWS)/HN(n+1)
	if(is.null(bias_init0)) {thetaSW_init0 <- sum(dataSW)/HN(n+1) ; bias_init0 <- (dataGC*thetaSW_init0)/((1-dataGC)*thetaWS_init0)}
	vectR_init0 <- rep(1,ncat-1)
	Cnorm <- 0 # Normalisation constante for R_i parameters
	for(j in BoundInf[1]:BoundSup[1]){
		Cnorm <- Cnorm + (dataNeutral[j])*j/(BoundSup[1]-BoundInf[1]+1)
	}
	for(i in 2:ncat){
		Rinit <- 0
		for(j in BoundInf[i]:BoundSup[i]){
			Rinit <- Rinit + (dataNeutral[j])*j/(BoundSup[i]-BoundInf[i]+1)	
		}
		vectR_init0[i-1] <- Rinit/Cnorm	
	}
	init0 <- c(theta_init0,thetaWS_init0,bias_init0,vectR_init0)
	inf0 <- c(ZERO,ZERO,ZERO,rep(ZERO,ncat-1))
	sup0 <- c(Inf,Inf,Inf,rep(Inf,ncat-1))
	ml0 <- optim(init0,lnL0,SFSneutral=dataNeutral,SFSws=dataWS,SFSsw=dataSW,GC=dataGC,CatInf=BoundInf,CatSup=BoundSup,lower=inf0,upper=sup0,method="L-BFGS-B",control=list(maxit=MAXIT,factr=FACTR,lmm=LMM))
	warning <- FALSE
	# Output of model 0
	if(ml0$convergence==0) info <- "OK" else {info <- ml0$message; warning <- TRUE}
	if(verbose==0){
		out0 <- list(lnL=-ml0$value,convergence=info)
	} else {
		if(verbose==1){
			out0 <- list(lnL=-ml0$value,theta=ml0$par[1],thetaWS=ml0$par[2],bias=ml0$par[3],convergence=info)
		} else {
			out0 <- list(lnL=-ml0$value,theta=ml0$par[1],thetaWS=ml0$par[2],bias=ml0$par[3],vectR=ml0$par[4:(2+ncat)],convergence=info)
		}
	}
  if(PrintInit){
    param <- init0[1:3]
    names(param) <- c("theta","thetaWS","bias")
    out0 <-c(out0,init=list(param))
  }
	
	# Optimization of model 1
	if(is.null(theta_init1)) theta_init1 <- ml0$par[1]	# Estimates of the previous model are used as starting values
	if(is.null(thetaWS_init1)) thetaWS_init1 <- ml0$par[2]
	x<-c(1:n)/(n+1)
	y<-log((dataWS+1)/(dataSW+1)) # 1 is added to avoid division by 0 and log(0)
	if(is.null(B_init1)) B_init1 <- lm(y~x)$coefficients[[2]] 	# For the starting value, we used the property: log(SFSws/SFSsw) = B*j/n in stable populations
	if(is.null(bias_init1)) bias_init1 <- ml0$par[3]*exp(B_init1/5) # We used the fact that the bias is underestimated without gBGC, so that initial value is slightly increased compared to model 0
	vectR_init1 <- vectR_init0 # ml0$par[4:(2+ncat)]
	init1 <- c(theta_init1,thetaWS_init1,bias_init1,B_init1,vectR_init1)
	inf1 <- c(ZERO,ZERO,ZERO,-BMAX,rep(ZERO,ncat-1))
	sup1 <- c(Inf,Inf,Inf,BMAX,rep(Inf,ncat-1))
	ml1 <- optim(init1,lnL1,SFSneutral=dataNeutral,SFSws=dataWS,SFSsw=dataSW,GC=dataGC,CatInf=BoundInf,CatSup=BoundSup,lower=inf1,upper=sup1,method="L-BFGS-B",control=list(maxit=MAXIT,factr=FACTR,lmm=LMM))
	# Output of model 1
	if(ml1$convergence==0) info <- "OK" else {info <- ml1$message; warning <- TRUE}
	if(verbose==0){
		out1 <- list(lnL=-ml1$value,B=ml1$par[4],convergence=info)
	} else {
		if(verbose==1){
			out1 <- list(lnL=-ml1$value,theta=ml1$par[1],thetaWS=ml1$par[2],bias=ml1$par[3],B=ml1$par[4],convergence=info)
			} else {
				out1 <- list(lnL=-ml1$value,theta=ml1$par[1],thetaWS=ml1$par[2],bias=ml1$par[3],B=ml1$par[4],vectR=ml1$par[5:(3+ncat)],convergence=info)
		}
	}
	if(PrintInit){
	  param <- init1[1:4]
	  names(param) <- c("theta","thetaWS","bias","B")
	  out1 <-c(out1,init=list(param))
	}  
  
  
	# Optional models
	
	# Optimization of model 2a
	if (models[1] || models[2]){
	  if(is.null(theta_init2a)) theta_init2a <- ml1$par[1] # Estimates of the previous model are used as starting values
	  if(is.null(thetaWS_init2a)) thetaWS_init2a <- ml1$par[2]
	  if(is.null(bias_init2a)) bias_init2a <- ml1$par[3]
	  if(is.null(f_init2a)) f_init2a <- 0.25
	  if(is.null(B_init2a)) B_init2a <- ml1$par[4]/f_init2a
		vectR_init2a <- ml1$par[5:(3+ncat)]
		init2a <- c(theta_init2a,thetaWS_init2a,bias_init2a,B_init2a,f_init2a,vectR_init2a)
		inf2a <- c(ZERO,ZERO,ZERO,-BMAX,ZERO,rep(ZERO,ncat-1))
		sup2a <- c(Inf,Inf,Inf,BMAX,1,rep(Inf,ncat-1))
		ml2a <- optim(init2a,lnL2a,SFSneutral=dataNeutral,SFSws=dataWS,SFSsw=dataSW,GC=dataGC,CatInf=BoundInf,CatSup=BoundSup,lower=inf2a,upper=sup2a,method="L-BFGS-B",control=list(maxit=MAXIT,factr=FACTR,lmm=LMM))
		# Output of model 2a
		if(ml2a$convergence==0) info <- "OK" else {info <- ml2a$message; warning <- TRUE}
		if(verbose==0){
			out2a <- list(lnL=-ml2a$value,B=ml2a$par[4],f=ml2a$par[5],convergence=info)
		} else {
			if(verbose==1){
				out2a <- list(lnL=-ml2a$value,theta=ml2a$par[1],thetaWS=ml2a$par[2],bias=ml2a$par[3],B=ml2a$par[4],f=ml2a$par[5],convergence=info)
				} else {
					out2a <- list(lnL=-ml2a$value,theta=ml2a$par[1],thetaWS=ml2a$par[2],bias=ml2a$par[3],B=ml2a$par[4],f=ml2a$par[5],vectR=ml2a$par[6:(4+ncat)],convergence=info)
			}
		}
	  if(PrintInit){
	    param <- init2a[1:5]
	    names(param) <- c("theta","thetaWS","bias","B","f")
	    out2a <-c(out2a,init=list(param))
	  }
	}
 
	
	# Optimization of model 2b
	if (models[2]){
		if(is.null(theta_init2b)) theta_init2b <- ml2a$par[1] # Estimates of the previous model are used as starting values
		if(is.null(thetaWS_init2b)) thetaWS_init2b <- ml2a$par[2]
		if(is.null(bias_init2b)) bias_init2b <- ml2a$par[3]
		if(is.null(B0_init2b)) B0_init2b <- 0.1
		if(is.null(B1_init2b)) B1_init2b <- ml2a$par[4]
		if(is.null(f_init2b)) f_init2b <- 0.9*ml2a$par[5]
		vectR_init2b <- ml2a$par[6:(4+ncat)]
		init2b <- c(theta_init2b,thetaWS_init2b,bias_init2b,B0_init2b,B1_init2b,f_init2b,vectR_init2b)
		inf2b <- c(ZERO,ZERO,ZERO,-BMAX,-BMAX,ZERO,rep(ZERO,ncat-1))
		sup2b <- c(Inf,Inf,Inf,BMAX,BMAX,1/2,rep(Inf,ncat-1))
		ml2b <- optim(init2b,lnL2b,SFSneutral=dataNeutral,SFSws=dataWS,SFSsw=dataSW,GC=dataGC,CatInf=BoundInf,CatSup=BoundSup,lower=inf2b,upper=sup2b,method="L-BFGS-B",control=list(maxit=MAXIT,factr=FACTR,lmm=LMM))
		# Output of model 2b
		if(ml2b$convergence==0) info <- "OK" else {info <- ml2b$message; warning <- TRUE}
		if(verbose==0){
			out2b <- list(lnL=-ml2b$value,B0=ml2b$par[4],B1=ml2b$par[5],f=ml2b$par[6],convergence=info)
		} else {
			if(verbose==1){
				out2b <- list(lnL=-ml2b$value,theta=ml2b$par[1],thetaWS=ml2b$par[2],bias=ml2b$par[3],B0=ml2b$par[4],B1=ml2b$par[5],f=ml2b$par[6],convergence=info)
				} else {
					out2b <- list(lnL=-ml2b$value,theta=ml2b$par[1],thetaWS=ml2b$par[2],bias=ml2b$par[3],B0=ml2b$par[4],B1=ml2b$par[5],f=ml2b$par[6],vectR=ml2b$par[7:(5+ncat)],convergence=info)
			}
		}
		if(PrintInit){
		  param <- init2b[1:6]
		  names(param) <- c("theta","thetaWS","bias","B0","B1","f")
		  out2b <-c(out2b,init=list(param))
		}
	}
 
	
	if(PrintMessage & warning) print("One optimization (or more) may have failed to converge. Try increasing MAXIT or starting with other initial values")
	
	output <- list(model0=out0,model1=out1)
	if(models[1] || models[2]) output <-c(output,list(model2a=out2a))
	if(models[2]) output <-c(output,list(model2b=out2b))

	return(output)
	
}




##############################################
# Same functions with polarization errors ####
##############################################



optimizeBGCwith3error <- function(dataNeutral,dataWS,dataSW,dataGC,
                                  BoundInf=c(1:length(dataNeutral)),BoundSup=BoundInf,models=c(FALSE,FALSE),
                                  MAXIT=1000,FACTR=10^7,LMM=20,BMAX=100,
                                  verbose=0,PrintMessage=TRUE,PrintInit=FALSE,
                                  theta_init0=NULL,thetaWS_init0=NULL,bias_init0=NULL,err_init0=NULL,errWS_init0=NULL,errSW_init0=NULL, #initial parameters for model0
                                  theta_init1=NULL,thetaWS_init1=NULL,bias_init1=NULL,B_init1=NULL,err_init1=NULL,errWS_init1=NULL,errSW_init1=NULL, #initial parameters for model1
                                  theta_init2a=NULL,thetaWS_init2a=NULL,bias_init2a=NULL,B_init2a=NULL,f_init2a=NULL,err_init2a=NULL,errWS_init2a=NULL,errSW_init2a=NULL, #initial parameters of model2a
                                  theta_init2b=NULL,thetaWS_init2b=NULL,bias_init2b=NULL,B0_init2b=NULL,B1_init2b=NULL,f_init2b=NULL,err_init2b=NULL,errWS_init2b=NULL,errSW_init2b=NULL #initial parameters of model2b     
) {
  
  n <- length(dataNeutral)
  # Check the length of the vectors of data given to the function
  nws <- length(dataWS)
  nsw <- length(dataSW)
  if(n != nws || n != nsw) return( print("Failure: The frequency spectrum do not have compatible lengths!") ) 
  
  # Symetrize the vectors of bounadries
  BoundInf <- symetrize(BoundInf,BoundSup,n+1)$CinfSym
  BoundSup <- symetrize(BoundInf,BoundSup,n+1)$CsupSym
  
  ncat <- length(BoundInf)
  
  
  # Optimization of model 0
  if(is.null(theta_init0)) theta_init0 <- sum(dataNeutral)/HN(n+1)	#Watterson's estimate of theta is used as the starting value
  if(is.null(thetaWS_init0)) thetaWS_init0 <- sum(dataWS)/HN(n+1)
  if(is.null(bias_init0)) {thetaSW_init0 <- sum(dataSW)/HN(n+1); bias_init0 <- (dataGC*thetaSW_init0)/((1-dataGC)*thetaWS_init0)}
  if(is.null(err_init0)) err_init0 <- c(0.01) # 1% of orientation error for starting values
  if(is.null(errWS_init0)) errWS_init0 <- c(0.01) # 1% of orientation error for starting values
  if(is.null(errSW_init0)) errSW_init0 <- c(0.01) # 1% of orientation error for starting values
  vectR_init0 <- rep(1,ncat-1)
  Cnorm <- 0 # Normalisation constante for R_i parameters
  for(j in BoundInf[1]:BoundSup[1]){
    Cnorm <- Cnorm + (dataNeutral[j])*j/(BoundSup[1]-BoundInf[1]+1)
  }
  for(i in 2:ncat){
    Rinit <- 0
    for(j in BoundInf[i]:BoundSup[i]){
      Rinit <- Rinit + (dataNeutral[j])*j/(BoundSup[i]-BoundInf[i]+1)	
    }
    vectR_init0[i-1] <- Rinit/Cnorm	
  }
  init0 <- c(theta_init0,thetaWS_init0,bias_init0,vectR_init0,err_init0,errWS_init0,errSW_init0)
  inf0 <- c(ZERO,ZERO,ZERO,rep(ZERO,ncat-1),ZERO,ZERO,ZERO)
  sup0 <- c(Inf,Inf,Inf,rep(Inf,ncat-1),1/2,1/2,1/2)
  ml0 <- optim(init0,lnL0err3,SFSneutral=dataNeutral,SFSws=dataWS,SFSsw=dataSW,GC=dataGC,CatInf=BoundInf,CatSup=BoundSup,lower=inf0,upper=sup0,method="L-BFGS-B",control=list(maxit=MAXIT,factr=FACTR,lmm=LMM))
  warning <- FALSE
  # Output of model 0
  if(ml0$convergence==0) info <- "OK" else {info <- ml0$message; warning <- TRUE}
  if(verbose==0){
    out0 <- list(lnL=-ml0$value,convergence=info)
  } else {
    if(verbose==1){
      out0 <- list(lnL=-ml0$value,theta=ml0$par[1],thetaWS=ml0$par[2],bias=ml0$par[3],err=ml0$par[ncat+3],errWS=ml0$par[ncat+4],errSW=ml0$par[ncat+5],convergence=info)
    } else {
      out0 <- list(lnL=-ml0$value,theta=ml0$par[1],thetaWS=ml0$par[2],bias=ml0$par[3],err=ml0$par[ncat+3],errWS=ml0$par[ncat+4],errSW=ml0$par[ncat+5],vectR=ml0$par[4:(2+ncat)],convergence=info)
    }
  }
  if(PrintInit){
    param <- init0[-4:-(ncat+2)]
    names(param) <- c("theta","thetaWS","bias","err","errWS","errSW")
    out0 <-c(out0,init=list(param))
  }
  
  # Optimization of model 1
  if(is.null(theta_init1)) theta_init1 <- ml0$par[1]	# Estimates of the previous model are used as starting values
  if(is.null(thetaWS_init1)) thetaWS_init1 <- ml0$par[2]
  if(is.null(B_init1)){
    x<-c(1:n)/(n+1)
    y<-log((dataWS+1)/(dataSW+1)) # 1 is added to avoid division by 0 and log(0)
    B_init1 <- lm(y~x)$coefficients[[2]] 	# For the starting value, we used the property: log(SFSws/SFSsw) = B*j/n in stable populations
  }
  if(is.null(bias_init1))  bias_init1 <- ml0$par[3]*exp(B_init1/5) # We used the fact that the bias is underestimated without gBGC, so that initial value is slightly increased compared to model 0
  if(is.null(err_init1))   err_init1 <- ml0$par[ncat+3]
  if(is.null(errWS_init1))   errWS_init1 <- ml0$par[ncat+4]
  if(is.null(errSW_init1))   errSW_init1 <- ml0$par[ncat+5]
  vectR_init1 <- vectR_init0 # ml0$par[4:(2+ncat)]
  init1 <- c(theta_init1,thetaWS_init1,bias_init1,B_init1,vectR_init1,err_init1,errWS_init1,errSW_init1)
  inf1 <- c(ZERO,ZERO,ZERO,-BMAX,rep(ZERO,ncat-1),ZERO,ZERO,ZERO)
  sup1 <- c(Inf,Inf,Inf,BMAX,rep(Inf,ncat-1),1/2,1/2,1/2)
  ml1 <- optim(init1,lnL1err3,SFSneutral=dataNeutral,SFSws=dataWS,SFSsw=dataSW,GC=dataGC,CatInf=BoundInf,CatSup=BoundSup,lower=inf1,upper=sup1,method="L-BFGS-B",control=list(maxit=MAXIT,factr=FACTR,lmm=LMM))
  # Output of model 1
  if(ml1$convergence==0) info <- "OK" else {info <- ml1$message; warning <- TRUE}
  if(verbose==0){
    out1 <- list(lnL=-ml1$value,B=ml1$par[4],convergence=info)
  } else {
    if(verbose==1){
      out1 <- list(lnL=-ml1$value,theta=ml1$par[1],thetaWS=ml1$par[2],bias=ml1$par[3],B=ml1$par[4],err=ml1$par[ncat+4],errWS=ml1$par[ncat+5],errSW=ml1$par[ncat+6],convergence=info)
    } else {
      out1 <- list(lnL=-ml1$value,theta=ml1$par[1],thetaWS=ml1$par[2],bias=ml1$par[3],B=ml1$par[4],err=ml1$par[ncat+4],errWS=ml1$par[ncat+5],errSW=ml1$par[ncat+6],vectR=ml1$par[5:(3+ncat)],convergence=info)
    }
  }
  if(PrintInit){
    param <- init1[-5:-(ncat+3)]
    names(param) <- c("theta","thetaWS","bias","B","err","errWS","errSW")
    out1 <-c(out1,init=list(param))
  }
  
  
  # Optional models
  
  # Optimization of model 2a
  if (models[1] || models[2]){
    if(is.null(theta_init2a)) theta_init2a <- ml1$par[1] # Estimates of the previous model are used as starting values
    if(is.null(thetaWS_init2a)) thetaWS_init2a <- ml1$par[2]
    if(is.null(bias_init2a)) bias_init2a <- ml1$par[3]
    if(is.null(f_init2a)) f_init2a <- 0.9
    if(is.null(B_init2a)) B_init2a <- ml1$par[4]/f_init2a
    if(is.null(err_init2a)) err_init2a <- ml1$par[ncat+4]
    if(is.null(errWS_init2a)) errWS_init2a <- ml1$par[ncat+5]
    if(is.null(errSW_init2a)) errSW_init2a <- ml1$par[ncat+6]
    vectR_init2a <- ml1$par[5:(3+ncat)]
    init2a <- c(theta_init2a,thetaWS_init2a,bias_init2a,B_init2a,f_init2a,vectR_init2a,err_init2a,errWS_init2a,errSW_init2a)
    inf2a <- c(ZERO,ZERO,ZERO,-BMAX,ZERO,rep(ZERO,ncat-1),ZERO,ZERO,ZERO)
    sup2a <- c(Inf,Inf,Inf,BMAX,1,rep(Inf,ncat-1),1/2,1/2,1/2)
    ml2a <- optim(init2a,lnL2aerr3,SFSneutral=dataNeutral,SFSws=dataWS,SFSsw=dataSW,GC=dataGC,CatInf=BoundInf,CatSup=BoundSup,lower=inf2a,upper=sup2a,method="L-BFGS-B",control=list(maxit=MAXIT,factr=FACTR,lmm=LMM))
    # Output of model 2a
    if(ml2a$convergence==0) info <- "OK" else {info <- ml2a$message; warning <- TRUE}
    if(verbose==0){
      out2a <- list(lnL=-ml2a$value,B=ml2a$par[4],f=ml2a$par[5],convergence=info)
    } else {
      if(verbose==1){
        out2a <- list(lnL=-ml2a$value,theta=ml2a$par[1],thetaWS=ml2a$par[2],bias=ml2a$par[3],B=ml2a$par[4],f=ml2a$par[5],err=ml2a$par[ncat+5],errWS=ml2a$par[ncat+6],errSW=ml2a$par[ncat+7],convergence=info)
      } else {
        out2a <- list(lnL=-ml2a$value,theta=ml2a$par[1],thetaWS=ml2a$par[2],bias=ml2a$par[3],B=ml2a$par[4],f=ml2a$par[5],err=ml2a$par[ncat+5],errWS=ml2a$par[ncat+6],errSW=ml2a$par[ncat+7],vectR=ml2a$par[6:(4+ncat)],convergence=info)
      }
    }
    if(PrintInit){
      param <- init2a[-6:-(ncat+4)]
      names(param) <- c("theta","thetaWS","bias","B","f","err","errWS","errSW")
      out2a <-c(out2a,init=list(param))
    }
  }
  
  # Optimization of model 2b
  if (models[2]){
    if(is.null(theta_init2b)) theta_init2b <- ml2a$par[1] # Estimates of the previous model are used as starting values
    if(is.null(thetaWS_init2b)) thetaWS_init2b <- ml2a$par[2]
    if(is.null(bias_init2b)) bias_init2b <- ml2a$par[3]
    if(is.null(B0_init2b)) B0_init2b <- ml2a$par[4]*0.9
    if(is.null(B1_init2b)) B1_init2b <- ml2a$par[4]*1.1
    if(is.null(f_init2b)) f_init2b <- 0.9*ml2a$par[5]/2
    if(is.null(err_init2b)) err_init2b <- ml2a$par[ncat+5]
    if(is.null(errWS_init2b)) errWS_init2b <- ml2a$par[ncat+6]
    if(is.null(errSW_init2b)) errSW_init2b <- ml2a$par[ncat+7]
    vectR_init2b <- ml2a$par[6:(4+ncat)]
    init2b <- c(theta_init2b,thetaWS_init2b,bias_init2b,B0_init2b,B1_init2b,f_init2b,vectR_init2b,err_init2b,errWS_init2b,errSW_init2b)
    inf2b <- c(ZERO,ZERO,ZERO,-BMAX,-BMAX,ZERO,rep(ZERO,ncat-1),ZERO,ZERO,ZERO)
    sup2b <- c(Inf,Inf,Inf,BMAX,BMAX,1/2,rep(Inf,ncat-1),1/2,1/2,1/2)
    ml2b <- optim(init2b,lnL2berr3,SFSneutral=dataNeutral,SFSws=dataWS,SFSsw=dataSW,GC=dataGC,CatInf=BoundInf,CatSup=BoundSup,lower=inf2b,upper=sup2b,method="L-BFGS-B",control=list(maxit=MAXIT,factr=FACTR,lmm=LMM))
    # Output of model 2b
    if(ml2b$convergence==0) info <- "OK" else {info <- ml2b$message; warning <- TRUE}
    if(verbose==0){
      out2b <- list(lnL=-ml2b$value,B0=ml2b$par[4],B1=ml2b$par[5],f=ml2b$par[6],convergence=info)
    } else {
      if(verbose==1){
        out2b <- list(lnL=-ml2b$value,theta=ml2b$par[1],thetaWS=ml2b$par[2],bias=ml2b$par[3],B0=ml2b$par[4],B1=ml2b$par[5],f=ml2b$par[6],err=ml2b$par[ncat+6],errWS=ml2b$par[ncat+7],errSW=ml2b$par[ncat+8],convergence=info)
      } else {
        out2b <- list(lnL=-ml2b$value,theta=ml2b$par[1],thetaWS=ml2b$par[2],bias=ml2b$par[3],B0=ml2b$par[4],B1=ml2b$par[5],f=ml2b$par[6],err=ml2b$par[ncat+6],errWS=ml2b$par[ncat+7],errSW=ml2b$par[ncat+8],vectR=ml2b$par[7:(5+ncat)],convergence=info)
      }
    }
    if(PrintInit){
      param <- init2b[-7:-(ncat+5)]
      names(param) <- c("theta","thetaWS","bias","B0","B1","f","err","errWS","errSW")
      out2b <-c(out2b,init=list(param))
    }
  }
  
  
  if(PrintMessage & warning) print("One optimization (or more) may have failed to converge. Try increasing MAXIT or starting with other initial values")
  
  output <- list(model0=out0,model1=out1)
  if(models[1] || models[2]) output <-c(output,list(model2a=out2a))
  if(models[2]) output <-c(output,list(model2b=out2b))
  
  return(output)
  
}
