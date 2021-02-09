# Estimation of GC-biased gene conversion intensity from frequency spectrum data
# Sylvain Glémin




####################
# Useful functions #
####################

library("HMP") #library needed to use the dirichlet-multinomial distribution

library("VGAM")  #library needed for the use of the lerch function
# We use the property: Lerch(1,s,v) = HurwitzZeta(s,v)
# Because the continuity in 1 of the Lerch function is not implemented in R we used ONE = 1-10^(-10)
ZERO <- 10^(-10) # Expression used when functions are not defined in 0
ONE <- 1-ZERO


# UNFOLDED SPECTRUM ####

# Several functions giving the expected Number of SNP (ENS) in frequency j/n under different models


# Neutrality
ENS_neutral <- function(theta,j) {
	theta/j
}

# Constant gBGC with intensity B. B can be positive (W->S mutations) or negative (S->W) mutations
ENS_constant <- function(theta,B,n,j) {
	if (B==0) return(ENS_neutral(theta,j))
	((n*theta)/(j*(n-j))) * ((1-exp(-B*(1-j/n)))/(1-exp(-B)))
}

# Hotspot1: gBGC with intensity B in hotspot (fraction f) and 0 otherwise (fraction 1-f)
ENS_hotspot1 <- function(theta,B,f,n,j) {
	if (B==0) return(ENS_neutral(theta,j))
	f * ENS_constant(theta,B,n,j) + (1-f) * ENS_neutral(theta,j)
}

# Hotspot2: gBGC with intensity B1 in hotspot (fraction f) and B0 otherwise (fraction 1-f)
ENS_hotspot2 <- function(theta,B0,B1,f,n,j) {
	if (B0==0) return(ENS_hotpsot1(theta,B1,f,n,j))
	if (B1==0) return(ENS_hotpsot1(theta,B0,1-f,n,j))
	if (B0==0 & B1==0) return(ENS_neutral(theta,j))
	f * ENS_constant(theta,B1,n,j) + (1-f) * ENS_constant(theta,B0,n,j)
}




# UNFOLDED SPECTRUM WITH ORIENTATION ERROR ####

# Neutrality
ENS_neutral_err <- function(theta1,theta2,e1,e2,n,j) {
  (1-e1)*ENS_neutral(theta1,j) + e2*ENS_neutral(theta2,n-j)
}

# Constant gBGC with intensity B. B can be positive (W->S mutations) or negative (S->W) mutations
ENS_constant_err <- function(theta1,theta2,B,e1,e2,n,j) {
  (1-e1)*ENS_constant(theta1,B,n,j) + e2*ENS_constant(theta2,-B,n,n-j)
}

# Hotspot1: gBGC with intensity B in hotspot (fraction f) and 0 otherwise (fraction 1-f)
ENS_hotspot1_err <- function(theta1,theta2,B,f,e1,e2,n,j) {
  (1-e1)*ENS_hotspot1(theta1,B,f,n,j) + e2*ENS_hotspot1(theta2,-B,f,n,n-j)
}

# Hotspot2: gBGC with intensity B1 in hotspot (fraction f) and B0 otherwise (fraction 1-f)
ENS_hotspot2_err <- function(theta1,theta2,B0,B1,f,e1,e2,n,j) {
  (1-e1)*ENS_hotspot2(theta1,B0,B1,f,n,j) + e2*ENS_hotspot2(theta2,-B0,-B1,f,n,n-j)
}


# Gamma distribution of gBGC intensity: mean: Bmean and shape: beta.
# Here B>0 and B<0 must be written in two separate functions and a positive B must be sent to each function
# The log of the function is used to avoid too large number
ENS_gamma_pos_err <- function(theta1,theta2,Bmean,beta,e1,e2,n,j) {
  (1-e1)*ENS_gamma_pos(theta1,Bmean,beta,n,j) + e2*ENS_gamma_neg(theta2,Bmean,beta,n,n-j)
}

ENS_gamma_neg_err <- function(theta1,theta2,Bmean,beta,e1,e2,n,j) {
  (1-e1)*ENS_gamma_neg(theta1,Bmean,beta,n,j) + e2*ENS_gamma_pos(theta2,Bmean,beta,n,n-j)
}



# OTHER FUNCTIONS ####

# Log Poisson
# We redefine the expression instead of using dpois(x,lambda,log=TRUE) that works only for interger values
# We also make the continuous extension of the function in (0,0)
log_Poisson <- function(k,lambda){
  if(lambda==0 & k==0) return(0)
	-lambda + k*log(lambda) - lfactorial(k)
}


# Harmonic number
HN <- function(n) sum(1/c(1:n))

# Multinomial log-likelihood
multilnL <- function(vect,BoundInf=c(1:length(vect)),BoundSup=BoundInf){
	if(identical(BoundInf,c(1:length(vect)))){
		N <- sum(vect)
		lnL <- lfactorial(N) + sum(vect*log(vect/N)-lfactorial(vect))
	} else {
		ncat <- length(BoundInf)
		vectcat <- rep(0,ncat)
		for(i in 1:ncat){
			for(j in BoundInf[i]:BoundSup[i]){
				vectcat[i] <- vectcat[i] + vect[j]
			}		
		}
		N <- sum(vectcat)
		lnL <- lfactorial(N) + sum(vectcat*log(vectcat/N)-lfactorial(vectcat))
	}
return(lnL)
}


# Likelihood of the saturated Poisson model
saturatedlnL <- function(vect,BoundInf=c(1:length(vect)),BoundSup=BoundInf){
  lnL <- 0
  if(identical(BoundInf,c(1:length(vect)))){
    for (i in 1:length(vect)) lnL <- lnL + log_Poisson(vect[i],vect[i])
  } else {
    ncat <- length(BoundInf)
    vectcat <- rep(0,ncat)
    for(i in 1:ncat){
      for(j in BoundInf[i]:BoundSup[i]){
        vectcat[i] <- vectcat[i] + vect[j]
      }	
      lnL <- lnL + log_Poisson(vectcat[i],vectcat[i])
    }
  }
  return(lnL)
}



# Function to compute the skweness of a frequency spectrum and testing the asymetry
# Typically for testing the asymetry of the GC spectrum

spectrum_skewness <- function(sfs) {
  Nclass <- length(sfs)+1
  Nobs <- sum(sfs)
  freq <- c(1:(Nclass-1))/Nclass
  mean <- sum(sfs*freq)/Nobs
  skew <- (sqrt(Nobs*(Nobs-1))/(Nobs-2))*sum(sfs*(freq-mean)^3/Nobs)/(sum(sfs*(freq-mean)^2)/Nobs)^3/2
  SES <- sqrt(6*Nobs*(Nobs-1)/((Nobs-2)*(Nobs+1)*(Nobs+3))) # Standard error of skewness
  pval <- 2*pnorm(abs(skew/SES),0,1,lower.tail=F)
  return(list(skewness=skew,p.value=pval))
}


# Fonction to reduce the SFS (as a vector) from n to m < n chromosomes
# Return the reduced SFS vector


reduceSFS <- function(sfs,m){
	n <- length(sfs)
	output <- rep(0,m)
	index<-c(1:m)
	for(k in 1:n){
		output <- output + sfs[k]*choose(k,index)*choose(n-k,m-index)/choose(n,m)
	}
	return(output)
}

# Function to symetrize the vector of boundaries

symetrize <- function(CatInf,CatSup,Nchrom){
  Cinf <- union(CatInf[CatInf<=floor(Nchrom/2)],Nchrom-rev(CatSup[CatSup<=floor(Nchrom/2)]))
  Csup <- union(CatSup[CatSup<=floor(Nchrom/2)],Nchrom-rev(CatInf[CatInf<=floor(Nchrom/2)]))  
  return(list(CinfSym=Cinf,CsupSym=Csup))
}

