
########################
# Likelihood functions #
########################

# Poisson distribution
######################

# Estimation of the orientation error is included


# source("/Users/sylvain/Bin/R/estimBGC/basic_functions.R")



# Null model: no gBGC
# theta,thetaWS,bias,vectR,errN,errWS,errSW: parameters to be optimized
# SFSneutral,SFSws,SFSsw,GC: data
# CatInf and CatSup: optional vectors to chose the SNPs and the grouping of SNPs to be used
# However, because of the symetry of the error scheme, groups must be symetric in i / n-i. This constrain is included in the function

lnL0err3 <- function(par,SFSneutral,SFSws,SFSsw,GC,CatInf=c(1:length(SFSneutral)),CatSup=CatInf){
  # Constraining categories to be symetric
  Nchrom <- length(SFSneutral)+1
  Cinf <- union(CatInf[CatInf<=floor(Nchrom/2)],Nchrom-rev(CatSup[CatSup<=floor(Nchrom/2)]))
  Csup <- union(CatSup[CatSup<=floor(Nchrom/2)],Nchrom-rev(CatInf[CatInf<=floor(Nchrom/2)]))
  Ncat <- length(Cinf)
	theta <- par[1]
	thetaWS <- par[2]
	bias <- par[3]
	vectR <- c(1,par[4:(2+Ncat)])
  errN <- par[Ncat+3]
	errWS <- par[Ncat+4]
	errSW <- par[Ncat+5] 
	if(identical(Cinf,Csup)){
	  sfsN <- SFSneutral[CatInf]
	  sfsWS <- SFSws[CatInf]
	  sfsSW <- SFSsw[CatInf]
	  lnLN <- log_Poisson(sfsN,ENS_neutral_err(theta*vectR,theta*rev(vectR),errN,errN,Nchrom,CatInf))
	  lnLWS <- log_Poisson(sfsWS,ENS_neutral_err((1-GC)*thetaWS*vectR,GC*bias*thetaWS*rev(vectR),errWS,errSW,Nchrom,CatInf))
	  lnLSW <- log_Poisson(sfsSW,ENS_neutral_err(GC*bias*thetaWS*vectR,(1-GC)*thetaWS*rev(vectR),errSW,errWS,Nchrom,CatInf))
	  lnL <- sum(lnLN + lnLWS + lnLSW)
	} else {
    lnL <- 0
		for (i in 1:Ncat){
			j1 <- Cinf[i]
			j2 <- Csup[i]
			Kn <- 0
			Ln <- 0
			Kws <- 0
			Lws <- 0
			Ksw <- 0
			Lsw <- 0
			for (k in j1:j2){
				Kn <- Kn + SFSneutral[k]
				Ln <- Ln + ENS_neutral_err(theta*vectR[i],theta*vectR[Ncat+1-i],errN,errN,Nchrom,k)
				Kws <- Kws + SFSws[k]
				Lws <- Lws + ENS_neutral_err((1-GC)*thetaWS*vectR[i],GC*bias*thetaWS*vectR[Ncat+1-i],errWS,errSW,Nchrom,k)
				Ksw <- Ksw + SFSsw[k]	
				Lsw	<- Lsw + ENS_neutral_err(GC*bias*thetaWS*vectR[i],(1-GC)*thetaWS*vectR[Ncat+1-i],errSW,errWS,Nchrom,k)	
			}
			lnL <- lnL + log_Poisson(Kn,Ln)
			lnL <- lnL + log_Poisson(Kws,Lws)
			lnL <- lnL + log_Poisson(Ksw,Lsw)
		}		
	}
	return(-lnL)
}

# Model 1: constant gBGC
# theta,thetaWS,bias,B,vectR errN,eerWs,errSW: parameters to be optimized
# SFSneutral,SFSws,SFSsw,GC: data
lnL1err3 <- function(par,SFSneutral,SFSws,SFSsw,GC,CatInf=c(1:length(SFSneutral)),CatSup=CatInf){
  # Constraining categories to be symetric
  Nchrom <- length(SFSneutral)+1
  Cinf <- union(CatInf[CatInf<=floor(Nchrom/2)],Nchrom-rev(CatSup[CatSup<=floor(Nchrom/2)]))
  Csup <- union(CatSup[CatSup<=floor(Nchrom/2)],Nchrom-rev(CatInf[CatInf<=floor(Nchrom/2)]))
  Ncat <- length(Cinf)
  theta <- par[1]
  thetaWS <- par[2]
  bias <- par[3]
  B <- par[4]
  vectR <- c(1,par[5:(3+Ncat)])
  errN <- par[Ncat+4]
  errWS <- par[Ncat+5]
  errSW <- par[Ncat+6] 
  if(identical(Cinf,Csup)){
    sfsN <- SFSneutral[CatInf]
    sfsWS <- SFSws[CatInf]
    sfsSW <- SFSsw[CatInf]
    lnLN <- log_Poisson(sfsN,ENS_neutral_err(theta*vectR,theta*rev(vectR),errN,errN,Nchrom,CatInf))
    lnLWS <- log_Poisson(sfsWS,ENS_constant_err((1-GC)*thetaWS*vectR,GC*bias*thetaWS*rev(vectR),B,errWS,errSW,Nchrom,CatInf))
    lnLSW <- log_Poisson(sfsSW,ENS_constant_err(GC*bias*thetaWS*vectR,(1-GC)*thetaWS*rev(vectR),-B,errSW,errWS,Nchrom,CatInf))
    lnL <- sum(lnLN + lnLWS + lnLSW)
  } else {
    lnL <- 0
    for (i in 1:Ncat){
      j1 <- Cinf[i]
      j2 <- Csup[i]
      Kn <- 0
      Ln <- 0
      Kws <- 0
      Lws <- 0
      Ksw <- 0
      Lsw <- 0
      for (k in j1:j2){
        Kn <- Kn + SFSneutral[k]
        Ln <- Ln + ENS_neutral_err(theta*vectR[i],theta*vectR[Ncat+1-i],errN,errN,Nchrom,k)
        Kws <- Kws + SFSws[k]
        Lws <- Lws + ENS_constant_err((1-GC)*thetaWS*vectR[i],GC*bias*thetaWS*vectR[Ncat+1-i],B,errWS,errSW,Nchrom,k)
        Ksw <- Ksw + SFSsw[k]	
        Lsw	<- Lsw + ENS_constant_err(GC*bias*thetaWS*vectR[i],(1-GC)*thetaWS*vectR[Ncat+1-i],-B,errSW,errWS,Nchrom,k)	
      }
      lnL <- lnL + log_Poisson(Kn,Ln)
      lnL <- lnL + log_Poisson(Kws,Lws)
      lnL <- lnL + log_Poisson(Ksw,Lsw)
    }		
  }
  return(-lnL)
}


# Model 2a: hotspot gBGC
# theta,thetaWS,bias,B,f,vectR errN,eerWs,errSW: parameters to be optimized
# SFSneutral,SFSws,SFSsw,GC: data
lnL2aerr3 <- function(par,SFSneutral,SFSws,SFSsw,GC,CatInf=c(1:length(SFSneutral)),CatSup=CatInf){
  # Constraining categories to be symetric
  Nchrom <- length(SFSneutral)+1
  Cinf <- union(CatInf[CatInf<=floor(Nchrom/2)],Nchrom-rev(CatSup[CatSup<=floor(Nchrom/2)]))
  Csup <- union(CatSup[CatSup<=floor(Nchrom/2)],Nchrom-rev(CatInf[CatInf<=floor(Nchrom/2)]))
  Ncat <- length(Cinf)
  theta <- par[1]
  thetaWS <- par[2]
  bias <- par[3]
  B <- par[4]
  f <- par[5]
  vectR <- c(1,par[6:(4+Ncat)])
  errN <- par[Ncat+5]
  errWS <- par[Ncat+6]
  errSW <- par[Ncat+7] 
  if(identical(Cinf,Csup)){
    sfsN <- SFSneutral[CatInf]
    sfsWS <- SFSws[CatInf]
    sfsSW <- SFSsw[CatInf]
    lnLN <- log_Poisson(sfsN,ENS_neutral_err(theta*vectR,theta*rev(vectR),errN,errN,Nchrom,CatInf))
    lnLWS <- log_Poisson(sfsWS,ENS_hotspot1_err((1-GC)*thetaWS*vectR,GC*bias*thetaWS*rev(vectR),B,f,errWS,errSW,Nchrom,CatInf))
    lnLSW <- log_Poisson(sfsSW,ENS_hotspot1_err(GC*bias*thetaWS*vectR,(1-GC)*thetaWS*rev(vectR),-B,f,errSW,errWS,Nchrom,CatInf))
    lnL <- sum(lnLN + lnLWS + lnLSW)
  } else {
    lnL <- 0
    for (i in 1:Ncat){
      j1 <- Cinf[i]
      j2 <- Csup[i]
      Kn <- 0
      Ln <- 0
      Kws <- 0
      Lws <- 0
      Ksw <- 0
      Lsw <- 0
      for (k in j1:j2){
        Kn <- Kn + SFSneutral[k]
        Ln <- Ln + ENS_neutral_err(theta*vectR[i],theta*vectR[Ncat+1-i],errN,errN,Nchrom,k)
        Kws <- Kws + SFSws[k]
        Lws <- Lws + ENS_hotspot1_err((1-GC)*thetaWS*vectR[i],GC*bias*thetaWS*vectR[Ncat+1-i],B,f,errWS,errSW,Nchrom,k)
        Ksw <- Ksw + SFSsw[k]  
        Lsw	<- Lsw + ENS_hotspot1_err(GC*bias*thetaWS*vectR[i],(1-GC)*thetaWS*vectR[Ncat+1-i],-B,f,errSW,errWS,Nchrom,k)	
      }
      lnL <- lnL + log_Poisson(Kn,Ln)
      lnL <- lnL + log_Poisson(Kws,Lws)
      lnL <- lnL + log_Poisson(Ksw,Lsw)
    }		
  }
  return(-lnL)
}


# Model 2b: two gBGC categories
# theta,thetaWS,bias,B0,B1,f,vectR errN,eerWs,errSW: parameters to be optimized
# SFSneutral,SFSws,SFSsw,GC: data
lnL2berr3 <- function(par,SFSneutral,SFSws,SFSsw,GC,CatInf=c(1:length(SFSneutral)),CatSup=CatInf){
  # Constraining categories to be symetric
  Nchrom <- length(SFSneutral)+1
  Cinf <- union(CatInf[CatInf<=floor(Nchrom/2)],Nchrom-rev(CatSup[CatSup<=floor(Nchrom/2)]))
  Csup <- union(CatSup[CatSup<=floor(Nchrom/2)],Nchrom-rev(CatInf[CatInf<=floor(Nchrom/2)]))
  Ncat <- length(Cinf)
  theta <- par[1]
  thetaWS <- par[2]
  bias <- par[3]
  B0 <- par[4]
  B1 <- par[5]
  f <- par[6]
  vectR <- c(1,par[7:(5+Ncat)])
  errN <- par[Ncat+6]
  errWS <- par[Ncat+7]
  errSW <- par[Ncat+8] 
  if(identical(Cinf,Csup)){
    sfsN <- SFSneutral[CatInf]
    sfsWS <- SFSws[CatInf]
    sfsSW <- SFSsw[CatInf]
    lnLN <- log_Poisson(sfsN,ENS_neutral_err(theta*vectR,theta*rev(vectR),errN,errN,Nchrom,CatInf))
    lnLWS <- log_Poisson(sfsWS,ENS_hotspot2_err((1-GC)*thetaWS*vectR,GC*bias*thetaWS*rev(vectR),B0,B1,f,errWS,errSW,Nchrom,CatInf))
    lnLSW <- log_Poisson(sfsSW,ENS_hotspot2_err(GC*bias*thetaWS*vectR,(1-GC)*thetaWS*rev(vectR),-B0,-B1,f,errSW,errWS,Nchrom,CatInf))
    lnL <- sum(lnLN + lnLWS + lnLSW)
  } else {
    lnL <- 0
    for (i in 1:Ncat){
      j1 <- Cinf[i]
      j2 <- Csup[i]
      Kn <- 0
      Ln <- 0
      Kws <- 0
      Lws <- 0
      Ksw <- 0
      Lsw <- 0
      for (k in j1:j2){
        Kn <- Kn + SFSneutral[k]
        Ln <- Ln + ENS_neutral_err(theta*vectR[i],theta*vectR[Ncat+1-i],errN,errN,Nchrom,k)
        Kws <- Kws + SFSws[k]
        Lws <- Lws + ENS_hotspot2_err((1-GC)*thetaWS*vectR[i],GC*bias*thetaWS*vectR[Ncat+1-i],B0,B1,f,errWS,errSW,Nchrom,k)
        Ksw <- Ksw + SFSsw[k]  
        Lsw  <- Lsw + ENS_hotspot2_err(GC*bias*thetaWS*vectR[i],(1-GC)*thetaWS*vectR[Ncat+1-i],-B0,-B1,f,errSW,errWS,Nchrom,k)	
      }
      lnL <- lnL + log_Poisson(Kn,Ln)
      lnL <- lnL + log_Poisson(Kws,Lws)
      lnL <- lnL + log_Poisson(Ksw,Lsw)
    }		
  }
  return(-lnL)
}
