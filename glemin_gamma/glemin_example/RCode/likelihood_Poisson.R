
########################
# Likelihood functions #
########################

# Poisson distribution
######################


# source("/Users/sylvain/Bin/R/estimBGC/basic_functions.R")

# Null model: no gBGC
# theta,thetaWS,bias,vectR: parameters to be optimized
# SFSneutral,SFSws,SFSsw,GC: data
# CatInf and CatSup: optional vectors to chose the SNPs and the grouping of SNPs to be used
lnL0 <- function(par,SFSneutral,SFSws,SFSsw,GC,CatInf=c(1:length(SFSneutral)),CatSup=CatInf){
	Ncat <- length(CatInf)
	theta<-par[1]
	thetaWS<-par[2]
	bias<-par[3]
	vectR<-c(1,par[4:(2+Ncat)])
	if(identical(CatInf,CatSup)){
		sfsN <- SFSneutral[CatInf]
		sfsWS <- SFSws[CatInf]
		sfsSW <- SFSsw[CatInf]
		lnLN <- log_Poisson(sfsN,ENS_neutral(theta*vectR,CatInf))
		lnLWS <- log_Poisson(sfsWS,ENS_neutral((1-GC)*thetaWS*vectR,CatInf))
		lnLSW <- log_Poisson(sfsSW,ENS_neutral(GC*bias*thetaWS*vectR,CatInf))
		lnL <- sum(lnLN + lnLWS + lnLSW)
	} else {
    lnL <- 0
		for (i in 1:Ncat){
			j1 <- CatInf[i]
			j2 <- CatSup[i]
			Kn <- 0
			Ln <- 0
			Kws <- 0
			Lws <- 0
			Ksw <- 0
			Lsw <- 0
			for (k in j1:j2){
				Kn <- Kn + SFSneutral[k]
				Ln <- Ln + ENS_neutral(theta*vectR[i],k)
				Kws <- Kws + SFSws[k]
				Lws <- Lws + ENS_neutral((1-GC)*thetaWS*vectR[i],k)
				Ksw <- Ksw + SFSsw[k]	
				Lsw	<- Lsw + ENS_neutral(GC*bias*thetaWS*vectR[i],k)	
			}
			lnL <- lnL + log_Poisson(Kn,Ln)
			lnL <- lnL + log_Poisson(Kws,Lws)
			lnL <- lnL + log_Poisson(Ksw,Lsw)
		}		
	}
	return(-lnL)
}

# Model 1: constant gBGC
# theta,thetaWS,bias,B,vectR: parameters to be optimized
# SFSneutral,SFSws,SFSsw,GC: data
lnL1 <- function(par,SFSneutral,SFSws,SFSsw,GC,CatInf=c(1:length(SFSneutral)),CatSup=CatInf){
	Ncat <- length(CatInf)
	Nchrom <- length(SFSneutral) +1
	theta <- par[1]
	thetaWS <- par[2]
	bias <- par[3]
	B <- par[4]
	vectR <- c(1,par[5:(3+Ncat)])
	if(identical(CatInf,CatSup)){
	  sfsN <- SFSneutral[CatInf]
	  sfsWS <- SFSws[CatInf]
	  sfsSW <- SFSsw[CatInf]
	  lnLN <- log_Poisson(sfsN,ENS_neutral(theta*vectR,CatInf))
	  lnLWS <- log_Poisson(sfsWS,ENS_constant((1-GC)*thetaWS*vectR,B,Nchrom,CatInf))
	  lnLSW <- log_Poisson(sfsSW,ENS_constant(GC*bias*thetaWS*vectR,-B,Nchrom,CatInf))
	  lnL <- sum(lnLN + lnLWS + lnLSW)
	} else {
    lnL <- 0
		for (i in 1:Ncat){
			j1 <- CatInf[i]
			j2 <- CatSup[i]
			Kn <- 0
			Ln <- 0
			Kws <- 0
			Lws <- 0
			Ksw <- 0
			Lsw <- 0
			for (k in j1:j2){
				Kn <- Kn + SFSneutral[k]
				Ln <- Ln + ENS_neutral(theta*vectR[i],k)
				Kws <- Kws + SFSws[k]
				Lws <- Lws + ENS_constant((1-GC)*thetaWS*vectR[i],B,Nchrom,k)
				Ksw <- Ksw + SFSsw[k]	
				Lsw	<- Lsw + ENS_constant(GC*bias*thetaWS*vectR[i],-B,Nchrom,k)	
			}
			lnL <- lnL + log_Poisson(Kn,Ln)
			lnL <- lnL + log_Poisson(Kws,Lws)
			lnL <- lnL + log_Poisson(Ksw,Lsw)
		}		
	} 
	return(-lnL)
}

# Model 2a: hotspot 1
# theta,thetaWS,bias,B,f,vectR: parameters to be optimized
# SFSneutral,SFSws,SFSsw,GC: data
lnL2a <- function(par,SFSneutral,SFSws,SFSsw,GC,CatInf=c(1:length(SFSneutral)),CatSup=CatInf){
	Ncat <- length(CatInf)
	Nchrom <- length(SFSneutral) +1
	theta<-par[1]
	thetaWS<-par[2]
	bias<-par[3]
	B<-par[4]
	f<-par[5]
	vectR<-c(1,par[6:(4+Ncat)])
	if(identical(CatInf,CatSup)){
	  sfsN <- SFSneutral[CatInf]
	  sfsWS <- SFSws[CatInf]
	  sfsSW <- SFSsw[CatInf]
	  lnLN <- log_Poisson(sfsN,ENS_neutral(theta*vectR,CatInf))
	  lnLWS <- log_Poisson(sfsWS,ENS_hotspot1((1-GC)*thetaWS*vectR,B,f,Nchrom,CatInf))
	  lnLSW <- log_Poisson(sfsSW,ENS_hotspot1(GC*bias*thetaWS*vectR,-B,f,Nchrom,CatInf))
	  lnL <- sum(lnLN + lnLWS + lnLSW)
	} else {
    lnL <- 0
		for (i in 1:Ncat){
			j1 <- CatInf[i]
			j2 <- CatSup[i]
			Kn <- 0
			Ln <- 0
			Kws <- 0
			Lws <- 0
			Ksw <- 0
			Lsw <- 0
			for (k in j1:j2){
				Kn <- Kn + SFSneutral[k]
				Ln <- Ln + ENS_neutral(theta*vectR[i],k)
				Kws <- Kws + SFSws[k]
				Lws <- Lws + ENS_hotspot1((1-GC)*thetaWS*vectR[i],B,f,Nchrom,k)
				Ksw <- Ksw + SFSsw[k]	
				Lsw	<- Lsw + ENS_hotspot1(GC*bias*thetaWS*vectR[i],-B,f,Nchrom,k)	
			}
			lnL <- lnL + log_Poisson(Kn,Ln)
			lnL <- lnL + log_Poisson(Kws,Lws)
			lnL <- lnL + log_Poisson(Ksw,Lsw)
		}		
	} 
	return(-lnL)
}

# Model 2b: hotspot 2
# theta,thetaWS,bias,B0,B1,f,vectR: parameters to be optimized
# SFSneutral,SFSws,SFSsw,GC: data
lnL2b <- function(par,SFSneutral,SFSws,SFSsw,GC,CatInf=c(1:length(SFSneutral)),CatSup=CatInf){
	Ncat <- length(CatInf)
	Nchrom <- length(SFSneutral) +1
	theta<-par[1]
	thetaWS<-par[2]
	bias<-par[3]
	B0<-par[4]
	B1<-par[5]
	f<-par[6]
	vectR<-c(1,par[7:(5+Ncat)])
	if(identical(CatInf,CatSup)){
	  sfsN <- SFSneutral[CatInf]
	  sfsWS <- SFSws[CatInf]
	  sfsSW <- SFSsw[CatInf]
	  lnLN <- log_Poisson(sfsN,ENS_neutral(theta*vectR,CatInf))
	  lnLWS <- log_Poisson(sfsWS,ENS_hotspot2((1-GC)*thetaWS*vectR,B0,B1,f,Nchrom,CatInf))
	  lnLSW <- log_Poisson(sfsSW,ENS_hotspot2(GC*bias*thetaWS*vectR,-B0,-B1,f,Nchrom,CatInf))
	  lnL <- sum(lnLN + lnLWS + lnLSW)
	} else {
    lnL <- 0
		for (i in 1:Ncat){
			j1 <- CatInf[i]
			j2 <- CatSup[i]
			Kn <- 0
			Ln <- 0
			Kws <- 0
			Lws <- 0
			Ksw <- 0
			Lsw <- 0
			for (k in j1:j2){
				Kn <- Kn + SFSneutral[k]
				Ln <- Ln + ENS_neutral(theta*vectR[i],k)
				Kws <- Kws + SFSws[k]
				Lws <- Lws + ENS_hotspot2((1-GC)*thetaWS*vectR[i],B0,B1,f,Nchrom,k)
				Ksw <- Ksw + SFSsw[k]	
				Lsw	<- Lsw + ENS_hotspot2(GC*bias*thetaWS*vectR[i],-B0,-B1,f,Nchrom,k)	
			}
			lnL <- lnL + log_Poisson(Kn,Ln)
			lnL <- lnL + log_Poisson(Kws,Lws)
			lnL <- lnL + log_Poisson(Ksw,Lsw)
		}		
	} 
	return(-lnL)
}