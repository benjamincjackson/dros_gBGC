# How to use the esimBGC program ####

# Two specific library are needed: "HMP" and "VGAM" 

# Loading the code

source("/Users/sylvain/Documents/Boulot/Recherche/Thematiques/GC/EstimationBGC/BGCSpectre/Article/V11/SupFiles/RCode/optimizeBGC.R")


# Dataset ####
# Simulation with B=1 thetaN=1000, thetaWS=2000 mut.bias=2
GC <- 0.5
Sn <- c(102,52,34,24,22,18,14,10)
Sws <- c(214,111,77,57,49,42,42,37)
Ssw <- c(386,172,113,80,55,51,39,34)




# Lilekilhood of the saturated model
lnLsat <- saturatedlnL(Sn) + saturatedlnL(Sws) + saturatedlnL(Ssw)


# Models without orientation error : Default usage
result1 <- optimizeBGC(Sn,Sws,Ssw,GC,PrintInit=T)



# Models with orientation errors
result2 <- optimizeBGCwith3error(Sn,Sws,Ssw,GC,verbose=1,MAXIT=500)

# Comparison of the likelihoods
result1$model0$lnL
result1$model1$lnL
result2$model0$lnL
result2$model1$lnL


# Hotspot models

# Without error
result3 <- optimizeBGC(Sn,Sws,Ssw,GC,models=c(T,T),MAXIT=200,verbose=1)
result3$model0$lnL
result3$model1$lnL
result3$model2a$lnL
result3$model2b$lnL



# With errors
result4 <- optimizeBGCwith3error(Sn,Sws,Ssw,GC,models=c(T,T))
result4$model0$lnL
result4$model1$lnL
result4$model2a$lnL
result4$model2b$lnL







