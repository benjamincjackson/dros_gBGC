setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# carry out the glemin gamma model fitting and comparison for simulans, 
# for each bin in turn

# need to load the SFSes and bin indices and the alignments:
load('../data/SI_binindices_GCcontent.RData')
load('../data/est-sfs_SFSs.RData')
load('../data/SI_alignments.RData')

# and get a GC content for the polymorphism samples after concatenating all sites
# in a bin

source('../scripts/R_helper_functions.R')
library(ape)

GCcontent.Dsim.MD.SI.A.concat.mean.5GCbins <- sapply(
  get_concat_seqs(alignments.intersection.MD21.SI.A,
                  GCbinindex.5bins.mean.SI.A),
  FUN = function(x){
    GC.content(as.DNAbin(editAlignment(x, 'MD')))
  })

GCcontent.Dmel.ZI69.SI.A.concat.mean.5GCbins <- sapply(
  get_concat_seqs(alignments.intersection.ZI69.SI.A,
                  GCbinindex.5bins.mean.SI.A),
  FUN = function(x){
    GC.content(as.DNAbin(editAlignment(x, 'ZI')))
  })

GCcontent.Dmel.ZI21.SI.A.concat.mean.5GCbins <- sapply(
  get_concat_seqs(alignments.intersection.ZI21.SI.A,
                  GCbinindex.5bins.mean.SI.A),
  FUN = function(x){
    GC.content(as.DNAbin(editAlignment(x, 'ZI')))
  })

GCcontent.Dmel.RG21.SI.A.concat.mean.5GCbins <- sapply(
  get_concat_seqs(alignments.intersection.RG21.SI.A,
                  GCbinindex.5bins.mean.SI.A),
  FUN = function(x){
    GC.content(as.DNAbin(editAlignment(x, 'RG')))
  })

GCcontent.Dsim.MD.SI.A.concat.sim.5GCbins <- sapply(
  get_concat_seqs(alignments.intersection.MD21.SI.A,
                  GCbinindex.5bins.Dsim.SI.A),
  FUN = function(x){
    GC.content(as.DNAbin(editAlignment(x, 'MD')))
  })

GCcontent.Dmel.ZI69.SI.A.concat.mel.5GCbins <- sapply(
  get_concat_seqs(alignments.intersection.ZI69.SI.A,
                  GCbinindex.5bins.Dmel.SI.A),
  FUN = function(x){
    GC.content(as.DNAbin(editAlignment(x, 'ZI')))
  })

GCcontent.Dmel.ZI21.SI.A.concat.mel.5GCbins <- sapply(
  get_concat_seqs(alignments.intersection.ZI21.SI.A,
                  GCbinindex.5bins.Dmel.SI.A),
  FUN = function(x){
    GC.content(as.DNAbin(editAlignment(x, 'ZI')))
  })

GCcontent.Dmel.RG21.SI.A.concat.mel.5GCbins <- sapply(
  get_concat_seqs(alignments.intersection.RG21.SI.A,
                  GCbinindex.5bins.Dmel.SI.A),
  FUN = function(x){
    GC.content(as.DNAbin(editAlignment(x, 'RG')))
  })
################################################################################

# The Gl?min method requires the R libraries "HMP" and "VGAM" 
# install.packages('HMP')
# install.packages('VGAM')

library(HMP)
library(VGAM)

# need to source the 'optimizeBGC' functions, and others
source('../glemin_gamma/paper/Glemin_supplement/RCode/optimizeBGC_DIFF_ALGORITHMS.R')
source('../glemin_gamma/paper/Glemin_supplement/RCode/basic_functions.R')
source('../glemin_gamma/paper/Glemin_supplement/RCode/likelihood_Poisson.R')
source('../glemin_gamma/paper/Glemin_supplement/RCode/likelihood_Error_Poisson.R')

# A function to take the 3 input SFSes and the 1 GC content and perform
# the tests, and save all the data, probably just in a list
# IN THIS VERSION, WE CAN PASS THE METHOD USED BY THE optim() FUNCTION: BY CHANGING THIS,
get_test_results_no_error <- function(SFS_Neu, SFS_S_2_W, SFS_W_2_S, GC,
                                      algorithmToUse, maxit = 1000){
  
  result_no_error <- optimizeBGCDiffAl(dataNeutral = SFS_Neu,
                                       dataWS = SFS_W_2_S,
                                       dataSW = SFS_S_2_W,
                                       dataGC = GC,
                                       PrintInit=T,
                                       optimAlgorithm = algorithmToUse,
                                       MAXIT = maxit,
                                       verbose = 2)
  
  # MAKE A CONVERGENCE TABLE, TO CHECK WHAT HAS/HASN'T WORKED MORE EASILY?
  conv_table <- matrix(data = c(result_no_error$model0$convergence,
                                result_no_error$model1$convergence),
                       nrow = 2, ncol = 1,
                       dimnames = list(c('M_0', 'M_1'),
                                       c('no_error')))
  
  lnlk_table <- matrix(data = c(result_no_error$model0$lnL,
                                result_no_error$model1$lnL),
                       nrow = 2, ncol = 1,
                       dimnames = list(c('M_0', 'M_1'),
                                       c('no_error')))
  
  #return the full results, but also the LnLks in a more convenient form?
  return(list('without_error' = result_no_error,
              'lnlk_table' = lnlk_table,
              'conv_table' = conv_table))
}

# 
results_MD.SI.A.5GCbins.mean_no_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    results_MD.SI.A.5GCbins.mean_no_error[[i]] <- get_test_results_no_error(SFS_Neu = N_to_N_ALL_SFSes.MD.SI.A.mean.5GCbins[[i]][2:21],
                                                               SFS_S_2_W = GC_to_AT_ALL_SFSes.MD.SI.A.mean.5GCbins[[i]][2:21],
                                                               SFS_W_2_S = AT_to_GC_ALL_SFSes.MD.SI.A.mean.5GCbins[[i]][2:21],
                                                               GC = GCcontent.Dsim.MD.SI.A.concat.mean.5GCbins[i],
                                                               algorithmToUse = "Nelder-Mead",
                                                               maxit=1000)
  }
)

results_MD.SI.A.5GCbins.sim_no_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    results_MD.SI.A.5GCbins.sim_no_error[[i]] <- get_test_results_no_error(SFS_Neu = N_to_N_ALL_SFSes.MD.SI.A.sim.5GCbins[[i]][2:21],
                                                                            SFS_S_2_W = GC_to_AT_ALL_SFSes.MD.SI.A.sim.5GCbins[[i]][2:21],
                                                                            SFS_W_2_S = AT_to_GC_ALL_SFSes.MD.SI.A.sim.5GCbins[[i]][2:21],
                                                                            GC = GCcontent.Dsim.MD.SI.A.concat.sim.5GCbins[i],
                                                                            algorithmToUse = "Nelder-Mead",
                                                                            maxit=1000)
  }
)

results_ZI69.SI.A.5GCbins.mean_no_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    results_ZI69.SI.A.5GCbins.mean_no_error[[i]] <- get_test_results_no_error(SFS_Neu = N_to_N_ALL_SFSes.ZI69.SI.A.mean.5GCbins[[i]][2:69],
                                                               SFS_S_2_W = GC_to_AT_ALL_SFSes.ZI69.SI.A.mean.5GCbins[[i]][2:69],
                                                               SFS_W_2_S = AT_to_GC_ALL_SFSes.ZI69.SI.A.mean.5GCbins[[i]][2:69],
                                                               GC = GCcontent.Dmel.ZI69.SI.A.concat.mean.5GCbins[i],
                                                               algorithmToUse = "Nelder-Mead",
                                                               maxit=1000)
  }
)

results_ZI69.SI.A.5GCbins.mel_no_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    results_ZI69.SI.A.5GCbins.mel_no_error[[i]] <- get_test_results_no_error(SFS_Neu = N_to_N_ALL_SFSes.ZI69.SI.A.mel.5GCbins[[i]][2:69],
                                                                              SFS_S_2_W = GC_to_AT_ALL_SFSes.ZI69.SI.A.mel.5GCbins[[i]][2:69],
                                                                              SFS_W_2_S = AT_to_GC_ALL_SFSes.ZI69.SI.A.mel.5GCbins[[i]][2:69],
                                                                              GC = GCcontent.Dmel.ZI69.SI.A.concat.mel.5GCbins[i],
                                                                              algorithmToUse = "Nelder-Mead",
                                                                              maxit=1000)
  }
)


results_ZI21.SI.A.5GCbins.mean_no_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    results_ZI21.SI.A.5GCbins.mean_no_error[[i]] <- get_test_results_no_error(SFS_Neu = N_to_N_ALL_SFSes.ZI21.SI.A.mean.5GCbins[[i]][2:21],
                                                                         SFS_S_2_W = GC_to_AT_ALL_SFSes.ZI21.SI.A.mean.5GCbins[[i]][2:21],
                                                                         SFS_W_2_S = AT_to_GC_ALL_SFSes.ZI21.SI.A.mean.5GCbins[[i]][2:21],
                                                                         GC = GCcontent.Dmel.ZI21.SI.A.concat.mean.5GCbins[i],
                                                                         algorithmToUse = "Nelder-Mead",
                                                                         maxit=1000)
  }
)

results_ZI21.SI.A.5GCbins.mel_no_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    results_ZI21.SI.A.5GCbins.mel_no_error[[i]] <- get_test_results_no_error(SFS_Neu = N_to_N_ALL_SFSes.ZI21.SI.A.mel.5GCbins[[i]][2:21],
                                                                              SFS_S_2_W = GC_to_AT_ALL_SFSes.ZI21.SI.A.mel.5GCbins[[i]][2:21],
                                                                              SFS_W_2_S = AT_to_GC_ALL_SFSes.ZI21.SI.A.mel.5GCbins[[i]][2:21],
                                                                              GC = GCcontent.Dmel.ZI21.SI.A.concat.mel.5GCbins[i],
                                                                              algorithmToUse = "Nelder-Mead",
                                                                              maxit=1000)
  }
)

results_RG21.SI.A.5GCbins.mean_no_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    results_RG21.SI.A.5GCbins.mean_no_error[[i]] <- get_test_results_no_error(SFS_Neu = N_to_N_ALL_SFSes.RG21.SI.A.mean.5GCbins[[i]][2:21],
                                                                         SFS_S_2_W = GC_to_AT_ALL_SFSes.RG21.SI.A.mean.5GCbins[[i]][2:21],
                                                                         SFS_W_2_S = AT_to_GC_ALL_SFSes.RG21.SI.A.mean.5GCbins[[i]][2:21],
                                                                         GC = GCcontent.Dmel.RG21.SI.A.concat.mean.5GCbins[i],
                                                                         algorithmToUse = "Nelder-Mead",
                                                                         maxit=1000)
  }
)

results_RG21.SI.A.5GCbins.mel_no_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    results_RG21.SI.A.5GCbins.mel_no_error[[i]] <- get_test_results_no_error(SFS_Neu = N_to_N_ALL_SFSes.RG21.SI.A.mel.5GCbins[[i]][2:21],
                                                                              SFS_S_2_W = GC_to_AT_ALL_SFSes.RG21.SI.A.mel.5GCbins[[i]][2:21],
                                                                              SFS_W_2_S = AT_to_GC_ALL_SFSes.RG21.SI.A.mel.5GCbins[[i]][2:21],
                                                                              GC = GCcontent.Dmel.RG21.SI.A.concat.mel.5GCbins[i],
                                                                              algorithmToUse = "Nelder-Mead",
                                                                              maxit=1000)
  }
)

# check for convergence:
lapply(results_MD.SI.A.5GCbins.mean_no_error, FUN = function(x){x$conv_table})
lapply(results_ZI69.SI.A.5GCbins.mean_no_error, FUN = function(x){x$conv_table})
lapply(results_ZI21.SI.A.5GCbins.mean_no_error, FUN = function(x){x$conv_table})
lapply(results_RG21.SI.A.5GCbins.mean_no_error, FUN = function(x){x$conv_table})
lapply(results_MD.SI.A.5GCbins.sim_no_error, FUN = function(x){x$conv_table})
lapply(results_ZI69.SI.A.5GCbins.mel_no_error, FUN = function(x){x$conv_table})
lapply(results_ZI21.SI.A.5GCbins.mel_no_error, FUN = function(x){x$conv_table})
lapply(results_RG21.SI.A.5GCbins.mel_no_error, FUN = function(x){x$conv_table})


# then run the models with polarisation error, using the ML estimates of parameters from
# the models without error as intiation values where possible
get_test_results_error <- function(SFS_Neu, SFS_S_2_W, SFS_W_2_S, GC,
                                   algorithmToUse, maxit = 1000,
                                   MYtheta_init0=NULL,MYthetaWS_init0=NULL,MYbias_init0=NULL,MYerr_init0=NULL,MYerrWS_init0=NULL,MYerrSW_init0=NULL, #initial parameters for model0
                                   MYtheta_init1=NULL,MYthetaWS_init1=NULL,MYbias_init1=NULL,MYB_init1=NULL,MYerr_init1=NULL,MYerrWS_init1=NULL,MYerrSW_init1=NULL #initial parameters for model1
){
  
  result_error <- optimizeBGCwith3errorDiffAl(dataNeutral = SFS_Neu,
                                              dataWS = SFS_W_2_S,
                                              dataSW = SFS_S_2_W,
                                              dataGC = GC,
                                              optimAlgorithm = algorithmToUse,
                                              MAXIT = maxit,
                                              verbose = 2,
                                              theta_init0=MYtheta_init0,thetaWS_init0=MYthetaWS_init0,bias_init0=MYbias_init0,err_init0=MYerr_init0,errWS_init0=MYerrWS_init0,errSW_init0=MYerrSW_init0, #initial parameters for model0
                                              theta_init1=MYtheta_init1,thetaWS_init1=MYthetaWS_init1,bias_init1=MYbias_init1,B_init1=MYB_init1,err_init1=MYerr_init1,errWS_init1=MYerrWS_init1,errSW_init1=MYerrSW_init1 #initial parameters for model1
  )
  
  # MAKE A CONVERGENCE TABLE, TO CHECK WHAT HAS/HASN'T WORKED MORE EASILY?
  conv_table <- matrix(data = c(result_error$model0$convergence,
                                result_error$model1$convergence),
                       nrow = 2, ncol = 1,
                       dimnames = list(c('M_0', 'M_1'),
                                       c('error')))
  
  lnlk_table <- matrix(data = c(result_error$model0$lnL,
                                result_error$model1$lnL),
                       nrow = 2, ncol = 1,
                       dimnames = list(c('M_0', 'M_1'),
                                       c('error')))
  
  #return the full results, but also the LnLks in a more convenient form?
  return(list('with_error' = result_error,
              'lnlk_table' = lnlk_table,
              'conv_table' = conv_table))
}

results_MD.SI.A.5GCbins.mean_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    # this is just to make the arguments list in the next part less long:
    x <- results_MD.SI.A.5GCbins.mean_no_error[[i]]$without_error$model0
    y <- results_MD.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1
    
    results_MD.SI.A.5GCbins.mean_error[[i]] <- get_test_results_error(SFS_Neu = N_to_N_ALL_SFSes.MD.SI.A.mean.5GCbins[[i]][2:21],
                                                                  SFS_S_2_W = GC_to_AT_ALL_SFSes.MD.SI.A.mean.5GCbins[[i]][2:21],
                                                                  SFS_W_2_S = AT_to_GC_ALL_SFSes.MD.SI.A.mean.5GCbins[[i]][2:21],
                                                                  GC = GCcontent.Dsim.MD.SI.A.concat.mean.5GCbins[i],
                                                         algorithmToUse = "Nelder-Mead",
                                                         maxit=10000,
                                                         MYtheta_init0=x$theta,MYthetaWS_init0=x$thetaWS,MYbias_init0=x$bias,MYerr_init0=NULL,MYerrWS_init0=NULL,MYerrSW_init0=NULL, #initial parameters for model0
                                                         MYtheta_init1=y$theta,MYthetaWS_init1=y$thetaWS,MYbias_init1=y$bias,MYB_init1=y$B,MYerr_init1=NULL,MYerrWS_init1=NULL,MYerrSW_init1=NULL)
  })

results_MD.SI.A.5GCbins.sim_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    # this is just to make the arguments list in the next part less long:
    x <- results_MD.SI.A.5GCbins.sim_no_error[[i]]$without_error$model0
    y <- results_MD.SI.A.5GCbins.sim_no_error[[i]]$without_error$model1
    
    results_MD.SI.A.5GCbins.sim_error[[i]] <- get_test_results_error(SFS_Neu = N_to_N_ALL_SFSes.MD.SI.A.sim.5GCbins[[i]][2:21],
                                                                      SFS_S_2_W = GC_to_AT_ALL_SFSes.MD.SI.A.sim.5GCbins[[i]][2:21],
                                                                      SFS_W_2_S = AT_to_GC_ALL_SFSes.MD.SI.A.sim.5GCbins[[i]][2:21],
                                                                      GC = GCcontent.Dsim.MD.SI.A.concat.sim.5GCbins[i],
                                                                      algorithmToUse = "Nelder-Mead",
                                                                      maxit=10000,
                                                                      MYtheta_init0=x$theta,MYthetaWS_init0=x$thetaWS,MYbias_init0=x$bias,MYerr_init0=NULL,MYerrWS_init0=NULL,MYerrSW_init0=NULL, #initial parameters for model0
                                                                      MYtheta_init1=y$theta,MYthetaWS_init1=y$thetaWS,MYbias_init1=y$bias,MYB_init1=y$B,MYerr_init1=NULL,MYerrWS_init1=NULL,MYerrSW_init1=NULL)
  })

results_ZI69.SI.A.5GCbins.mean_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    # this is just to make the arguments list in the next part less long:
    x <- results_ZI69.SI.A.5GCbins.mean_no_error[[i]]$without_error$model0
    y <- results_ZI69.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1
    
    results_ZI69.SI.A.5GCbins.mean_error[[i]] <- get_test_results_error(SFS_Neu = N_to_N_ALL_SFSes.ZI69.SI.A.mean.5GCbins[[i]][2:69],
                                                                  SFS_S_2_W = GC_to_AT_ALL_SFSes.ZI69.SI.A.mean.5GCbins[[i]][2:69],
                                                                  SFS_W_2_S = AT_to_GC_ALL_SFSes.ZI69.SI.A.mean.5GCbins[[i]][2:69],
                                                                  GC = GCcontent.Dmel.ZI69.SI.A.concat.mean.5GCbins[i],
                                                                  algorithmToUse = "Nelder-Mead",
                                                                  maxit=1000,
                                                                  MYtheta_init0=x$theta,MYthetaWS_init0=x$thetaWS,MYbias_init0=x$bias,MYerr_init0=NULL,MYerrWS_init0=NULL,MYerrSW_init0=NULL, #initial parameters for model0
                                                                  MYtheta_init1=y$theta,MYthetaWS_init1=y$thetaWS,MYbias_init1=y$bias,MYB_init1=y$B,MYerr_init1=NULL,MYerrWS_init1=NULL,MYerrSW_init1=NULL)
  })

results_ZI69.SI.A.5GCbins.mel_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    # this is just to make the arguments list in the next part less long:
    x <- results_ZI69.SI.A.5GCbins.mel_no_error[[i]]$without_error$model0
    y <- results_ZI69.SI.A.5GCbins.mel_no_error[[i]]$without_error$model1
    
    results_ZI69.SI.A.5GCbins.mel_error[[i]] <- get_test_results_error(SFS_Neu = N_to_N_ALL_SFSes.ZI69.SI.A.mel.5GCbins[[i]][2:69],
                                                                        SFS_S_2_W = GC_to_AT_ALL_SFSes.ZI69.SI.A.mel.5GCbins[[i]][2:69],
                                                                        SFS_W_2_S = AT_to_GC_ALL_SFSes.ZI69.SI.A.mel.5GCbins[[i]][2:69],
                                                                        GC = GCcontent.Dmel.ZI69.SI.A.concat.mel.5GCbins[i],
                                                                        algorithmToUse = "Nelder-Mead",
                                                                        maxit=1000,
                                                                        MYtheta_init0=x$theta,MYthetaWS_init0=x$thetaWS,MYbias_init0=x$bias,MYerr_init0=NULL,MYerrWS_init0=NULL,MYerrSW_init0=NULL, #initial parameters for model0
                                                                        MYtheta_init1=y$theta,MYthetaWS_init1=y$thetaWS,MYbias_init1=y$bias,MYB_init1=y$B,MYerr_init1=NULL,MYerrWS_init1=NULL,MYerrSW_init1=NULL)
  })

results_ZI21.SI.A.5GCbins.mean_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    # this is just to make the arguments list in the next part less long:
    x <- results_ZI21.SI.A.5GCbins.mean_no_error[[i]]$without_error$model0
    y <- results_ZI21.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1
    
    results_ZI21.SI.A.5GCbins.mean_error[[i]] <- get_test_results_error(SFS_Neu = N_to_N_ALL_SFSes.ZI21.SI.A.mean.5GCbins[[i]][2:21],
                                                                   SFS_S_2_W = GC_to_AT_ALL_SFSes.ZI21.SI.A.mean.5GCbins[[i]][2:21],
                                                                   SFS_W_2_S = AT_to_GC_ALL_SFSes.ZI21.SI.A.mean.5GCbins[[i]][2:21],
                                                                   GC = GCcontent.Dmel.ZI21.SI.A.concat.mean.5GCbins[i],
                                                                   algorithmToUse = "Nelder-Mead",
                                                                   maxit=1000,
                                                                   MYtheta_init0=x$theta,MYthetaWS_init0=x$thetaWS,MYbias_init0=x$bias,MYerr_init0=NULL,MYerrWS_init0=NULL,MYerrSW_init0=NULL, #initial parameters for model0
                                                                   MYtheta_init1=y$theta,MYthetaWS_init1=y$thetaWS,MYbias_init1=y$bias,MYB_init1=y$B,MYerr_init1=NULL,MYerrWS_init1=NULL,MYerrSW_init1=NULL)
  })

results_ZI21.SI.A.5GCbins.mel_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    # this is just to make the arguments list in the next part less long:
    x <- results_ZI21.SI.A.5GCbins.mel_no_error[[i]]$without_error$model0
    y <- results_ZI21.SI.A.5GCbins.mel_no_error[[i]]$without_error$model1
    
    results_ZI21.SI.A.5GCbins.mel_error[[i]] <- get_test_results_error(SFS_Neu = N_to_N_ALL_SFSes.ZI21.SI.A.mel.5GCbins[[i]][2:21],
                                                                        SFS_S_2_W = GC_to_AT_ALL_SFSes.ZI21.SI.A.mel.5GCbins[[i]][2:21],
                                                                        SFS_W_2_S = AT_to_GC_ALL_SFSes.ZI21.SI.A.mel.5GCbins[[i]][2:21],
                                                                        GC = GCcontent.Dmel.ZI21.SI.A.concat.mel.5GCbins[i],
                                                                        algorithmToUse = "Nelder-Mead",
                                                                        maxit=1000,
                                                                        MYtheta_init0=x$theta,MYthetaWS_init0=x$thetaWS,MYbias_init0=x$bias,MYerr_init0=NULL,MYerrWS_init0=NULL,MYerrSW_init0=NULL, #initial parameters for model0
                                                                        MYtheta_init1=y$theta,MYthetaWS_init1=y$thetaWS,MYbias_init1=y$bias,MYB_init1=y$B,MYerr_init1=NULL,MYerrWS_init1=NULL,MYerrSW_init1=NULL)
  })


results_RG21.SI.A.5GCbins.mean_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    # this is just to make the arguments list in the next part less long:
    x <- results_RG21.SI.A.5GCbins.mean_no_error[[i]]$without_error$model0
    y <- results_RG21.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1
    
    results_RG21.SI.A.5GCbins.mean_error[[i]] <- get_test_results_error(SFS_Neu = N_to_N_ALL_SFSes.RG21.SI.A.mean.5GCbins[[i]][2:21],
                                                                   SFS_S_2_W = GC_to_AT_ALL_SFSes.RG21.SI.A.mean.5GCbins[[i]][2:21],
                                                                   SFS_W_2_S = AT_to_GC_ALL_SFSes.RG21.SI.A.mean.5GCbins[[i]][2:21],
                                                                   GC = GCcontent.Dmel.RG21.SI.A.concat.mean.5GCbins[i],
                                                                   algorithmToUse = "Nelder-Mead",
                                                                   maxit=1000,
                                                                   MYtheta_init0=x$theta,MYthetaWS_init0=x$thetaWS,MYbias_init0=x$bias,MYerr_init0=NULL,MYerrWS_init0=NULL,MYerrSW_init0=NULL, #initial parameters for model0
                                                                   MYtheta_init1=y$theta,MYthetaWS_init1=y$thetaWS,MYbias_init1=y$bias,MYB_init1=y$B,MYerr_init1=NULL,MYerrWS_init1=NULL,MYerrSW_init1=NULL)
  })

results_RG21.SI.A.5GCbins.mel_error <- vector('list', length = 5)
system.time(
  for(i in 1:5){
    # this is just to make the arguments list in the next part less long:
    x <- results_RG21.SI.A.5GCbins.mel_no_error[[i]]$without_error$model0
    y <- results_RG21.SI.A.5GCbins.mel_no_error[[i]]$without_error$model1
    
    results_RG21.SI.A.5GCbins.mel_error[[i]] <- get_test_results_error(SFS_Neu = N_to_N_ALL_SFSes.RG21.SI.A.mel.5GCbins[[i]][2:21],
                                                                        SFS_S_2_W = GC_to_AT_ALL_SFSes.RG21.SI.A.mel.5GCbins[[i]][2:21],
                                                                        SFS_W_2_S = AT_to_GC_ALL_SFSes.RG21.SI.A.mel.5GCbins[[i]][2:21],
                                                                        GC = GCcontent.Dmel.RG21.SI.A.concat.mel.5GCbins[i],
                                                                        algorithmToUse = "Nelder-Mead",
                                                                        maxit=1000,
                                                                        MYtheta_init0=x$theta,MYthetaWS_init0=x$thetaWS,MYbias_init0=x$bias,MYerr_init0=NULL,MYerrWS_init0=NULL,MYerrSW_init0=NULL, #initial parameters for model0
                                                                        MYtheta_init1=y$theta,MYthetaWS_init1=y$thetaWS,MYbias_init1=y$bias,MYB_init1=y$B,MYerr_init1=NULL,MYerrWS_init1=NULL,MYerrSW_init1=NULL)
  })

# check for convergence:
lapply(results_MD.SI.A.5GCbins.mean_error, FUN = function(x){x$conv_table})
lapply(results_ZI69.SI.A.5GCbins.mean_error, FUN = function(x){x$conv_table})
lapply(results_ZI21.SI.A.5GCbins.mean_error, FUN = function(x){x$conv_table})
lapply(results_RG21.SI.A.5GCbins.mean_error, FUN = function(x){x$conv_table})
lapply(results_MD.SI.A.5GCbins.sim_error, FUN = function(x){x$conv_table})
lapply(results_ZI69.SI.A.5GCbins.mel_error, FUN = function(x){x$conv_table})
lapply(results_ZI21.SI.A.5GCbins.mel_error, FUN = function(x){x$conv_table})
lapply(results_RG21.SI.A.5GCbins.mel_error, FUN = function(x){x$conv_table})

# so save the results
save(results_MD.SI.A.5GCbins.mean_no_error,
     results_ZI69.SI.A.5GCbins.mean_no_error,
     results_ZI21.SI.A.5GCbins.mean_no_error,
     results_RG21.SI.A.5GCbins.mean_no_error,
     results_MD.SI.A.5GCbins.mean_error,
     results_ZI69.SI.A.5GCbins.mean_error,
     results_ZI21.SI.A.5GCbins.mean_error,
     results_RG21.SI.A.5GCbins.mean_error,
     results_MD.SI.A.5GCbins.sim_no_error,
     results_ZI69.SI.A.5GCbins.mel_no_error,
     results_ZI21.SI.A.5GCbins.mel_no_error,
     results_RG21.SI.A.5GCbins.mel_no_error,
     results_MD.SI.A.5GCbins.sim_error,
     results_ZI69.SI.A.5GCbins.mel_error,
     results_ZI21.SI.A.5GCbins.mel_error,
     results_RG21.SI.A.5GCbins.mel_error,
     file = '../data/glemin_results.RData')



#