# carry out the glemin gamma model fitting and comparison for simulans, 
# for each bin in turn.

# need to load the SFSes and the GC content of the bins:
load('SFSes.RData')
load('GC_contents.RData')

# The Gl?min method requires the R libraries "HMP" and "VGAM" 
# install.packages('HMP')
# install.packages('VGAM')

library(HMP)
library(VGAM)

# need to source the 'optimizeBGC' functions, and others
source('RCode/optimizeBGC.R')
source('RCode/basic_functions.R')
source('RCode/likelihood_Poisson.R')
source('RCode/likelihood_Error_Poisson.R')

# want to write a function to take the 3 input SFSes and the 1 GC content and perform
# the tests, and save all the data, probably just in a list
get_test_results <- function(SFS_Neu, SFS_S_2_W, SFS_W_2_S, GC){
  
  result_no_error <- optimizeBGC(dataNeutral = SFS_Neu,
                                 dataWS = SFS_W_2_S,
                                 dataSW = SFS_S_2_W,
                                 dataGC = GC,
                                 PrintInit=T)
  
  result_error <- optimizeBGCwith3error(dataNeutral = SFS_Neu,
                                        dataWS = SFS_W_2_S,
                                        dataSW = SFS_S_2_W,
                                        dataGC = GC,
                                        verbose=1,
                                        MAXIT=500)
  
  # MAKE A CONVERGENCE TABLE, TO CHECK WHAT HAS/HASN'T WORKED MORE EASILY?
  conv_table <- matrix(data = c(result_no_error$model0$convergence,
                                result_no_error$model1$convergence,
                                result_error$model0$convergence,
                                result_error$model1$convergence),
                       nrow = 2, ncol = 2,
                       dimnames = list(c('M_0', 'M_1'),
                                       c('no_error', 'error')))
  
  lnlk_table <- matrix(data = c(result_no_error$model0$lnL,
                                result_no_error$model1$lnL,
                                result_error$model0$lnL,
                                result_error$model1$lnL),
                       nrow = 2, ncol = 2,
                       dimnames = list(c('M_0', 'M_1'),
                                       c('no_error', 'error')))
  
  #return the full results, but also the LnLks in a more convenient form?
  return(list('without_error' = result_no_error,
              'with_error' = result_error,
              'lnlk_table' = lnlk_table,
              'conv_table' = conv_table))
}

# now apply this function over all the bins:
results_4fold_A <- vector('list', length = 20)
system.time(
#   for(i in 1:20){ # DEAR SYLVAIN, USUALLY I WOULD RUN THIS LINE, NOT THE ONE BELOW, BUT IT TAKES LONGER...
  for(i in 1){
    results_4fold_A[[i]] <- get_test_results(SFS_Neu = SFSes.a[[i]]$SFS_N_2_N,
                                             SFS_S_2_W = SFSes.a[[i]]$SFS_S_2_W,
                                             SFS_W_2_S = SFSes.a[[i]]$SFS_W_2_S,
                                             GC = gc.content.sim.ref.mel.bins.a[i])
  })# takes 200s on APS desktop

lapply(results_4fold_A, FUN = function(x) x$conv_table )
