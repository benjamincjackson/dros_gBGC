setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# compare ln-likelihoods for stationary versus non stationary models

# a function to get the ln-likelihood from the output file of PAML
#a function to get out one sequence (observed or reconstructed) from one input file:
getLnLk <- function(inputFileString){
  myCon <- file(inputFileString)
  myFile <- readLines(con = myCon)
  close(myCon)
  
  # the lnlikelihood is on a line that looks like this [note the parameters]:
  # lnL(ntime:  2  np: 10): -144715.932743    +0.000000       ### (equilibrium, rooted [INCORRECT])
  # lnL(ntime:  4  np: 39): -144538.305890    +0.000000       ### (non-equilibrium)
  # lnL(ntime:  3  np: 11): -144704.594399    +0.000000       ### (equilibrium, UNrooted [CORRECT])
  
  # find the line
  lnL_line_number <- grep('ntime', myFile)
  # get the line
  lnL_line <- myFile[[lnL_line_number]]
  # split the line on white space
  lnL_clean <- strsplit(lnL_line, split = '\\s+')
  # get the log likelihood
  lnL <- as.numeric(lnL_clean[[1]][5])
  return(lnL)
}

lnlk.mel.noneq <- vector('numeric', 5)
for(i in 1:5){
  lnlk.mel.noneq[i] <- getLnLk(paste0('../PAML/mel_SI_A_GCbins/output_files_GC_bins/output_GC_bin_0', i, collapse = ''))
}

lnlk.sim.noneq <- vector('numeric', 5)
for(i in 1:5){
  lnlk.sim.noneq[i] <- getLnLk(paste0('../PAML/sim_SI_A_GCbins/output_files_GC_bins/output_GC_bin_0', i, collapse = ''))
}

lnlk.mean.noneq <- vector('numeric', 5)
for(i in 1:5){
  lnlk.mean.noneq[i] <- getLnLk(paste0('../PAML/mean_SI_A_GCbins/output_files_GC_bins/output_GC_bin_0', i, collapse = ''))
}


lnlk.mel.eq <- vector('numeric', 5)
for(i in 1:5){
  lnlk.mel.eq[i] <- getLnLk(paste0('../PAML/eq_check/mel_SI_A_GCbins/bin_0', i, '/run_1/output', collapse = ''))
}

lnlk.sim.eq <- vector('numeric', 5)
for(i in 1:5){
  lnlk.sim.eq[i] <- getLnLk(paste0('../PAML/eq_check/sim_SI_A_GCbins/bin_0', i, '/run_1/output', collapse = ''))
}

lnlk.mean.eq <- vector('numeric', 5)
for(i in 1:5){
  lnlk.mean.eq[i] <- getLnLk(paste0('../PAML/eq_check/mean_SI_A_GCbins/bin_0', i, '/run_1/output', collapse = ''))
}

## FOR THE UNROOTED (GTR) MODEL VS THE ROOTED (GTR-NH) MODEL,
# print chi squared values and associated p-values
for(i in 1:5){
  print(lnlk.mel.noneq[i] - lnlk.mel.eq[i])
  print(pchisq(2*(lnlk.mel.noneq[i] - lnlk.mel.eq[i]), df = 28, lower.tail = F))
}

for(i in 1:5){
  print(lnlk.sim.noneq[i] - lnlk.sim.eq[i])
  print(pchisq(2*(lnlk.sim.noneq[i] - lnlk.sim.eq[i]), df = 28, lower.tail = F))
}

for(i in 1:5){
  print(lnlk.mean.noneq[i] - lnlk.mean.eq[i])
  print(pchisq(2*(lnlk.mean.noneq[i] - lnlk.mean.eq[i]), df = 28, lower.tail = F))
}
