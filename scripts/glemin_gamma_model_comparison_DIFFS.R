setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load the results
load('../data/glemin_results_DIFFS.RData')

# do the models with gBGC fit better when we don't include polarisation error?
# sim, diff bins - T,F,T,T,T
for(i in 1:5){
  x0 <- results_MD.SI.A.5GCbins.diff_no_error[[i]]$without_error$model0$lnL
  x1 <- results_MD.SI.A.5GCbins.diff_no_error[[i]]$without_error$model1$lnL
  print(pchisq(q = 2*(x1 - x0),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(x0, x1)

# mel, ZI, 69 samples, diff bins
# T,T,T,T,T
for(i in 1:5){
  x0 <- results_ZI69.SI.A.5GCbins.diff_no_error[[i]]$without_error$model0$lnL
  x1 <- results_ZI69.SI.A.5GCbins.diff_no_error[[i]]$without_error$model1$lnL
  print(pchisq(q = 2*(x1 - x0),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(x0, x1)

###############################
# do the models with gBGC fit better when we include polarisation error?
# sim diff - T,F,T,T,T
for(i in 1:5){
  m0star <- results_MD.SI.A.5GCbins.diff_error[[i]]$with_error$model0$lnL
  m1star <- results_MD.SI.A.5GCbins.diff_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m0star),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(m0star, m1star)

# mel, ZI, 69 samples, diff bins
# T,T,T,T,T
for(i in 1:5){
  m0star <- results_ZI69.SI.A.5GCbins.diff_error[[i]]$with_error$model0$lnL
  m1star <- results_ZI69.SI.A.5GCbins.diff_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m0star),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(m0star, m1star)

##############################
# do the models with polarisation error and gBGC fit better than models with
# no polarisation error and gBGC?
# sim diff- No
for(i in 1:5){
  m1 <- results_MD.SI.A.5GCbins.diff_no_error[[i]]$without_error$model1$lnL
  m1star <- results_MD.SI.A.5GCbins.diff_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m1),
               df = 3,
               lower.tail = F) < 0.05)
}
rm(m1, m1star)


# mel, ZI, 69 lines diff - No
for(i in 1:5){
  m1 <- results_ZI69.SI.A.5GCbins.diff_no_error[[i]]$without_error$model1$lnL
  m1star <- results_ZI69.SI.A.5GCbins.diff_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m1),
               df = 3,
               lower.tail = F) < 0.05)
}
rm(m1, m1star)


