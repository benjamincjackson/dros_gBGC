setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load the results
load('../data/glemin_results.RData')

# do the models with gBGC fit better when we don't include polarisation error?
# sim, mean bins - YES
for(i in 1:5){
  x0 <- results_MD.SI.A.5GCbins.mean_no_error[[i]]$without_error$model0$lnL
  x1 <- results_MD.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1$lnL
  print(pchisq(q = 2*(x1 - x0),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(x0, x1)

# sim, sim bins - No, Yes yes yes yes
for(i in 1:5){
  x0 <- results_MD.SI.A.5GCbins.sim_no_error[[i]]$without_error$model0$lnL
  x1 <- results_MD.SI.A.5GCbins.sim_no_error[[i]]$without_error$model1$lnL
  print(pchisq(q = 2*(x1 - x0),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(x0, x1)

# mel, ZI, 69 samples, mean bins
# No for the first two bins, yes for the last three
for(i in 1:5){
  x0 <- results_ZI69.SI.A.5GCbins.mean_no_error[[i]]$without_error$model0$lnL
  x1 <- results_ZI69.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1$lnL
  print(pchisq(q = 2*(x1 - x0),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(x0, x1)

# mel, ZI, 69 samples, mel bins
# yes, no, yes, yes, yes
for(i in 1:5){
  x0 <- results_ZI69.SI.A.5GCbins.mel_no_error[[i]]$without_error$model0$lnL
  x1 <- results_ZI69.SI.A.5GCbins.mel_no_error[[i]]$without_error$model1$lnL
  print(pchisq(q = 2*(x1 - x0),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(x0, x1)

# mel, ZI, 21 samples
# No for the first bin, yes for the last four
for(i in 1:5){
  x0 <- results_ZI21.SI.A.5GCbins.mean_no_error[[i]]$without_error$model0$lnL
  x1 <- results_ZI21.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1$lnL
  print(pchisq(q = 2*(x1 - x0),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(x0, x1)

# mel, ZI, 21 samples
# yes, no, yes, yes, yes
for(i in 1:5){
  x0 <- results_ZI21.SI.A.5GCbins.mel_no_error[[i]]$without_error$model0$lnL
  x1 <- results_ZI21.SI.A.5GCbins.mel_no_error[[i]]$without_error$model1$lnL
  print(pchisq(q = 2*(x1 - x0),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(x0, x1)

# mel, RG, 21 samples mean bins
# No for the first two bins, yes for the last three
for(i in 1:5){
  x0 <- results_RG21.SI.A.5GCbins.mean_no_error[[i]]$without_error$model0$lnL
  x1 <- results_RG21.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1$lnL
  print(pchisq(q = 2*(x1 - x0),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(x0, x1)

# mel, RG, 21 samples mel bins
# yes, no, yes yes yes
for(i in 1:5){
  x0 <- results_RG21.SI.A.5GCbins.mel_no_error[[i]]$without_error$model0$lnL
  x1 <- results_RG21.SI.A.5GCbins.mel_no_error[[i]]$without_error$model1$lnL
  print(pchisq(q = 2*(x1 - x0),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(x0, x1)

###############################
# do the models with gBGC fit better when we include polarisation error?
# sim mean - no for the first bin, yes for the last 4
for(i in 1:5){
  m0star <- results_MD.SI.A.5GCbins.mean_error[[i]]$with_error$model0$lnL
  m1star <- results_MD.SI.A.5GCbins.mean_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m0star),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(m0star, m1star)

# sim sim - no for the first bin, yes for the last 4
for(i in 1:5){
  m0star <- results_MD.SI.A.5GCbins.sim_error[[i]]$with_error$model0$lnL
  m1star <- results_MD.SI.A.5GCbins.sim_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m0star),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(m0star, m1star)

# mel, ZI, 69 samples, mean bins
# No for the first two bins, yes for the last three
for(i in 1:5){
  m0star <- results_ZI69.SI.A.5GCbins.mean_error[[i]]$with_error$model0$lnL
  m1star <- results_ZI69.SI.A.5GCbins.mean_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m0star),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(m0star, m1star)

# mel, ZI, 69 samples, mel bins
# yes, no, yes yes yes
for(i in 1:5){
  m0star <- results_ZI69.SI.A.5GCbins.mel_error[[i]]$with_error$model0$lnL
  m1star <- results_ZI69.SI.A.5GCbins.mel_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m0star),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(m0star, m1star)

# mel, ZI, 21 samples mean bins
# No, no, yes, yes, no
for(i in 1:5){
  m0star <- results_ZI21.SI.A.5GCbins.mean_error[[i]]$with_error$model0$lnL
  m1star <- results_ZI21.SI.A.5GCbins.mean_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m0star),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(m0star, m1star)

# mel, ZI, 21 samples mel bins
# yes, no, yes yes yes
for(i in 1:5){
  m0star <- results_ZI21.SI.A.5GCbins.mel_error[[i]]$with_error$model0$lnL
  m1star <- results_ZI21.SI.A.5GCbins.mel_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m0star),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(m0star, m1star)

# mel, RG, 21 samples mean bins
# no, no, no, yes, yes
for(i in 1:5){
  m0star <- results_RG21.SI.A.5GCbins.mean_error[[i]]$with_error$model0$lnL
  m1star <- results_RG21.SI.A.5GCbins.mean_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m0star),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(m0star, m1star)

# mel, RG, 21 samples
# Yes, no, yes, yes, yes
for(i in 1:5){
  m0star <- results_RG21.SI.A.5GCbins.mel_error[[i]]$with_error$model0$lnL
  m1star <- results_RG21.SI.A.5GCbins.mel_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m0star),
               df = 1,
               lower.tail = F) < 0.05)
}
rm(m0star, m1star)
##############################
# do the models with polarisation error and gBGC fit better than models with
# no polarisation error and gBGC?
# sim mean- No
for(i in 1:5){
  m1 <- results_MD.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1$lnL
  m1star <- results_MD.SI.A.5GCbins.mean_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m1),
               df = 3,
               lower.tail = F) < 0.05)
}
rm(m1, m1star)

# sim sim- No
for(i in 1:5){
  m1 <- results_MD.SI.A.5GCbins.sim_no_error[[i]]$without_error$model1$lnL
  m1star <- results_MD.SI.A.5GCbins.sim_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m1),
               df = 3,
               lower.tail = F) < 0.05)
}
rm(m1, m1star)

# mel, ZI, 69 lines mean - No
for(i in 1:5){
  m1 <- results_ZI69.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1$lnL
  m1star <- results_ZI69.SI.A.5GCbins.mean_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m1),
               df = 3,
               lower.tail = F) < 0.05)
}
rm(m1, m1star)

# mel, ZI, 69 lines mel - No
for(i in 1:5){
  m1 <- results_ZI69.SI.A.5GCbins.mel_no_error[[i]]$without_error$model1$lnL
  m1star <- results_ZI69.SI.A.5GCbins.mel_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m1),
               df = 3,
               lower.tail = F) < 0.05)
}
rm(m1, m1star)

# mel, ZI, 21 samples
# No
for(i in 1:5){
  m1 <- results_ZI21.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1$lnL
  m1star <- results_ZI21.SI.A.5GCbins.mean_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m1),
               df = 3,
               lower.tail = F) < 0.05)
}
rm(m1, m1star)

# mel, ZI, 21 samples mel
# No
for(i in 1:5){
  m1 <- results_ZI21.SI.A.5GCbins.mel_no_error[[i]]$without_error$model1$lnL
  m1star <- results_ZI21.SI.A.5GCbins.mel_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m1),
               df = 3,
               lower.tail = F) < 0.05)
}
rm(m1, m1star)

# mel, RG, 21 samples mean - No
for(i in 1:5){
  m1 <- results_RG21.SI.A.5GCbins.mean_no_error[[i]]$without_error$model1$lnL
  m1star <- results_RG21.SI.A.5GCbins.mean_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m1),
               df = 3,
               lower.tail = F) < 0.05)
}
rm(m1, m1star)

# mel, RG, 21 samples mel - No
for(i in 1:5){
  m1 <- results_RG21.SI.A.5GCbins.mel_no_error[[i]]$without_error$model1$lnL
  m1star <- results_RG21.SI.A.5GCbins.mel_error[[i]]$with_error$model1$lnL
  print(pchisq(q = 2*(m1star - m1),
               df = 3,
               lower.tail = F) < 0.05)
}
rm(m1, m1star)

