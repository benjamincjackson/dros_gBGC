setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# plotting new analyses

# load the bin indices and GC contents
load('../data/SI_binindices_GCcontent.RData')
load('../data/SI_binindices_GCcontent_DIFFERENCES.RData')
source('../scripts/PAML_parse_output_file.R')

# load the bootstrap data
load('../data/PAML_bootstraps_noPoly.RData')

myOutputFiles.mean <- list.files(path = '../PAML/no_poly/mean_SI_A_GCbins/output_files_GC_bins',
                              pattern = 'output_GC_bin', full.names = T)
myOutputFiles.diff <- list.files(path = '../PAML/no_poly/diff_SI_A_GCbins/output_files_GC_bins',
                                 pattern = 'output_GC_bin', full.names = T)
myOutputFiles.mel <- list.files(path = '../PAML/no_poly/mel_SI_A_GCbins/output_files_GC_bins',
                                pattern = 'output_GC_bin', full.names = T)
myOutputFiles.sim <- list.files(path = '../PAML/no_poly/sim_SI_A_GCbins/output_files_GC_bins',
                                pattern = 'output_GC_bin', full.names = T)

results.mean <- lapply(myOutputFiles.mean, readOutput)
results.diff <- lapply(myOutputFiles.diff, readOutput)
results.mel <- lapply(myOutputFiles.mel, readOutput)
results.sim <- lapply(myOutputFiles.sim, readOutput)

kappa.mean.mel.obs <- sapply(results.mean, FUN = function(x){x$mel_kappa})
kappa.mean.sim.obs <- sapply(results.mean, FUN = function(x){x$sim_kappa})
kappa.diff.mel.obs <- sapply(results.diff, FUN = function(x){x$mel_kappa})
kappa.diff.sim.obs <- sapply(results.diff, FUN = function(x){x$sim_kappa})
kappa.mel.mel.obs <- sapply(results.mel, FUN = function(x){x$mel_kappa})
kappa.sim.sim.obs <- sapply(results.sim, FUN = function(x){x$sim_kappa})
# also the ratio down the opposing lineage
kappa.mel.sim.obs <- sapply(results.sim, FUN = function(x){x$mel_kappa})
kappa.sim.mel.obs <- sapply(results.mel, FUN = function(x){x$sim_kappa})

kappa.mean.sim.CIs <- sapply(results.mean.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_kappa
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

kappa.mean.mel.CIs <- sapply(results.mean.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_kappa
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

kappa.diff.sim.CIs <- sapply(results.diff.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_kappa
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

kappa.diff.mel.CIs <- sapply(results.diff.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_kappa
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

kappa.sim.sim.CIs <- sapply(results.sim.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_kappa
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

kappa.mel.mel.CIs <- sapply(results.mel.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_kappa
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

pdf('../plots/ratio_of_rates_SI_A_sim_meanbins_noPoly.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = kappa.mean.sim.obs,
     xlab = 'GC content',
     ylab = 'substitution rate ratio',
     main = 'simulans (mean bins,\n polymorphism excluded)',
     ylim = c(0, 10),
     axes = F
)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.sim.meanbins.SI.A),
     las = 2)
axis(2)
arrows(x0 = 1:5,
       y0 = kappa.mean.sim.CIs[1,],
       y1 = kappa.mean.sim.CIs[2,],
       code = 3,
       angle = 90,
       length = 0.1)
box()
dev.off()

pdf('../plots/ratio_of_rates_SI_A_sim_diffbins_noPoly.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = kappa.diff.sim.obs,
     xlab = 'GC content',
     ylab = 'substitution rate ratio',
     main = 'simulans (diff bins,\n polymorphism excluded)',
     ylim = c(0, 10),
     axes = F
)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      signif(aggregate(GCcontent.diff, by = list(GCbinindex.5bins.diff.SI.A), FUN = mean)$x,
                             2)),
     las = 2)
axis(2)
arrows(x0 = 1:5,
       y0 = kappa.diff.sim.CIs[1,],
       y1 = kappa.diff.sim.CIs[2,],
       code = 3,
       angle = 90,
       length = 0.1)
box()
dev.off()

pdf('../plots/ratio_of_rates_SI_A_mel_diffbins_noPoly.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = kappa.diff.mel.obs,
     xlab = expression(GC[diff]),
     ylab = 'substitution rate ratio',
     main = 'simulans (diff bins,\n polymorphism excluded)',
     ylim = c(0, 10),
     axes = F
)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      signif(aggregate(GCcontent.diff, by = list(GCbinindex.5bins.diff.SI.A), FUN = mean)$x,
                             2)),
     las = 2)
axis(2)
arrows(x0 = 1:5,
       y0 = kappa.diff.mel.CIs[1,],
       y1 = kappa.diff.mel.CIs[2,],
       code = 3,
       angle = 90,
       length = 0.1)
box()
dev.off()

pdf('../plots/ratio_of_rates_SI_A_mel_meanbins_noPoly.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = kappa.mean.mel.obs,
     xlab = 'GC content',
     ylab = 'substitution rate ratio',
     main = 'melanogaster (mean bins,\n polymorphism excluded)',
     ylim = c(0, 10),
     axes = F
)

axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2)
arrows(x0 = 1:5,
       y0 = kappa.mean.mel.CIs[1,],
       y1 = kappa.mean.mel.CIs[2,],
       code = 3,
       angle = 90,
       length = 0.1)
axis(2)
box()
dev.off()


################################################################################

pdf('../plots/ratio_of_rates_SI_A_sim_simbins_noPoly.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = kappa.sim.sim.obs,
     xlab = 'GC content',
     ylab = 'substitution rate ratio',
     main = 'simulans (sim bins,\n polymorphism excluded)',
     ylim = c(0, 10),
     axes = F
)

axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.sim.simbins.SI.A),
     las = 2)
axis(2)
arrows(x0 = 1:5,
       y0 = kappa.sim.sim.CIs[1,],
       y1 = kappa.sim.sim.CIs[2,],
       code = 3,
       angle = 90,
       length = 0.1)
box()
dev.off()

pdf('../plots/ratio_of_rates_SI_A_mel_melbins_noPoly.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = kappa.mel.mel.obs,
     xlab = 'GC content',
     ylab = 'substitution rate ratio',
     main = 'melanogaster (mel bins,\n polymorphism excluded)',
     ylim = c(0, 10),
     axes = F
)
# mtext('substitution count ratio',
#       side = 2,
#       line = 1.2,
#       cex = 0.75)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.mel.melbins.SI.A),
     las = 2)
arrows(x0 = 1:5,
       y0 = kappa.mel.mel.CIs[1,],
       y1 = kappa.mel.mel.CIs[2,],
       code = 3,
       angle = 90,
       length = 0.1)
axis(2)
box()
dev.off()

# ######################################################################################
# 
# pdf('~/Dropbox/Biology_projects/dros_X_A_mut/plots/ratio_of_rates_SI_A_sim_melbins_noPoly.pdf')
# par(lwd = 2,
#     cex = 1.6)
# plot(x = 1:5,
#      y = kappa.sim.mel.obs,
#      xlab = 'GC content',
#      ylab = 'substitution rate ratio',
#      main = 'simulans (mel bins,\n polymorphism excluded)',
#      ylim = c(0, 10),
#      axes = F
# )
# 
# axis(1, at = 1:5,
#      labels = sprintf("%.2f",
#                       GCcontent.concat.5bins.melbins.SI.A),
#      las = 2)
# axis(2)
# box()
# dev.off()
# 
# pdf('~/Dropbox/Biology_projects/dros_X_A_mut/plots/ratio_of_rates_SI_A_mel_simbins_noPoly.pdf')
# par(lwd = 2,
#     cex = 1.6)
# plot(x = 1:5,
#      y = kappa.mel.sim.obs,
#      xlab = 'GC content',
#      ylab = 'substitution rate ratio',
#      main = 'melanogaster (sim bins,\n polymorphism excluded)',
#      ylim = c(0, 10),
#      axes = F
# )
# # mtext('substitution count ratio',
# #       side = 2,
# #       line = 1.2,
# #       cex = 0.75)
# axis(1, at = 1:5,
#      labels = sprintf("%.2f",
#                       GCcontent.concat.5bins.simbins.SI.A),
#      las = 2)
# axis(2)
# box()
# dev.off()