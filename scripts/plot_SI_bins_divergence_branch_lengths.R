# plotting new analyses

# load the bin indices and GC contents
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_binindices_GCcontent.RData')
source('~/Dropbox/Biology_projects/dros_X_A_mut/scripts/parse_PAML_output_file.R')

# load the bootstrap data
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/PAML_bootstraps.RData')

myOutputFiles.mean <- list.files(path = '~/Dropbox/Biology_projects/dros_X_A_mut/PAML/mean_SI_A_GCbins/output_files_GC_bins',
                                 pattern = 'output_GC_bin', full.names = T)
myOutputFiles.mel <- list.files(path = '~/Dropbox/Biology_projects/dros_X_A_mut/PAML/mel_SI_A_GCbins/output_files_GC_bins',
                                pattern = 'output_GC_bin', full.names = T)
myOutputFiles.sim <- list.files(path = '~/Dropbox/Biology_projects/dros_X_A_mut/PAML/sim_SI_A_GCbins/output_files_GC_bins',
                                pattern = 'output_GC_bin', full.names = T)

results.mean <- lapply(myOutputFiles.mean, readOutput)
results.mel <- lapply(myOutputFiles.mel, readOutput)
results.sim <- lapply(myOutputFiles.sim, readOutput)

blength.mean.mel.obs <- sapply(results.mean, FUN = function(x){x$mel_blength})
blength.mean.sim.obs <- sapply(results.mean, FUN = function(x){x$sim_blength})
blength.mel.mel.obs <- sapply(results.mel, FUN = function(x){x$mel_blength})
blength.sim.sim.obs <- sapply(results.sim, FUN = function(x){x$sim_blength})

blength.mean.sim.CIs <- sapply(results.mean.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_blength
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

blength.mean.mel.CIs <- sapply(results.mean.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_blength
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

blength.sim.sim.CIs <- sapply(results.sim.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_blength
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

blength.mel.mel.CIs <- sapply(results.mel.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_blength
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

pdf('~/Dropbox/Biology_projects/dros_X_A_mut/plots/blength_SI_A_sim_meanbins_ALT.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = blength.mean.sim.obs,
     xlab = 'GC content',
     ylab = 'branch length',
     main = 'simulans (mean bins)',
     ylim = c(0, 0.12),
     axes = F
)

axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.meanbins.SI.A),
     las = 2)
axis(2)
arrows(x0 = 1:5,
       y0 = blength.mean.sim.CIs[1,],
       y1 = blength.mean.sim.CIs[2,],
       code = 3,
       angle = 90,
       length = 0.1)
abline(h = 0.058167, col = 'red')
box()
dev.off()

pdf('~/Dropbox/Biology_projects/dros_X_A_mut/plots/blength_SI_A_mel_meanbins_ALT.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = blength.mean.mel.obs,
     xlab = 'GC content',
     ylab = 'branch length',
     main = 'melanogaster (mean bins)',
     ylim = c(0, 0.12),
     axes = F
)

axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.meanbins.SI.A),
     las = 2)
arrows(x0 = 1:5,
       y0 = blength.mean.mel.CIs[1,],
       y1 = blength.mean.mel.CIs[2,],
       code = 3,
       angle = 90,
       length = 0.1)
abline(h = 0.06478, col = 'red')
axis(2)
box()
dev.off()


################################################################################

pdf('~/Dropbox/Biology_projects/dros_X_A_mut/plots/blength_SI_A_sim_simbins_ALT.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = blength.sim.sim.obs,
     xlab = 'GC content',
     ylab = 'branch length',
     main = 'simulans (sim bins)',
     ylim = c(0, 0.12),
     axes = F
)

axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.simbins.SI.A),
     las = 2)
axis(2)
arrows(x0 = 1:5,
       y0 = blength.sim.sim.CIs[1,],
       y1 = blength.sim.sim.CIs[2,],
       code = 3,
       angle = 90,
       length = 0.1)
abline(h = 0.058167, col = 'red')
box()
dev.off()

pdf('~/Dropbox/Biology_projects/dros_X_A_mut/plots/blength_SI_A_mel_melbins_ALT.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = blength.mel.mel.obs,
     xlab = 'GC content',
     ylab = 'branch length',
     main = 'melanogaster',
     ylim = c(0, 0.12),
     axes = F
)
# mtext('substitution count ratio',
#       side = 2,
#       line = 1.2,
#       cex = 0.75)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.melbins.SI.A),
     las = 2)
arrows(x0 = 1:5,
       y0 = blength.mel.mel.CIs[1,],
       y1 = blength.mel.mel.CIs[2,],
       code = 3,
       angle = 90,
       length = 0.1)
abline(h = 0.06478, col = 'red')
axis(2)
box()
dev.off()
