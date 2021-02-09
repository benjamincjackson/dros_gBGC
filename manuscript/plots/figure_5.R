setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# plotting figure 5

# load the bin indices and GC contents
load('../../data/SI_binindices_GCcontent.RData')
source('../../scripts/PAML_parse_output_file.R')

# load the bootstrap data
load('../../data/PAML_bootstraps.RData')

myOutputFiles.mean <- list.files(path = '../../PAML/mean_SI_A_GCbins/output_files_GC_bins',
                                 pattern = 'output_GC_bin', full.names = T)
myOutputFiles.mel <- list.files(path = '../../PAML/mel_SI_A_GCbins/output_files_GC_bins',
                                pattern = 'output_GC_bin', full.names = T)
myOutputFiles.sim <- list.files(path = '../../PAML/sim_SI_A_GCbins/output_files_GC_bins',
                                pattern = 'output_GC_bin', full.names = T)

results.mean <- lapply(myOutputFiles.mean, readOutput)
results.mel <- lapply(myOutputFiles.mel, readOutput)
results.sim <- lapply(myOutputFiles.sim, readOutput)

kappa.mean.mel.obs <- sapply(results.mean, FUN = function(x){x$mel_kappa})
kappa.mean.sim.obs <- sapply(results.mean, FUN = function(x){x$sim_kappa})
kappa.mel.mel.obs <- sapply(results.mel, FUN = function(x){x$mel_kappa})
kappa.sim.sim.obs <- sapply(results.sim, FUN = function(x){x$sim_kappa})

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


pdf('figure_5.pdf',
    width = 10,
    height = 10)
par(lwd = 2,
    cex = 8,
    mfrow = c(2,2),
    oma = c(2,2,1,1))
plot(x = 1:5,
     y = kappa.mean.sim.obs,
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(0, 10),
     axes = F,
     type = 'n')
text(x = 1.1, y = 9.6,
     labels = c("A"),
     cex = 3)
points(x = 1:5,
       y = kappa.mean.sim.obs,
       cex = 2)
title(main = 'simulans',
      line = 1,
      cex.main = 2)
mtext('substitution rate ratio',
      side = 2,
      line = 3,
      cex = 1.3)
mtext('GC content',
      side = 1,
      line = 4,
      cex = 1.3)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.sim.meanbins.SI.A),
     las = 2, cex.axis = 1.4)
axis(2, cex.axis = 1.4)
box()
arrows(x0 = 1:5,
       y0 = kappa.mean.sim.CIs[1,],
       y1 = kappa.mean.sim.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1)


plot(x = 1:5,
     y = kappa.mean.mel.obs,
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(0, 10),
     axes = F,
     type = 'n'
)
text(x = 1.1, y = 9.6,
     labels = c("B"),
     cex = 3)
points(x = 1:5,
       y = kappa.mean.mel.obs,
       cex = 2)
title(main = 'melanogaster',
      line = 1,
      cex.main = 2)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2, cex.axis = 1.4)
axis(2, cex.axis = 1.4)
box()
legend('topright',
       legend = expression(r[S %->% W]/r[W %->% S]),
       pch = '',
       bty = 'n',
       cex = 2)
mtext('GC content',
      side = 1,
      line = 4,
      cex = 1.3)
arrows(x0 = 1:5,
       y0 = kappa.mean.mel.CIs[1,],
       y1 = kappa.mean.mel.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1)


plot(x = 1:5,
     y = kappa.sim.sim.obs,
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(0, 10),
     axes = F,
     type = 'n'
)
text(x = 4.9, y = 9.6,
     labels = c("C"),
     cex = 3)
points(x = 1:5,
       y = kappa.sim.sim.obs,
       cex = 2)
title(main = 'simulans',
      line = 1,
      cex.main = 2)
mtext('substitution count ratio',
      side = 2,
      line = 3,
      cex = 1.3)
mtext('GC content',
      side = 1,
      line = 4,
      cex = 1.3)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.sim.simbins.SI.A),
     las = 2, cex.axis = 1.4)
axis(2, cex.axis = 1.4)
box()
arrows(x0 = 1:5,
       y0 = kappa.sim.sim.CIs[1,],
       y1 = kappa.sim.sim.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1)


plot(x = 1:5,
     y = kappa.mel.mel.obs,
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(0, 10),
     axes = F,
     type = 'n'
)
text(x = 4.9, y = 9.6,
     labels = c("D"),
     cex = 3)
points(x = 1:5,
       y = kappa.mel.mel.obs,
       cex = 2)
title(main = 'melanogaster',
      line = 1,
      cex.main = 2)
mtext('GC content',
      side = 1,
      line = 4,
      cex = 1.3)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.mel.melbins.SI.A),
     las = 2, cex.axis = 1.4)
axis(2, cex.axis = 1.4)
box()
arrows(x0 = 1:5,
       y0 = kappa.mel.mel.CIs[1,],
       y1 = kappa.mel.mel.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1)
dev.off()

