setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# plotting new analyses

# load the bin indices and GC contents
load('../data/SI_binindices_GCcontent.RData')
source('../scripts/PAML_parse_output_file.R')

# load the bootstrap data
load('../data/PAML_bootstraps.RData')

myOutputFiles.mean <- list.files(path = '../PAML/mean_SI_A_GCbins/output_files_GC_bins',
                                 pattern = 'output_GC_bin', full.names = T)
myOutputFiles.mel <- list.files(path = '../PAML/mel_SI_A_GCbins/output_files_GC_bins',
                                 pattern = 'output_GC_bin', full.names = T)
myOutputFiles.sim <- list.files(path = '../PAML/sim_SI_A_GCbins/output_files_GC_bins',
                                 pattern = 'output_GC_bin', full.names = T)

results.mean <- lapply(myOutputFiles.mean, readOutput)
results.mel <- lapply(myOutputFiles.mel, readOutput)
results.sim <- lapply(myOutputFiles.sim, readOutput)

rates.AT_2_GC.mean.mel.obs <- sapply(results.mean, FUN = function(x){x$mel_AT_2_CG_rate})
rates.GC_2_AT.mean.mel.obs <- sapply(results.mean, FUN = function(x){x$mel_GC_2_AT_rate})
rates.AT_2_GC.mean.sim.obs <- sapply(results.mean, FUN = function(x){x$sim_AT_2_CG_rate})
rates.GC_2_AT.mean.sim.obs <- sapply(results.mean, FUN = function(x){x$sim_GC_2_AT_rate})

rates.GC_2_AT.mel.mel.obs <- sapply(results.mel, FUN = function(x){x$mel_GC_2_AT_rate})
rates.AT_2_GC.mel.mel.obs <- sapply(results.mel, FUN = function(x){x$mel_AT_2_CG_rate})
rates.GC_2_AT.sim.sim.obs <- sapply(results.sim, FUN = function(x){x$sim_GC_2_AT_rate})
rates.AT_2_GC.sim.sim.obs <- sapply(results.sim, FUN = function(x){x$sim_AT_2_CG_rate})

rates.AT_2_GC.mean.mel.CIs <- sapply(results.mean.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_AT_2_CG_rate
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

rates.GC_2_AT.mean.mel.CIs <- sapply(results.mean.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_GC_2_AT_rate
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

rates.AT_2_GC.mean.sim.CIs <- sapply(results.mean.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_AT_2_CG_rate
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

rates.GC_2_AT.mean.sim.CIs <- sapply(results.mean.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_GC_2_AT_rate
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})


rates.AT_2_GC.mel.mel.CIs <- sapply(results.mel.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_AT_2_CG_rate
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

rates.GC_2_AT.mel.mel.CIs <- sapply(results.mel.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_GC_2_AT_rate
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

rates.AT_2_GC.sim.sim.CIs <- sapply(results.sim.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_AT_2_CG_rate
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

rates.GC_2_AT.sim.sim.CIs <- sapply(results.sim.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_GC_2_AT_rate
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})


ratios.mean.sim <- sapply(results.mean, FUN = function(x){x$sim_ratio_counts}) 
ratios.mean.mel <- sapply(results.mean, FUN = function(x){x$mel_ratio_counts}) 
ratios.sim.sim <- sapply(results.sim, FUN = function(x){x$sim_ratio_counts}) 
ratios.mel.mel <- sapply(results.mel, FUN = function(x){x$mel_ratio_counts}) 

ratios.mean.sim.CIs <- sapply(results.mean.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_ratio_counts
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

ratios.mean.mel.CIs <- sapply(results.mean.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_ratio_counts
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

ratios.sim.sim.CIs <- sapply(results.sim.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$sim_ratio_counts
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

ratios.mel.mel.CIs <- sapply(results.mel.bs, FUN = function(x){
  vec <- sapply(x, FUN = function(y){
    y$mel_ratio_counts
  })
  CIs <- quantile(vec, probs = c(0.025, 0.975))
})

#### plotting it ####

# one file, with four panels (mel and sim, A and X), panel labels in the figure?
# colour-blind friendly colour scheme, different symbols. CMYK colour. Tiff, 600dpi,
# 180mm width (but square?)

# install.packages('RColorBrewer')
library(RColorBrewer)

colr <- brewer.pal(n = 3, name = 'Dark2')[1:2]
pnt <- c(1, 2)


pdf('../plots/rates_SI_A_mel_meanbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = rates.AT_2_GC.mean.mel.obs,
     xlab = 'GC content',
     ylab = 'substitution rate',
     main = 'melanogaster (mean bins)',
     col = colr[1],
     ylim = c(0, 0.2),
     axes = F,
     pch = pnt[1])
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2)
axis(2)
box()
points(x = 1:5,
       y = rates.GC_2_AT.mean.mel.obs,
       col = colr[2],
       pch = pnt[2])
legend('topright',
       legend = c(expression(r[W%->%S]), expression(r[S%->%W])),
       col = colr,
       pch = pnt,
       bty = 'n',
       cex = 0.8)
arrows(x0 = 1:5,
       y0 = rates.AT_2_GC.mean.mel.CIs[1,],
       y1 = rates.AT_2_GC.mean.mel.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1, 
       col = colr[1])
arrows(x0 = 1:5,
       y0 = rates.GC_2_AT.mean.mel.CIs[1,],
       y1 = rates.GC_2_AT.mean.mel.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1, 
       col = colr[2])
dev.off()
#######

pdf('../plots/rates_SI_A_sim_meanbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = rates.AT_2_GC.mean.sim.obs,
     xlab = 'GC content',
     ylab = 'substitution rate',
     main = 'simulans (mean bins)',
     col = colr[1],
     ylim = c(0, 0.2),
     axes = F,
     pch = pnt[1])
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.sim.meanbins.SI.A),
     las = 2)
axis(2)
box()
points(x = 1:5,
       y = rates.GC_2_AT.mean.sim.obs,
       col = colr[2],
       pch = pnt[2])
legend('topright',
       legend = c(expression(r[W%->%S]), expression(r[S%->%W])),
       col = colr,
       pch = pnt,
       bty = 'n',
       cex = 0.8)
arrows(x0 = 1:5,
       y0 = rates.AT_2_GC.mean.sim.CIs[1,],
       y1 = rates.AT_2_GC.mean.sim.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1, 
       col = colr[1])
arrows(x0 = 1:5,
       y0 = rates.GC_2_AT.mean.sim.CIs[1,],
       y1 = rates.GC_2_AT.mean.sim.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1, 
       col = colr[2])
dev.off()

###################################################################################

pdf('../plots/rates_SI_A_mel_melbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = rates.AT_2_GC.mel.mel.obs,
     xlab = 'GC content',
     ylab = 'substitution rate',
     main = 'melanogaster (mel bins)',
     col = colr[1],
     ylim = c(0, 0.2),
     axes = F,
     pch = pnt[1])
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.melbins.SI.A),
     las = 2)
axis(2)
# mtext('substitution rate',
#       side = 2,
#       line = 2.2,
#       cex = 1.5)
# mtext('GC content',
#       side = 1,
#       line = 3,
#       cex = 1.5)
box()
points(x = 1:5,
       y = rates.GC_2_AT.mel.mel.obs,
       col = colr[2],
       pch = pnt[2])
legend('topright',
       legend = c(expression(r[W%->%S]), expression(r[S%->%W])),
       col = colr,
       pch = pnt,
       bty = 'n',
       cex = 0.8)
arrows(x0 = 1:5,
       y0 = rates.AT_2_GC.mel.mel.CIs[1,],
       y1 = rates.AT_2_GC.mel.mel.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1, 
       col = colr[1])
arrows(x0 = 1:5,
       y0 = rates.GC_2_AT.mel.mel.CIs[1,],
       y1 = rates.GC_2_AT.mel.mel.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1, 
       col = colr[2])
dev.off()
#######

pdf('../plots/rates_SI_A_sim_simbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = rates.AT_2_GC.sim.sim.obs,
     xlab = 'GC content',
     ylab = 'substitution rate',
     main = 'simulans (sim bins)',
     col = colr[1],
     ylim = c(0, 0.2),
     axes = F,
     pch = pnt[1])
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.sim.simbins.SI.A),
     las = 2)
axis(2)
# mtext('substitution rate',
#       side = 2,
#       line = 2.2,
#       cex = 1.5)
# mtext('GC content',
#       side = 1,
#       line = 3,
#       cex = 1.5)
box()
points(x = 1:5,
       y = rates.GC_2_AT.sim.sim.obs,
       col = colr[2],
       pch = pnt[2])
legend('topright',
       legend = c(expression(r[W%->%S]), expression(r[S%->%W])),
       col = colr,
       pch = pnt,
       bty = 'n',
       cex = 0.8)
arrows(x0 = 1:5,
       y0 = rates.AT_2_GC.sim.sim.CIs[1,],
       y1 = rates.AT_2_GC.sim.sim.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1, 
       col = colr[1])
arrows(x0 = 1:5,
       y0 = rates.GC_2_AT.sim.sim.CIs[1,],
       y1 = rates.GC_2_AT.sim.sim.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1, 
       col = colr[2])
dev.off()


#### RATE END ####

#### RATIO START ####


pdf('../plots/ratios_SI_A_sim_meanbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = ratios.mean.sim,
     xlab = 'GC content',
     ylab = 'substitution count ratio',
     main = 'simulans (mean bins)',
     ylim = c(0, 2),
     axes = F
)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.sim.meanbins.SI.A),
     las = 2)
axis(2)
box()
legend('topright',
       legend = expression(N[W %->% S]/N[S %->% W]),
       pch = '',
       bty = 'n',
       cex = 0.8)
arrows(x0 = 1:5,
       y0 = ratios.mean.sim.CIs[1,],
       y1 = ratios.mean.sim.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1)
dev.off()

pdf('../plots/ratios_SI_A_mel_meanbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = ratios.mean.mel,
     xlab = 'GC content',
     ylab = 'substitution count ratio',
     main = 'melanogaster (mean bins)',
     ylim = c(0, 2),
     axes = F
)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2)
axis(2)
box()
legend('topright',
       legend = expression(N[W %->% S]/N[S %->% W]),
       pch = '',
       bty = 'n',
       cex = 0.8)
arrows(x0 = 1:5,
       y0 = ratios.mean.mel.CIs[1,],
       y1 = ratios.mean.mel.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1)
dev.off()

###############################################

pdf('../plots/ratios_SI_A_sim_simbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = ratios.sim.sim,
     xlab = 'GC content',
     ylab = 'substitution count ratio',
     main = 'simulans (sim bins)',
     ylim = c(0, 2),
     axes = F
)
# mtext('substitution count ratio',
#       side = 2,
#       line = 1.2,
#       cex = 0.75)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.sim.simbins.SI.A),
     las = 2)
axis(2)
box()
legend('topright',
       legend = expression(N[W %->% S]/N[S %->% W]),
       pch = '',
       bty = 'n',
       cex = 0.8)
arrows(x0 = 1:5,
       y0 = ratios.sim.sim.CIs[1,],
       y1 = ratios.sim.sim.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1)
dev.off()

pdf('../plots/ratios_SI_A_mel_melbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = ratios.mel.mel,
     xlab = 'GC content',
     ylab = 'substitution count ratio',
     main = 'melanogaster (mel bins)',
     ylim = c(0, 2),
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
axis(2)
box()
legend('topright',
       legend = expression(N[W %->% S]/N[S %->% W]),
       pch = '',
       bty = 'n',
       cex = 0.8)
arrows(x0 = 1:5,
       y0 = ratios.mel.mel.CIs[1,],
       y1 = ratios.mel.mel.CIs[2,],
       code = 3,
       angle = 90, 
       length = 0.1)
dev.off()
