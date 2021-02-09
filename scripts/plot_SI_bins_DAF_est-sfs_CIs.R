setwd(dirname(rstudioapi::getSourceEditorContext()$path))

load('../data/SI_binindices_GCcontent.RData')

load('../data/DAF_CI_bootstraps.RData')
load('../data/est-sfs_SFSs.RData')

N_to_N_ALL_DAF.MD.mean <- sapply(N_to_N_ALL_SFSes.MD.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
AT_to_GC_ALL_DAF.MD.mean <- sapply(AT_to_GC_ALL_SFSes.MD.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
GC_to_AT_ALL_DAF.MD.mean <- sapply(GC_to_AT_ALL_SFSes.MD.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })

N_to_N_ALL_DAF.MD.sim <- sapply(N_to_N_ALL_SFSes.MD.SI.A.sim.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
AT_to_GC_ALL_DAF.MD.sim <- sapply(AT_to_GC_ALL_SFSes.MD.SI.A.sim.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
GC_to_AT_ALL_DAF.MD.sim <- sapply(GC_to_AT_ALL_SFSes.MD.SI.A.sim.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })

N_to_N_ALL_DAF.ZI69.mean <- sapply(N_to_N_ALL_SFSes.ZI69.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:69] * 1:68/69) /  sum(x[2:69]) })
AT_to_GC_ALL_DAF.ZI69.mean <- sapply(AT_to_GC_ALL_SFSes.ZI69.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:69] * 1:68/69) /  sum(x[2:69]) })
GC_to_AT_ALL_DAF.ZI69.mean <- sapply(GC_to_AT_ALL_SFSes.ZI69.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:69] * 1:68/69) /  sum(x[2:69]) })

N_to_N_ALL_DAF.ZI69.mel <- sapply(N_to_N_ALL_SFSes.ZI69.SI.A.mel.5GCbins, FUN = function(x){ sum(x[2:69] * 1:68/69) /  sum(x[2:69]) })
AT_to_GC_ALL_DAF.ZI69.mel <- sapply(AT_to_GC_ALL_SFSes.ZI69.SI.A.mel.5GCbins, FUN = function(x){ sum(x[2:69] * 1:68/69) /  sum(x[2:69]) })
GC_to_AT_ALL_DAF.ZI69.mel <- sapply(GC_to_AT_ALL_SFSes.ZI69.SI.A.mel.5GCbins, FUN = function(x){ sum(x[2:69] * 1:68/69) /  sum(x[2:69]) })

N_to_N_ALL_DAF.ZI21.mean <- sapply(N_to_N_ALL_SFSes.ZI21.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
AT_to_GC_ALL_DAF.ZI21.mean <- sapply(AT_to_GC_ALL_SFSes.ZI21.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
GC_to_AT_ALL_DAF.ZI21.mean <- sapply(GC_to_AT_ALL_SFSes.ZI21.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })

N_to_N_ALL_DAF.ZI21.mel <- sapply(N_to_N_ALL_SFSes.ZI21.SI.A.mel.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
AT_to_GC_ALL_DAF.ZI21.mel <- sapply(AT_to_GC_ALL_SFSes.ZI21.SI.A.mel.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
GC_to_AT_ALL_DAF.ZI21.mel <- sapply(GC_to_AT_ALL_SFSes.ZI21.SI.A.mel.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })


N_to_N_ALL_DAF.RG21.mean <- sapply(N_to_N_ALL_SFSes.RG21.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
AT_to_GC_ALL_DAF.RG21.mean <- sapply(AT_to_GC_ALL_SFSes.RG21.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
GC_to_AT_ALL_DAF.RG21.mean <- sapply(GC_to_AT_ALL_SFSes.RG21.SI.A.mean.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })

N_to_N_ALL_DAF.RG21.mel <- sapply(N_to_N_ALL_SFSes.RG21.SI.A.mel.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
AT_to_GC_ALL_DAF.RG21.mel <- sapply(AT_to_GC_ALL_SFSes.RG21.SI.A.mel.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })
GC_to_AT_ALL_DAF.RG21.mel <- sapply(GC_to_AT_ALL_SFSes.RG21.SI.A.mel.5GCbins, FUN = function(x){ sum(x[2:21] * 1:20/21) /  sum(x[2:21]) })


library(RColorBrewer)
col <- brewer.pal(n = 3, name = 'Dark2')
pnt <- c(1, 2, 0)


pdf('../plots/DAF_SI_A_MD_est-sfs_meanbins_CIs.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = AT_to_GC_ALL_DAF.MD.mean,
     main = 'MD (mean bins)',
     xlab = 'GC content',
     ylab = 'derived allele frequency',
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.35),
     axes = F,
     col = col[1],
     pch = pnt[1])
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.sim.meanbins.SI.A),
     las = 2)
axis(2)
box()
points(x = 1:5 - 0.13,
       y = GC_to_AT_ALL_DAF.MD.mean, col = col[2], pch = pnt[2])
points(x = 1:5 + 0.13,
       y = N_to_N_ALL_DAF.MD.mean, col = col[3], pch = pnt[3])
legend('topleft',
       legend = c(expression(DAF[W %->% S]), expression(DAF[neu]),
                  expression(DAF[S %->% W])),
       col = c(col[c(1,3,2)]),
       pch = c(pnt[c(1,3,2)]),
       bty = 'n')
for(i in 1:5){
  temp <- lapply(get(paste0('DAFCIs.MD21.A.mean.GCbin0', i)),
                 quantile, c(0.025, 0.975))
  arrows(x0 = i + 0.13,
         x1 = i + 0.13,
         y0 = temp[[1]][1],
         y1 = temp[[1]][2],
         col = col[3],
         angle = 90,
         code = 3,
         length = 0.1)
  
  arrows(x0 = i,
         x1 = i,
         y0 = temp[[2]][1],
         y1 = temp[[2]][2],
         col = col[1],
         angle = 90,
         code = 3,
         length = 0.1)
  
  arrows(x0 = i - 0.13,
         x1 = i - 0.13,
         y0 = temp[[3]][1],
         y1 = temp[[3]][2],
         col = col[2],
         angle = 90,
         code = 3,
         length = 0.1)
}
dev.off()

pdf('../plots/DAF_SI_A_MD_est-sfs_simbins_CIs.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = AT_to_GC_ALL_DAF.MD.sim,
     main = 'MD (sim bins)',
     xlab = 'GC content',
     ylab = 'derived allele frequency',
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.35),
     axes = F,
     col = col[1],
     pch = pnt[1])
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.sim.simbins.SI.A),
     las = 2)
axis(2)
box()
points(x = 1:5 - 0.13,
       y = GC_to_AT_ALL_DAF.MD.sim, col = col[2], pch = pnt[2])
points(x = 1:5 + 0.13,
       y = N_to_N_ALL_DAF.MD.sim, col = col[3], pch = pnt[3])
legend('topleft',
       legend = c(expression(DAF[W %->% S]), expression(DAF[neu]),
                  expression(DAF[S %->% W])),
       col = c(col[c(1,3,2)]),
       pch = c(pnt[c(1,3,2)]),
       bty = 'n')
for(i in 1:5){
  temp <- lapply(get(paste0('DAFCIs.MD21.A.sim.GCbin0', i)),
                 quantile, c(0.025, 0.975))
  arrows(x0 = i + 0.13,
         x1 = i + 0.13,
         y0 = temp[[1]][1],
         y1 = temp[[1]][2],
         col = col[3],
         angle = 90,
         code = 3,
         length = 0.1)
  
  arrows(x0 = i,
         x1 = i,
         y0 = temp[[2]][1],
         y1 = temp[[2]][2],
         col = col[1],
         angle = 90,
         code = 3,
         length = 0.1)
  
  arrows(x0 = i - 0.13,
         x1 = i - 0.13,
         y0 = temp[[3]][1],
         y1 = temp[[3]][2],
         col = col[2],
         angle = 90,
         code = 3,
         length = 0.1)
}
dev.off()


pdf('../plots/DAF_SI_A_ZI69_est-sfs_meanbins_CIs.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = AT_to_GC_ALL_DAF.ZI69.mean,
     main = 'ZI69 (mean bins)',
     xlab = 'GC content',
     ylab = 'derived allele frequency',
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.35),
     axes = F,
     col = col[1],
     pch = pnt[1])
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2)
axis(2)
box()
points(x = 1:5 - 0.13,
       y = GC_to_AT_ALL_DAF.ZI69.mean, col = col[2], pch = pnt[2])
points(x = 1:5 + 0.13,
       y = N_to_N_ALL_DAF.ZI69.mean, col = col[3], pch = pnt[3])
legend('topleft',
       legend = c(expression(DAF[W %->% S]), expression(DAF[neu]),
                  expression(DAF[S %->% W])),
       col = c(col[c(1,3,2)]),
       pch = c(pnt[c(1,3,2)]),
       bty = 'n')
for(i in 1:5){
  temp <- lapply(get(paste0('DAFCIs.ZI69.A.mean.GCbin0', i)),
                 quantile, c(0.025, 0.975))
  arrows(x0 = i + 0.13,
         x1 = i + 0.13,
         y0 = temp[[1]][1],
         y1 = temp[[1]][2],
         col = col[3],
         angle = 90,
         code = 3,
         length = 0.1)
  
  arrows(x0 = i,
         x1 = i,
         y0 = temp[[2]][1],
         y1 = temp[[2]][2],
         col = col[1],
         angle = 90,
         code = 3,
         length = 0.1)
  
  arrows(x0 = i - 0.13,
         x1 = i - 0.13,
         y0 = temp[[3]][1],
         y1 = temp[[3]][2],
         col = col[2],
         angle = 90,
         code = 3,
         length = 0.1)
}
dev.off()


pdf('../plots/DAF_SI_A_ZI69_est-sfs_melbins_CIs.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = AT_to_GC_ALL_DAF.ZI69.mel,
     main = 'ZI69 (mel bins)',
     xlab = 'GC content',
     ylab = 'derived allele frequency',
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.35),
     axes = F,
     col = col[1],
     pch = pnt[1])
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.melbins.SI.A),
     las = 2)
axis(2)
box()
points(x = 1:5 - 0.13,
       y = GC_to_AT_ALL_DAF.ZI69.mel, col = col[2], pch = pnt[2])
points(x = 1:5 + 0.13,
       y = N_to_N_ALL_DAF.ZI69.mel, col = col[3], pch = pnt[3])
legend('topleft',
       legend = c(expression(DAF[W %->% S]), expression(DAF[neu]),
                  expression(DAF[S %->% W])),
       col = c(col[c(1,3,2)]),
       pch = c(pnt[c(1,3,2)]),
       bty = 'n')
for(i in 1:5){
  temp <- lapply(get(paste0('DAFCIs.ZI69.A.mel.GCbin0', i)),
                 quantile, c(0.025, 0.975))
  arrows(x0 = i + 0.13,
         x1 = i + 0.13,
         y0 = temp[[1]][1],
         y1 = temp[[1]][2],
         col = col[3],
         angle = 90,
         code = 3,
         length = 0.1)
  
  arrows(x0 = i,
         x1 = i,
         y0 = temp[[2]][1],
         y1 = temp[[2]][2],
         col = col[1],
         angle = 90,
         code = 3,
         length = 0.1)
  
  arrows(x0 = i - 0.13,
         x1 = i - 0.13,
         y0 = temp[[3]][1],
         y1 = temp[[3]][2],
         col = col[2],
         angle = 90,
         code = 3,
         length = 0.1)
}
dev.off()
