setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# plotting figure 2

load('../../data/SI_binindices_GCcontent.RData')

load('../../data/DAF_CI_bootstraps.RData')
load('../../data/est-sfs_SFSs.RData')

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


library(RColorBrewer)
col <- brewer.pal(n = 3, name = 'Dark2')
pnt <- c(1, 2, 0)

pdf('figure_2.pdf',
    width = 10,
    height = 10)
par(lwd = 2,
    cex = 8,
    mfrow = c(2,2),
    oma = c(2,2,1,1))
plot(x = 1:5,
     y = AT_to_GC_ALL_DAF.MD.mean,
     xlab = '',
     ylab = '',
     main = '',
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.35),
     axes = F,
     type = 'n')
text(x = 5, 
     y = 0.33,
     labels = c("A"),
     cex = 3)
legend('topleft',
       legend = c(expression(DAF[W %->% S]), expression(DAF[neu]),
                  expression(DAF[S %->% W])),
       col = c(col[c(1,3,2)]),
       pch = c(pnt[c(1,3,2)]),
       bty = 'n',
       cex = 1.8)
points(x = 1:5,
       y = AT_to_GC_ALL_DAF.MD.mean, col = col[1], pch = pnt[1], cex = 2)
points(x = 1:5 - 0.13,
       y = GC_to_AT_ALL_DAF.MD.mean, col = col[2], pch = pnt[2], cex = 2)
points(x = 1:5 + 0.13,
       y = N_to_N_ALL_DAF.MD.mean, col = col[3], pch = pnt[3], cex = 2)
title(main = 'simulans',
      line = 1,
      cex.main = 2)
mtext('derived allele frequency',
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


plot(x = 1:5,
     y = AT_to_GC_ALL_DAF.ZI69.mean,
     xlab = '',
     ylab = '',
     main = '',
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.35),
     axes = F,
     type = 'n'
)
text(x = 1, 
     y = 0.33,
     labels = c("B"),
     cex = 3)
points(x = 1:5,
       y = AT_to_GC_ALL_DAF.ZI69.mean, col = col[1], pch = pnt[1], cex = 2)
points(x = 1:5 - 0.13,
       y = GC_to_AT_ALL_DAF.ZI69.mean, col = col[2], pch = pnt[2], cex = 2)
points(x = 1:5 + 0.13,
       y = N_to_N_ALL_DAF.ZI69.mean, col = col[3], pch = pnt[3], cex = 2)
title(main = 'melanogaster',
      line = 1,
      cex.main = 2)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2, cex.axis = 1.4)
axis(2, cex.axis = 1.4)
box()
mtext('GC content',
      side = 1,
      line = 4,
      cex = 1.3)
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


plot(x = 1:5,
     y = AT_to_GC_ALL_DAF.MD.sim,
     xlab = '',
     ylab = '',
     main = '',
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.35),
     axes = F,
     type = 'n'
)
text(x = 1, 
     y = 0.33,
     labels = c("C"),
     cex = 3)
points(x = 1:5,
       y = AT_to_GC_ALL_DAF.MD.sim, col = col[1], pch = pnt[1], cex = 2)
points(x = 1:5 - 0.13,
       y = GC_to_AT_ALL_DAF.MD.sim, col = col[2], pch = pnt[2], cex = 2)
points(x = 1:5 + 0.13,
       y = N_to_N_ALL_DAF.MD.sim, col = col[3], pch = pnt[3], cex = 2)
title(main = 'simulans',
      line = 1,
      cex.main = 2)
mtext('derived allele frequency',
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


plot(x = 1:5,
     y = AT_to_GC_ALL_DAF.ZI69.mel,
     xlab = '',
     ylab = '',
     main = '',
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.35),
     axes = F,
     type = 'n'
)
text(x = 1, 
     y = 0.33,
     labels = c("D"),
     cex = 3)
points(x = 1:5,
       y = AT_to_GC_ALL_DAF.ZI69.mel, col = col[1], pch = pnt[1], cex = 2)
points(x = 1:5 - 0.13,
       y = GC_to_AT_ALL_DAF.ZI69.mel, col = col[2], pch = pnt[2], cex = 2)
points(x = 1:5 + 0.13,
       y = N_to_N_ALL_DAF.ZI69.mel, col = col[3], pch = pnt[3], cex = 2)
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
