setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# get pi + 95% CIs per GC content bin

# make an SFS per intron and bootstrap inton SFSes within bins to get CIs

# these are the data:
load('../data/SI_binindices_GCcontent.RData')
load('../data/SI_binindices_GCcontent_DIFFERENCES.RData')
load('../data/SI_alignments.RData')

# these are functions for dealing with alignments:
library('seqinr')
source('../scripts/R_helper_functions.R')
source('../scripts/sumstats.R')

# get rid of the reference sequences
SFSes.MD21 <- lapply(alignments.intersection.MD21.SI.A, FUN = function(x){
  SFS <- get_SFS_from_alignment(x, "MD")
})

SFSes.ZI69 <- lapply(alignments.intersection.ZI69.SI.A, FUN = function(x){
  SFS <- get_SFS_from_alignment(x, "ZI")
})

# aggregate SFSes by GC content bins
SFSes.MD21.aggregate.mean <- split(SFSes.MD21, GCbinindex.5bins.mean.SI.A)
SFSes.MD21.aggregate.sim <- split(SFSes.MD21, GCbinindex.5bins.Dsim.SI.A)
SFSes.MD21.aggregate.diff <- split(SFSes.MD21, GCbinindex.5bins.diff.SI.A)
SFSes.ZI69.aggregate.mean <- split(SFSes.ZI69, GCbinindex.5bins.mean.SI.A)
SFSes.ZI69.aggregate.mel <- split(SFSes.ZI69, GCbinindex.5bins.Dmel.SI.A)
SFSes.ZI69.aggregate.diff <- split(SFSes.ZI69, GCbinindex.5bins.diff.SI.A)

# get point estimates of pi
MD21.mean.pi.obs <- sapply(SFSes.MD21.aggregate.mean, FUN = function(x) {
  pi_from_uSFS(sum_SFS(x))
})

MD21.sim.pi.obs <- sapply(SFSes.MD21.aggregate.sim, FUN = function(x) {
  pi_from_uSFS(sum_SFS(x))
})

MD21.diff.pi.obs <- sapply(SFSes.MD21.aggregate.diff, FUN = function(x) {
  pi_from_uSFS(sum_SFS(x))
})

ZI69.mean.pi.obs <- sapply(SFSes.ZI69.aggregate.mean, FUN = function(x) {
  pi_from_uSFS(sum_SFS(x))
})

ZI69.mel.pi.obs <- sapply(SFSes.ZI69.aggregate.mel, FUN = function(x) {
  pi_from_uSFS(sum_SFS(x))
})

ZI69.diff.pi.obs <- sapply(SFSes.ZI69.aggregate.diff, FUN = function(x) {
  pi_from_uSFS(sum_SFS(x))
})

# bootstrap the data
MD21.mean.pi.bs <- lapply(SFSes.MD21.aggregate.mean, FUN = function(x) {
  bootstrap_sumstats(x)$pi
})

MD21.sim.pi.bs <- lapply(SFSes.MD21.aggregate.sim, FUN = function(x) {
  bootstrap_sumstats(x)$pi
})

MD21.diff.pi.bs <- lapply(SFSes.MD21.aggregate.diff, FUN = function(x) {
  bootstrap_sumstats(x)$pi
})

ZI69.mean.pi.bs <- lapply(SFSes.ZI69.aggregate.mean, FUN = function(x) {
  bootstrap_sumstats(x)$pi
})

ZI69.mel.pi.bs <- lapply(SFSes.ZI69.aggregate.mel, FUN = function(x) {
  bootstrap_sumstats(x)$pi
})

ZI69.diff.pi.bs <- lapply(SFSes.ZI69.aggregate.diff, FUN = function(x) {
  bootstrap_sumstats(x)$pi
})

# and get the CIs
MD21.mean.pi.CIs <- lapply(MD21.mean.pi.bs, FUN = quantile, c(0.025, 0.975))
MD21.sim.pi.CIs <- lapply(MD21.sim.pi.bs, FUN = quantile, c(0.025, 0.975))
MD21.diff.pi.CIs <- lapply(MD21.diff.pi.bs, FUN = quantile, c(0.025, 0.975))
ZI69.mean.pi.CIs <- lapply(ZI69.mean.pi.bs, FUN = quantile, c(0.025, 0.975))
ZI69.mel.pi.CIs <- lapply(ZI69.mel.pi.bs, FUN = quantile, c(0.025, 0.975))
ZI69.diff.pi.CIs <- lapply(ZI69.diff.pi.bs, FUN = quantile, c(0.025, 0.975))

# need to get the diff bins GC contents because these arent anywhere else:
GCcontent.concat.5bins.sim.diffbins.SI.A <- lapply(get_concat_seqs(alignments.intersection.MD21.SI.A,
                                                  GCbinindex.5bins.diff.SI.A),
                      FUN = function(al) {
                        ref <- editAlignment(al, "sim")
                        GC <- GC(s2c(ref$seq[[1]]))
                      })

GCcontent.concat.5bins.mel.diffbins.SI.A <- lapply(get_concat_seqs(alignments.intersection.ZI69.SI.A,
                                                                   GCbinindex.5bins.diff.SI.A),
                                                   FUN = function(al) {
                                                     ref <- editAlignment(al, "mel")
                                                     GC <- GC(s2c(ref$seq[[1]]))
                                                   })


# plot things
pdf('../plots/pi_SI_A_MD_meanbins_CIs.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = MD21.mean.pi.obs,
     main = 'MD (mean bins)',
     xlab = 'GC content',
     ylab = expression(pi),
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.06),
     axes = F)
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.sim.meanbins.SI.A),
     las = 2)
axis(2)
box()
for(i in 1:5){
  arrows(x0 = i,
         x1 = i,
         y0 = MD21.mean.pi.CIs[[i]][1],
         y1 = MD21.mean.pi.CIs[[i]][2],
         angle = 90,
         code = 3,
         length = 0.1)
}
points(x = 1:5,
       y = MD21.mean.pi.obs)
dev.off()


pdf('../plots/pi_SI_A_MD_simbins_CIs.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = MD21.sim.pi.obs,
     main = 'MD (sim bins)',
     xlab = 'GC content',
     ylab = expression(pi),
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.06),
     axes = F)
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.sim.simbins.SI.A),
     las = 2)
axis(2)
box()
for(i in 1:5){
  arrows(x0 = i,
         x1 = i,
         y0 = MD21.sim.pi.CIs[[i]][1],
         y1 = MD21.sim.pi.CIs[[i]][2],
         angle = 90,
         code = 3,
         length = 0.1)
}
points(x = 1:5,
       y = MD21.sim.pi.obs)
dev.off()


pdf('../plots/pi_SI_A_MD_diffbins_CIs.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = MD21.diff.pi.obs,
     main = 'MD (diff bins)',
     xlab = 'GC content',
     ylab = expression(pi),
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.06),
     axes = F)
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.sim.diffbins.SI.A),
     las = 2)
axis(2)
box()
for(i in 1:5){
  arrows(x0 = i,
         x1 = i,
         y0 = MD21.diff.pi.CIs[[i]][1],
         y1 = MD21.diff.pi.CIs[[i]][2],
         angle = 90,
         code = 3,
         length = 0.1)
}
points(x = 1:5,
       y = MD21.diff.pi.obs)
dev.off()


pdf('../plots/pi_SI_A_ZI_meanbins_CIs.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = ZI69.mean.pi.obs,
     main = 'ZI (mean bins)',
     xlab = 'GC content',
     ylab = expression(pi),
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.06),
     axes = F)
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2)
axis(2)
box()
for(i in 1:5){
  arrows(x0 = i,
         x1 = i,
         y0 = ZI69.mean.pi.CIs[[i]][1],
         y1 = ZI69.mean.pi.CIs[[i]][2],
         angle = 90,
         code = 3,
         length = 0.1)
}
points(x = 1:5,
       y = ZI69.mean.pi.obs)
dev.off()


pdf('../plots/pi_SI_A_ZI_melbins_CIs.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = ZI69.mel.pi.obs,
     main = 'ZI (mel bins)',
     xlab = 'GC content',
     ylab = expression(pi),
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.06),
     axes = F)
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.melbins.SI.A),
     las = 2)
axis(2)
box()
for(i in 1:5){
  arrows(x0 = i,
         x1 = i,
         y0 = ZI69.mel.pi.CIs[[i]][1],
         y1 = ZI69.mel.pi.CIs[[i]][2],
         angle = 90,
         code = 3,
         length = 0.1)
}
points(x = 1:5,
       y = ZI69.mel.pi.obs)
dev.off()


pdf('../plots/pi_SI_A_ZI_diffbins_CIs.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = ZI69.diff.pi.obs,
     main = 'ZI (diff bins)',
     xlab = 'GC content',
     ylab = expression(pi),
     xlim = c(0.8, 5.2),
     ylim = c(0, 0.06),
     axes = F)
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.diffbins.SI.A),
     las = 2)
axis(2)
box()
for(i in 1:5){
  arrows(x0 = i,
         x1 = i,
         y0 = ZI69.diff.pi.CIs[[i]][1],
         y1 = ZI69.diff.pi.CIs[[i]][2],
         angle = 90,
         code = 3,
         length = 0.1)
}
points(x = 1:5,
       y = ZI69.diff.pi.obs)
dev.off()





#