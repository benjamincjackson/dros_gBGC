setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# plot the alignment scores for bins of introns

library(simpleboot)
library(boot)

scoresTable <- read.delim('../data/scores.txt', header = F)

load('../data/SI_binindices_GCcontent.RData')
load('../data/SI_alignments.RData')

LUT <- scoresTable[,2]
names(LUT) <- scoresTable[,1]

scores <- sapply(files.intersection, FUN = function(x){
  LUT[strsplit(x, '&')[[1]][3]]
})

names(scores) <- NULL

scores.mean <- vector('list', 5)
for(i in 1:5){
  scores.mean[[i]] <- scores[GCbinindex.5bins.mean.SI.A == i]
}

scores.mel <- vector('list', 5)
for(i in 1:5){
  scores.mel[[i]] <- scores[GCbinindex.5bins.Dmel.SI.A == i]
}

scores.sim <- vector('list', 5)
for(i in 1:5){
  scores.sim[[i]] <- scores[GCbinindex.5bins.Dsim.SI.A == i]
}


CIs.mean.lower <- vector('numeric', 5)
CIs.mean.higher <- vector('numeric', 5)
for(i in 1:5){
  CIs.mean.lower[i] <- boot.ci(one.boot(data = scores.mean[[i]],
                                    FUN = median,
                                    R = 1000),
                           type = 'basic')$basic[4]
  CIs.mean.higher[i] <- boot.ci(one.boot(data = scores.mean[[i]],
                                        FUN = median,
                                        R = 1000),
                               type = 'basic')$basic[5]
}

pdf('../plots/alignment_mean.pdf')
par(lwd = 2,
    cex = 1.6)
plot(sapply(scores.mean, median),
     ylim = c(2500000, 3500000),
     main = 'mean bins',
     xlab = 'GC content',
     ylab = 'alignment score',
     axes = F)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2)
axis(2)
box()
arrows(x0 = 1:5,
       x1 = 1:5,
       y0 = CIs.mean.lower,
       y1 = CIs.mean.higher,
       angle = 90,
       code = 3,
       length = 0.1)
dev.off()


CIs.mel.lower <- vector('numeric', 5)
CIs.mel.higher <- vector('numeric', 5)
for(i in 1:5){
  CIs.mel.lower[i] <- boot.ci(one.boot(data = scores.mel[[i]],
                                        FUN = median,
                                        R = 1000),
                               type = 'basic')$basic[4]
  CIs.mel.higher[i] <- boot.ci(one.boot(data = scores.mel[[i]],
                                         FUN = median,
                                         R = 1000),
                                type = 'basic')$basic[5]
}

pdf('../plots/alignment_mel.pdf')
par(lwd = 2,
    cex = 1.6)
plot(sapply(scores.mel, median),
     ylim = c(2500000, 3500000),
     main = 'mel bins',
     xlab = 'GC content',
     ylab = 'alignment score',
     axes = F)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.mel.melbins.SI.A),
     las = 2)
axis(2)
box()
arrows(x0 = 1:5,
       x1 = 1:5,
       y0 = CIs.mel.lower,
       y1 = CIs.mel.higher,
       angle = 90,
       code = 3,
       length = 0.1)
dev.off()


CIs.sim.lower <- vector('numeric', 5)
CIs.sim.higher <- vector('numeric', 5)
for(i in 1:5){
  CIs.sim.lower[i] <- boot.ci(one.boot(data = scores.sim[[i]],
                                       FUN = median,
                                       R = 1000),
                              type = 'basic')$basic[4]
  CIs.sim.higher[i] <- boot.ci(one.boot(data = scores.sim[[i]],
                                        FUN = median,
                                        R = 1000),
                               type = 'basic')$basic[5]
}

pdf('../plots/alignment_sim.pdf')
par(lwd = 2,
    cex = 1.6)
plot(sapply(scores.sim, median),
     ylim = c(2500000, 3500000),
     main = 'sim bins',
     xlab = 'GC content',
     ylab = 'alignment score',
     axes = F)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.sim.simbins.SI.A),
     las = 2)
axis(2)
box()
arrows(x0 = 1:5,
       x1 = 1:5,
       y0 = CIs.sim.lower,
       y1 = CIs.sim.higher,
       angle = 90,
       code = 3,
       length = 0.1)
dev.off()


