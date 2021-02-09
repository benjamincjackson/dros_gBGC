load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_alignments.RData')
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/recombination_Dmel.RData')
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_binindices_GCcontent.RData')
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_binindices_GCcontent_DIFFERENCES.RData')
source('~/Dropbox/Biology_projects/dros_X_A_mut/scripts/R_helper_functions.R')

################################################################################
# first, to sanity-check recombination, does it correlate with nucleotide diversity?
library(seqinr)
library(pegas)

piZI <- sapply(alignments.intersection.ZI69.SI.A, FUN = function(x){
  nuc.div(as.DNAbin(editAlignment(x, 'ZI')))
})


mydf = data.frame('pi' = piZI, 'r_raw' = recomb.comeron.dmel.SI.A, 'r_smooth' = recomb.comeronSmoothed.dmel.SI.A)

pdf('~/Desktop/ZI_SI_A_vs_recomb_comeron_raw.pdf')
mylm <- lm(pi ~ r_raw, mydf)
mylm_squared <- lm(pi ~ r_raw + I(r_raw^2), mydf)
plot(pi ~ r_raw, mydf,
     main = 'Zambian autosomal SIs',
     ylab = expression(pi),
     xlab = 'recombination rate (raw)')
curve(predict(mylm, newdata=data.frame(r_raw=x)), add = TRUE, col = 'red')
curve(predict(mylm_squared, newdata=data.frame(r_raw=x)), add = TRUE, col = 'blue')
legend('topright',
       legend = c('linear', 'squared'),
       pch = '-',
       col = c('red', 'blue'),
       bty = 'n')
dev.off()

pdf('~/Desktop/ZI_SI_A_vs_recomb_comeron_smooth.pdf')
mylm_smoothed <- lm(pi ~ r_smooth, mydf)
mylm_smoothed_squared <- lm(pi ~ r_smooth + I(r_smooth^2), mydf)
plot(pi ~ r_smooth, mydf,
     main = 'Zambian autosomal SIs',
     ylab = expression(pi),
     xlab = 'recombination rate (smoothed)')
curve(predict(mylm_smoothed, newdata=data.frame(r_smooth=x)), add = TRUE, col = 'red')
curve(predict(mylm_smoothed_squared, newdata=data.frame(r_smooth=x)), add = TRUE, col = 'blue')
legend('topleft',
       legend = c('linear', 'squared'),
       pch = '-',
       col = c('red', 'blue'),
       bty = 'n')
dev.off()

mylm <- lm(piZI ~ RRC$RRCmid)
plot(y = piZI,
     x = RRC$RRCmid)
abline(mylm, col = 'red')

################################################################################
# and then here are some plots of recombination rate against GC content: 
boxplot(
recomb.comeron.dmel.SI.A ~ GCbinindex.5bins.diff.SI.A,
ylim = c(0,5)
)

boxplot(
  recomb.comeron.dmel.SI.A ~ GCbinindex.5bins.mean.SI.A,
  ylim = c(0,5)
)

boxplot(
  recomb.comeron.dmel.SI.A ~ GCbinindex.5bins.Dmel.SI.A,
  ylim = c(0,5)
)

boxplot(
  recomb.comeronSmoothed.dmel.SI.A ~ GCbinindex.5bins.diff.SI.A
)

boxplot(
  recomb.comeronSmoothed.dmel.SI.A ~ GCbinindex.5bins.mean.SI.A
)

boxplot(
  recomb.comeronSmoothed.dmel.SI.A ~ GCbinindex.5bins.Dmel.SI.A
)


boxplot(
  RRC$RRCmid ~ GCbinindex.5bins.diff.SI.A
)

boxplot(
  RRC$RRCmid ~ GCbinindex.5bins.mean.SI.A
)

boxplot(
  RRC$RRCmid ~ GCbinindex.5bins.Dmel.SI.A
)


mylm <- lm(GCcontent.Dmel.SI.A ~ recomb.comeron.dmel.SI.A)
plot(y = GCcontent.Dmel.SI.A,
     x = recomb.comeron.dmel.SI.A)
abline(mylm, col = 'red')

mylm <- lm(GCcontent.diff ~ recomb.comeron.dmel.SI.A)
plot(y = GCcontent.diff,
     x = recomb.comeron.dmel.SI.A)
abline(mylm, col = 'red')

mylm <- lm(GCcontent.mean.SI.A ~ recomb.comeron.dmel.SI.A)
plot(y = GCcontent.mean.SI.A,
     x = recomb.comeron.dmel.SI.A)
abline(mylm, col = 'red')


mylm <- lm(GCcontent.Dmel.SI.A ~ recomb.comeronSmoothed.dmel.SI.A)
plot(y = GCcontent.Dmel.SI.A,
     x = recomb.comeronSmoothed.dmel.SI.A)
abline(mylm, col = 'red')

mylm <- lm(GCcontent.diff ~ recomb.comeronSmoothed.dmel.SI.A)
plot(y = GCcontent.diff,
     x = recomb.comeronSmoothed.dmel.SI.A)
abline(mylm, col = 'red')

mylm <- lm(GCcontent.mean.SI.A ~ recomb.comeronSmoothed.dmel.SI.A)
plot(y = GCcontent.mean.SI.A,
     x = recomb.comeronSmoothed.dmel.SI.A)
abline(mylm, col = 'red')


mylm <- lm(GCcontent.Dmel.SI.A ~ RRC$RRCmid)
plot(y = GCcontent.Dmel.SI.A,
     x = RRC$RRCmid)
abline(mylm, col = 'red')

mylm <- lm(GCcontent.diff ~ RRC$RRCmid)
plot(y = GCcontent.diff,
     x = RRC$RRCmid)
abline(mylm, col = 'red')

mylm <- lm(GCcontent.mean.SI.A ~ RRC$RRCmid)
plot(y = GCcontent.mean.SI.A,
     x = RRC$RRCmid)
abline(mylm, col = 'red')









#