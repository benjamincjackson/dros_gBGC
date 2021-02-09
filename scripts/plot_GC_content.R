setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# plot the relationship between GC contents

load('../data/SI_binindices_GCcontent.RData')

mylm <- lm(GCcontent.Dsim.SI.A ~ GCcontent.Dmel.SI.A)
summary(mylm)

cor(GCcontent.Dsim.SI.A, GCcontent.Dmel.SI.A, method = "kendall")

pdf('../plots/GCcontent.pdf')
plot(GCcontent.Dmel.SI.A, GCcontent.Dsim.SI.A,
     xlab = 'Dmel GC content',
     ylab = 'Dsim GC content')
# abline(mylm,
#        col = 'red')
# text(0.07, 0.75, expression(R^2 == 0.72))
dev.off()

##################################################################

plot(x = GCbinindex.5bins.Dmel.SI.A,
     y = GCbinindex.5bins.Dsim.SI.A)


plot(x = 1:5,
     y = aggregate(GCcontent.Dsim.SI.A,
                   by = list(GCbinindex.5bins.Dsim.SI.A),
                   FUN = max)$x,
     ylim = c(0,1))
points(x = 1:5,
       y = aggregate(GCcontent.Dsim.SI.A,
                     by = list(GCbinindex.5bins.Dsim.SI.A),
                     FUN = min)$x)

plot(x = 1:5,
     y = aggregate(GCcontent.Dmel.SI.A,
                   by = list(GCbinindex.5bins.Dmel.SI.A),
                   FUN = max)$x,
     ylim = c(0,1))
points(x = 1:5,
       y = aggregate(GCcontent.Dmel.SI.A,
                     by = list(GCbinindex.5bins.Dmel.SI.A),
                     FUN = min)$x)

plot(x = 1:5,
     y = aggregate(GCcontent.mean.SI.A,
                   by = list(GCbinindex.5bins.mean.SI.A),
                   FUN = max)$x,
     ylim = c(0,1))
points(x = 1:5,
       y = aggregate(GCcontent.mean.SI.A,
                     by = list(GCbinindex.5bins.mean.SI.A),
                     FUN = min)$x)

plot(x = 1:5,
     y = aggregate(GCcontent.Dsim.SI.A,
                   by = list(GCbinindex.5bins.Dmel.SI.A),
                   FUN = max)$x,
     ylim = c(0,1))
points(x = 1:5,
       y = aggregate(GCcontent.Dsim.SI.A,
                     by = list(GCbinindex.5bins.Dmel.SI.A),
                     FUN = min)$x)

plot(x = 1:5,
     y = aggregate(GCcontent.Dmel.SI.A,
                   by = list(GCbinindex.5bins.Dsim.SI.A),
                   FUN = max)$x,
     ylim = c(0,1))
points(x = 1:5,
       y = aggregate(GCcontent.Dmel.SI.A,
                     by = list(GCbinindex.5bins.Dsim.SI.A),
                     FUN = min)$x)

plot(x = 1:5,
     y = aggregate(GCcontent.Dmel.SI.A,
                   by = list(GCbinindex.5bins.mean.SI.A),
                   FUN = max)$x,
     ylim = c(0,1))
points(x = 1:5,
       y = aggregate(GCcontent.Dmel.SI.A,
                     by = list(GCbinindex.5bins.mean.SI.A),
                     FUN = min)$x)

plot(x = 1:5,
     y = aggregate(GCcontent.Dsim.SI.A,
                   by = list(GCbinindex.5bins.mean.SI.A),
                   FUN = max)$x,
     ylim = c(0,1))
points(x = 1:5,
       y = aggregate(GCcontent.Dsim.SI.A,
                     by = list(GCbinindex.5bins.mean.SI.A),
                     FUN = min)$x)

plot(GCcontent.Dmel.SI.A,
     GCcontent.Dsim.SI.A)

library(seqinr)
library(ape)
library(pegas)
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_alignments.RData')
source('~/Dropbox/Biology_projects/dros_X_A_mut/scripts/R_helper_functions.R')

GCrealsim <- sapply(alignments.intersection, FUN = function(x){
  GC.content(as.DNAbin(editAlignment(x, 'sim')))
})

GCrealmel <- sapply(alignments.intersection, FUN = function(x){
  GC.content(as.DNAbin(editAlignment(x, 'mel')))
})  

plot(GCcontent.Dmel.SI.A,
     GCrealmel) 

plot(GCcontent.Dsim.SI.A,
     GCrealsim)  

summary(lm(GCcontent.Dmel.SI.A ~ GCrealmel))
summary(lm(GCcontent.Dsim.SI.A ~ GCrealsim))

summary(lm(GCrealmel ~ GCrealsim))
summary(lm(GCcontent.Dmel.SI.A ~ GCcontent.Dsim.SI.A))

pdf('~/Dropbox/Biology_projects/dros_X_A_mut/plots/GCcontent_boxplots.pdf',
    onefile = T)
boxplot(GCcontent.Dmel.SI.A ~ GCbinindex.5bins.Dmel.SI.A,
        at = 0:4 * 3 + 1,
        xlim = c(0,16),
        xaxt = "n",
        ylab = 'GC content',
        xlab = 'bin',
        main = 'melanogaster GC content',
        col = 'orange')
boxplot(GCcontent.Dmel.SI.A ~ GCbinindex.5bins.mean.SI.A, add = T,
        at = 0:4 * 3 + 2,
        col = 'forestgreen')
legend('topleft',
       col = c('orange', 'forestgreen'),
       pch = 15,
       legend = c('mel bins', 'mean bins'),
       bty = 'n')

boxplot(GCcontent.Dsim.SI.A ~ GCbinindex.5bins.Dsim.SI.A,
        at = 0:4 * 3 + 1,
        xlim = c(0,16),
        xaxt = "n",
        ylab = 'GC content',
        xlab = 'bin',
        main = 'simulans GC content',
        col = 'orange')
boxplot(GCcontent.Dsim.SI.A ~ GCbinindex.5bins.mean.SI.A, add = T,
        at = 0:4 * 3 + 2,
        col = 'forestgreen')
legend('topleft',
       col = c('orange', 'forestgreen'),
       pch = 15,
       legend = c('sim bins', 'mean bins'),
       bty = 'n')


boxplot(GCcontent.Dsim.SI.A ~ GCbinindex.5bins.Dsim.SI.A,
        at = 0:4 * 3 + 1,
        xlim = c(0,16),
        xaxt = "n",
        ylab = 'GC content',
        xlab = 'bin',
        main = 'simulans bins',
        col = 'orange')
boxplot(GCcontent.Dmel.SI.A ~ GCbinindex.5bins.Dsim.SI.A, add = T,
        at = 0:4 * 3 + 2,
        col = 'forestgreen')
legend('topleft',
       col = c('orange', 'forestgreen'),
       pch = 15,
       legend = c('sim GC content', 'mel GC content'),
       bty = 'n')

boxplot(GCcontent.Dmel.SI.A ~ GCbinindex.5bins.Dmel.SI.A,
        at = 0:4 * 3 + 1,
        xlim = c(0,16),
        xaxt = "n",
        ylab = 'GC content',
        xlab = 'bin',
        main = 'melanogaster bins',
        col = 'orange')
boxplot(GCcontent.Dsim.SI.A ~ GCbinindex.5bins.Dmel.SI.A, add = T,
        at = 0:4 * 3 + 2,
        col = 'forestgreen')
legend('topleft',
       col = c('orange', 'forestgreen'),
       pch = 15,
       legend = c('mel GC content', 'sim GC content'),
       bty = 'n')
dev.off()




