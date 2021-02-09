setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# plotting figure 3

# load the results
load('../../data/glemin_results.RData')
# load the gc content and bin indices
load('../../data/SI_binindices_GCcontent.RData')

pdf('figure_3.pdf',
    width = 10,
    height = 10)
par(lwd = 2,
    cex = 8,
    mfrow = c(2,2),
    oma = c(2,2,1,1))
plot(x = 1:5,
     y = sapply(results_MD.SI.A.5GCbins.mean_error, FUN = function(x){
       x$with_error$model1$B
     }),
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(-1, 3),
     axes = F,
     type = 'n')
text(x = 1.1, 
     y = 2.8,
     labels = c("A"),
     cex = 3)
points(x = 1:5,
       y = sapply(results_MD.SI.A.5GCbins.mean_error, FUN = function(x){
         x$with_error$model1$B
       }),
       pch = c(1, 16, 16, 16, 16),
       cex = 2)
points(x = 1:5,
       y = sapply(results_MD.SI.A.5GCbins.mean_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(17, 17, 17, 17, 17),
       cex = 2)
title(main = 'simulans',
      line = 1,
      cex.main = 2)
mtext(expression(gamma),
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

plot(x = 1:5,
     y = sapply(results_ZI69.SI.A.5GCbins.mean_error, FUN = function(x){
       x$with_error$model1$B
     }),
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(-1, 3),
     axes = F,
     type = 'n'
)
legend('topright',
       bty = 'n',
       legend = c('model M1*, gamma != 0',
                  'model M1, gamma != 0',
                  'model M1*, gamma = 0',
                  'model M1, gamma = 0'),
       pch = c(16, 17, 1, 2),
       cex = 1.2)
text(x = 1.1, 
     y = 2.8,
     labels = c("B"),
     cex = 3)
points(x = 1:5,
       y = sapply(results_ZI69.SI.A.5GCbins.mean_error, FUN = function(x){
         x$with_error$model1$B
       }),
       pch = c(1, 1, 16, 16, 16),
       cex = 2)
points(x = 1:5,
       y = sapply(results_ZI69.SI.A.5GCbins.mean_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(2, 2, 17, 17, 17),
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
mtext('GC content',
      side = 1,
      line = 4,
      cex = 1.3)


plot(x = 1:5,
     y = sapply(results_MD.SI.A.5GCbins.sim_error, FUN = function(x){
       x$with_error$model1$B
     }),
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(-1, 3),
     axes = F,
     type = 'n')
text(x = 1.1, 
     y = 2.8,
     labels = c("C"),
     cex = 3)
points(x = 1:5,
       y = sapply(results_MD.SI.A.5GCbins.sim_error, FUN = function(x){
         x$with_error$model1$B
       }),
       pch = c(1, 16, 16, 16, 16),
       cex = 2)
points(x = 1:5,
       y = sapply(results_MD.SI.A.5GCbins.sim_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(2, 17, 17, 17, 17),
       cex = 2)
title(main = 'simulans',
      line = 1,
      cex.main = 2)
mtext(expression(gamma),
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


plot(x = 1:5,
     y = sapply(results_ZI69.SI.A.5GCbins.mean_error, FUN = function(x){
       x$with_error$model1$B
     }),
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(-1, 3),
     axes = F,
     type = 'n'
)
text(x = 1.1, 
     y = 2.8,
     labels = c("D"),
     cex = 3)
points(x = 1:5,
       y = sapply(results_ZI69.SI.A.5GCbins.mel_error, FUN = function(x){
         x$with_error$model1$B
       }),
       pch = c(16, 1, 16, 16, 16),
       cex = 2)
points(x = 1:5,
       y = sapply(results_ZI69.SI.A.5GCbins.mel_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(17, 2, 17, 17, 17),
       cex = 2)
title(main = 'melanogaster',
      line = 1,
      cex.main = 2)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      GCcontent.concat.5bins.mel.melbins.SI.A),
     las = 2, cex.axis = 1.4)
axis(2, cex.axis = 1.4)
box()
mtext('GC content',
      side = 1,
      line = 4,
      cex = 1.3)
dev.off()
