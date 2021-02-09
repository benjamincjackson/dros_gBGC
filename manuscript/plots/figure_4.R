setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# plotting figure 4

# load the results
load('../../data/glemin_results.RData')
# load the gc content and bin indices
load('../../data/SI_binindices_GCcontent.RData')


pdf('figure_4.pdf',
    width = 10,
    height = 10)
par(lwd = 2,
    cex = 8,
    mfrow = c(2,2),
    oma = c(2,2,1,1))
plot(x = 1:5,
     y = sapply(results_MD.SI.A.5GCbins.mean_no_error, FUN = function(x){
       x$without_error$model1$bias
     }),
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(0, 5),
     axes = F,
     type = 'n')
text(x = 1.1, 
     y = 4.7,
     labels = c("A"),
     cex = 3)
points(x = 1:5,
       y = sapply(results_MD.SI.A.5GCbins.mean_no_error, FUN = function(x){
         x$without_error$model1$bias
       }),
       cex = 2)
title(main = 'simulans',
      line = 1,
      cex.main = 2)
mtext(expression(kappa),
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
     y = sapply(results_ZI69.SI.A.5GCbins.mean_no_error, FUN = function(x){
       x$without_error$model1$bias
     }),
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(0, 5),
     axes = F,
     type = 'n')
text(x = 1.1, 
     y = 4.7,
     labels = c("B"),
     cex = 3)
points(x = 1:5,
       y = sapply(results_ZI69.SI.A.5GCbins.mean_no_error, FUN = function(x){
         x$without_error$model1$bias
       }),
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
                      GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2, cex.axis = 1.4)
axis(2, cex.axis = 1.4)
box()


plot(x = 1:5,
     y = sapply(results_MD.SI.A.5GCbins.sim_no_error, FUN = function(x){
       x$without_error$model1$bias
     }),
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(0, 5),
     axes = F,
     type = 'n')
text(x = 1.1, 
     y = 4.7,
     labels = c("C"),
     cex = 3)
points(x = 1:5,
       y = sapply(results_MD.SI.A.5GCbins.sim_no_error, FUN = function(x){
         x$without_error$model1$bias
       }),
       cex = 2)
title(main = 'simulans',
      line = 1,
      cex.main = 2)
mtext(expression(kappa),
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
     y = sapply(results_ZI69.SI.A.5GCbins.mel_no_error, FUN = function(x){
       x$without_error$model1$bias
     }),
     xlab = '',
     ylab = '',
     main = '',
     ylim = c(0, 5),
     axes = F,
     type = 'n'
)
text(x = 1.1, 
     y = 4.7,
     labels = c("D"),
     cex = 3)
points(x = 1:5,
       y = sapply(results_ZI69.SI.A.5GCbins.mel_no_error, FUN = function(x){
         x$without_error$model1$bias
       }),
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

