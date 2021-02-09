setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load the results
load('../data/glemin_results.RData')
# load the gc content and bin indices
load('../data/SI_binindices_GCcontent.RData')


pdf('../plots/gamma_glemin_MD_SI_A_meanbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = sapply(results_MD.SI.A.5GCbins.mean_error, FUN = function(x){
             x$with_error$model1$B
           }),
     main = 'MD (mean bins)',
     xlab = 'GC content',
     ylab = expression(gamma),
     axes = F,
     pch = c(1, 16, 16, 16, 16),
     ylim = c(-1, 3))
points(x = 1:5,
       y = sapply(results_MD.SI.A.5GCbins.mean_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(17, 17, 17, 17, 17))
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.sim.meanbins.SI.A),
     las = 2)
legend('topleft',
       bty = 'n',
       legend = c('accounting for polarisation error',
                  'not accounting for polarisation error',
                  'filled symbols - significantly different from 0',
                  'open symbols - not significanttly different from 0'),
       pch = c(16, 17, NA, NA),
       cex = 0.6)
axis(2)
box()
dev.off()

pdf('../plots/gamma_glemin_MD_SI_A_simbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = sapply(results_MD.SI.A.5GCbins.sim_error, FUN = function(x){
       x$with_error$model1$B
     }),
     main = 'MD (sim bins)',
     xlab = 'GC content',
     ylab = expression(gamma),
     axes = F,
     pch = c(1, 16, 16, 16, 16),
     ylim = c(-1, 3))
points(x = 1:5,
       y = sapply(results_MD.SI.A.5GCbins.sim_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(2, 17, 17, 17, 17))
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.sim.simbins.SI.A),
     las = 2)
legend('topleft',
       bty = 'n',
       legend = c('accounting for polarisation error',
                  'not accounting for polarisation error',
                  'filled symbols - significantly different from 0',
                  'open symbols - not significanttly different from 0'),
       pch = c(16, 17, NA, NA),
       cex = 0.6)
axis(2)
box()
dev.off()

pdf('../plots/gamma_glemin_ZI69_SI_A_meanbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = sapply(results_ZI69.SI.A.5GCbins.mean_error, FUN = function(x){
       x$with_error$model1$B
     }),
     main = 'ZI 69 lines (mean bins)',
     xlab = 'GC content',
     ylab = expression(gamma),
     axes = F,
     pch = c(1, 1, 16, 16, 16),
     ylim = c(-1, 3))
points(x = 1:5,
       y = sapply(results_ZI69.SI.A.5GCbins.mean_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(2, 2, 17, 17, 17))
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2)
legend('topleft',
       bty = 'n',
       legend = c('accounting for polarisation error',
                  'not accounting for polarisation error',
                  'filled symbols - significantly different from 0',
                  'open symbols - not significantly different from 0'),
       pch = c(16, 17, NA, NA),
       cex = 0.6)
axis(2)
box()
dev.off()

pdf('../plots/gamma_glemin_ZI69_SI_A_melbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = sapply(results_ZI69.SI.A.5GCbins.mel_error, FUN = function(x){
       x$with_error$model1$B
     }),
     main = 'ZI 69 lines (mel bins)',
     xlab = 'GC content',
     ylab = expression(gamma),
     axes = F,
     pch = c(16, 1, 16, 16, 16),
     ylim = c(-1, 3))
points(x = 1:5,
       y = sapply(results_ZI69.SI.A.5GCbins.mel_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(17, 2, 17, 17, 17))
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.melbins.SI.A),
     las = 2)
legend('topleft',
       bty = 'n',
       legend = c('accounting for polarisation error',
                  'not accounting for polarisation error',
                  'filled symbols - significantly different from 0',
                  'open symbols - not significantly different from 0'),
       pch = c(16, 17, NA, NA),
       cex = 0.6)
axis(2)
box()
dev.off()

pdf('../plots/gamma_glemin_ZI21_SI_A_meanbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = sapply(results_ZI21.SI.A.5GCbins.mean_error, FUN = function(x){
       x$with_error$model1$B
     }),
     main = 'ZI 21 lines, mean bins',
     xlab = 'GC content',
     ylab = expression(gamma),
     axes = F,
     pch = c(1, 1, 16, 16, 16),
     ylim = c(-1, 3))
points(x = 1:5,
       y = sapply(results_ZI21.SI.A.5GCbins.mean_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(2, 2, 17, 17, 17))
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2)
legend('topleft',
       bty = 'n',
       legend = c('accounting for polarisation error',
                  'not accounting for polarisation error',
                  'filled symbols - significantly different from 0',
                  'open symbols - not significantly different from 0'),
       pch = c(16, 17, NA, NA),
       cex = 0.6)
axis(2)
box()
dev.off()


pdf('../plots/gamma_glemin_RG21_SI_A_meanbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = sapply(results_RG21.SI.A.5GCbins.mean_error, FUN = function(x){
       x$with_error$model1$B
     }),
     main = 'RG 21 lines',
     xlab = 'maximum GC content',
     ylab = expression(gamma),
     axes = F,
     pch = c(16, 1, 1, 16, 16),
     ylim = c(-1, 3))
points(x = 1:5,
       y = sapply(results_RG21.SI.A.5GCbins.mean_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(2, 2, 17, 17, 17))
axis(1, at = 1:5,
     labels =  sprintf("%.2f",
                       GCcontent.concat.5bins.mel.meanbins.SI.A),
     las = 2)
legend('topleft',
       bty = 'n',
       legend = c('accounting for polarisation error',
                  'not accounting for polarisation error',
                  'filled symbols - significantly different from 0',
                  'open symbols - not significantly different from 0'),
       pch = c(16, 17, NA, NA),
       cex = 0.6)
axis(2)
box()
dev.off()






#