setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load the results
load('../data/glemin_results_DIFFS.RData')
# load the gc content and bin indices
load('../data/SI_binindices_GCcontent_DIFFERENCES.RData')


pdf('../plots/gamma_glemin_MD_SI_A_diffbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = sapply(results_MD.SI.A.5GCbins.diff_error, FUN = function(x){
       x$with_error$model1$B
     }),
     main = 'MD (diff bins)',
     xlab = expression(GC[italic(sim)]-GC[italic(mel)]),
     ylab = expression(gamma),
     axes = F,
     pch = c(16, 1, 16, 16, 16),
     ylim = c(-2, 4))
points(x = 1:5,
       y = sapply(results_MD.SI.A.5GCbins.diff_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(17, 2, 17, 17, 17))
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      signif(aggregate(GCcontent.diff, by = list(GCbinindex.5bins.diff.SI.A), FUN = mean)$x,
                             2)),
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


pdf('../plots/gamma_glemin_ZI69_SI_A_diffbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = sapply(results_ZI69.SI.A.5GCbins.diff_error, FUN = function(x){
       x$with_error$model1$B
     }),
     main = 'ZI 69 lines (diff bins)',
     xlab = expression(GC[italic(sim)]-GC[italic(mel)]),
     ylab = expression(gamma),
     axes = F,
     pch = c(16, 16, 16, 16, 16),
     ylim = c(-2, 4))
points(x = 1:5,
       y = sapply(results_ZI69.SI.A.5GCbins.diff_no_error, FUN = function(x){
         x$without_error$model1$B
       }),
       pch = c(17, 17, 17, 17, 17))
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      signif(aggregate(GCcontent.diff, by = list(GCbinindex.5bins.diff.SI.A), FUN = mean)$x,
                             2)),
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
