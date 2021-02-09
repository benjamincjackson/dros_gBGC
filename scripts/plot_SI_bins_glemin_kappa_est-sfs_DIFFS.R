setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load the results
load('../data/glemin_results_DIFFS.RData')
# load the gc content and bin indices
load('../data/SI_binindices_GCcontent_DIFFERENCES.RData')


pdf('../plots/kappa_glemin_MD_SI_A_diffbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = sapply(results_MD.SI.A.5GCbins.diff_no_error, FUN = function(x){
       x$without_error$model1$bias
     }),
     main = 'MD (diff bins)',
     xlab = expression(GC[italic(sim)]-GC[italic(mel)]),
     ylab = expression(kappa),
     axes = F,
     pch = 1,
     ylim = c(0, 5))
points(x = 1:5,
       y = sapply(results_MD.SI.A.5GCbins.diff_error, FUN = function(x){
         x$with_error$model1$bias
       }),
       pch = 2)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      signif(aggregate(GCcontent.diff, by = list(GCbinindex.5bins.diff.SI.A), FUN = mean)$x,
                             2)),
     las = 2)
# legend('topleft',
#        bty = 'n',
#        legend = c(expression(gamma!=0),
#                   expression(gamma==0)),
#        pch = c(1, 2))
axis(2)
box()
dev.off()


pdf('../plots/kappa_glemin_ZI69_SI_A_diffbins.pdf')
par(lwd = 2,
    cex = 1.6)
plot(x = 1:5,
     y = sapply(results_ZI69.SI.A.5GCbins.diff_no_error, FUN = function(x){
       x$without_error$model1$bias
     }),
     main = 'ZI 69 lines (diff bins)',
     xlab = expression(GC[italic(sim)]-GC[italic(mel)]),
     ylab = expression(kappa),
     axes = F,
     pch = 1,
     ylim = c(0, 5))
points(x = 1:5,
       y = sapply(results_ZI69.SI.A.5GCbins.diff_error, FUN = function(x){
         x$with_error$model1$bias
       }),
       pch = 2)
axis(1, at = 1:5,
     labels = sprintf("%.2f",
                      signif(aggregate(GCcontent.diff, by = list(GCbinindex.5bins.diff.SI.A), FUN = mean)$x,
                             2)),
     las = 2)
# legend('topleft',
#        bty = 'n',
#        legend = c(expression(gamma!=0),
#                   expression(gamma==0)),
#        pch = c(1, 2))
axis(2)
box()
dev.off()

