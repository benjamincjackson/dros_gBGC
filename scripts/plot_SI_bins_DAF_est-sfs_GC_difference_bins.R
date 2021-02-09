# plot DAFs defining ancestral sites using Peter's method

# load the matrices that we used to write the est-sfs input files with:
load('../data/est-sfs_matrices_DIFFS.RData')
load('../data/est-sfs_SFSs.RData')

# load the bin indices and GC contents
load('../data/SI_binindices_GCcontent_DIFFERENCES.RData')


library(RColorBrewer)
col <- brewer.pal(n = 3, name = 'Dark2')
pnt <- c(1, 2, 0)

pdf('~/Dropbox/Biology_projects/dros_X_A_mut/plots/DAF_SI_A_MD_est-sfs_diffbins.pdf')
plot(AT_to_GC_ALL_DAF.MD.diff,
     main = 'MD (diff bins)',
     xlab = 'GC content difference bin',
     ylab = 'derived allele frequency',
     ylim = c(0, 0.35),
     axes = F,
     col = col[1],
     pch = pnt[1])
axis(1, at = 1:5,
     labels = 1:5,
     las = 2)
axis(2)
box()
points(GC_to_AT_ALL_DAF.MD.diff, col = col[2], pch = pnt[2])
points(N_to_N_ALL_DAF.MD.diff, col = col[3], pch = pnt[3])
legend('topleft',
       legend = c(expression(DAF[W %->% S]), expression(DAF[neu]),
                  expression(DAF[S %->% W])),
       col = c(col[c(1,3,2)]),
       pch = c(pnt[c(1,3,2)]),
       bty = 'n')
dev.off()

#######

pdf('~/Dropbox/Biology_projects/dros_X_A_mut/plots/DAF_SI_A_ZI69_est-sfs_diffbins.pdf')
plot(AT_to_GC_ALL_DAF.ZI69.diff,
     main = 'ZI 69 lines (diff bins)',
     xlab = 'GC content difference bin',
     ylab = 'derived allele frequency',
     ylim = c(0, 0.35),
     axes = F,
     col = col[1],
     pch = pnt[1])
axis(1, at = 1:5,
     labels =  1:5,
     las = 2)
axis(2)
box()
points(GC_to_AT_ALL_DAF.ZI69.diff, col = col[2], pch = pnt[2])
points(N_to_N_ALL_DAF.ZI69.diff, col = col[3], pch = pnt[3])
legend('topleft',
       legend = c(expression(DAF[W %->% S]), expression(DAF[neu]),
                  expression(DAF[S %->% W])),
       col = c(col[c(1,3,2)]),
       pch = c(pnt[c(1,3,2)]),
       bty = 'n')
dev.off()




#
