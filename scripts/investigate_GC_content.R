# let's investigate the distribution of bin indices of the opposing species' bin
# with the mean compared to the species bins. How?

# load the bin index and GC content data
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_binindices_GCcontent.RData')

x <- GCcontent.Dsim.SI.A - GCcontent.Dmel.SI.A
hist(x, breaks = 24)
plot(x)

mean(GCcontent.Dsim.SI.A)
mean(GCcontent.Dmel.SI.A)

order(x)
# these are the introns that have moved highest towards GC in mel:
x[order(x)][1:1133]

# there are the introns that have moved highest towards GC in sim:
x[order(x, decreasing = T)][1:1695]
