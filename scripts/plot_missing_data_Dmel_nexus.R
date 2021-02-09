# plot the proportion of missing data in the Dmel nexus

library(seqinr)

# the table of missing values doens't include samples that have no data masked
# for admixture
myTable <- read.table('~/data/dros_X_A_mut/Dmel_nexus/proportion_admixture_masked.txt',
                      header = T)

tempZIal <- read.alignment(file = '/Users/ben/data/dros_X_A_mut/fasta_SIs/ZI/2L.intron_FBgn0002121:4_FBgn0002121:3.fa',
                          format = 'fasta')
allZIlines <- grep('ZI', tempZIal$nam, value = T)

tempRGal <- read.alignment(file = '/Users/ben/data/dros_X_A_mut/fasta_SIs/RG/2L.intron_FBgn0002121:4_FBgn0002121:3.fa',
                           format = 'fasta')
allRGlines <- grep('RG', tempRGal$nam, value = T)

extraRG <- setdiff(allRGlines, myTable$sample)
extraZI <- setdiff(allZIlines, myTable$sample)

ZI_zeros <- data.frame('sample' = extraZI, 'proportion_masked' = rep(0, length(extraZI)))
RG_zeros <- data.frame('sample' = extraRG, 'proportion_masked' = rep(0, length(extraRG)))

ZItable <- merge(droplevels(subset(myTable, subset = grepl("^ZI", myTable$sample))),
                 ZI_zeros,
                 all = T)

RGtable <- merge(droplevels(subset(myTable, subset = grepl("^RG", myTable$sample))),
                 RG_zeros,
                 all = T)


pdf('~/data/dros_X_A_mut/plots/ZI_mask_proportion.pdf',
    height = 7,
    width = 28)  
x <- barplot(as.vector(ZItable$proportion_masked[order(ZItable$proportion_masked, decreasing = T)]),
        ylab = 'proportion masked')
axis(1,
     labels = ZItable$sample[order(ZItable$proportion_masked, decreasing = T)],
     at = x[,1],
     las = 2)
dev.off()

RG17 <- c('RG18N', 'RG19', 'RG2', 'RG22', 'RG24', 'RG25', 'RG28', 'RG3', 'RG32N', 'RG33', 'RG34', 'RG36', 'RG38N', 'RG4N', 'RG5', 'RG7', 'RG9')
grepIndices <- sapply(RG17, FUN = function(x){
  grep(paste0('^', x, '$'), RGtable$sample[order(RGtable$proportion_masked, decreasing = T)])
})
core17_indx <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
core17_indx[grepIndices] <- 2
pdf('~/data/dros_X_A_mut/plots/RG_mask_proportion.pdf',
    height = 7,
    width = 28)  
x <- barplot(as.vector(RGtable$proportion_masked[order(RGtable$proportion_masked, decreasing = T)]),
             ylab = 'proportion masked')
axis(1,
     labels = RGtable$sample[order(RGtable$proportion_masked, decreasing = T)][core17_indx == 1],
     at = x[,1][core17_indx == 1],
     las = 2,
     col.axis = 'black')
axis(1,
     labels = RGtable$sample[order(RGtable$proportion_masked, decreasing = T)][core17_indx == 2],
     at = x[,1][core17_indx == 2],
     las = 2,
     col.axis = 'red')
dev.off()
