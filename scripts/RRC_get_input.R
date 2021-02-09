# write an input file for the drosophila recombination rate calculator.
# should be in the format:
# 2L:7500..9400
# 2L:210918..250151
# ...

load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_alignments.RData')

names.Dmel.SI.A <- sapply(files.intersection, FUN = function(x){
  strsplit(strsplit(x, '/')[[1]][8], '&')[[1]][3]
})
names(names.Dmel.SI.A) <- NULL

myTable <- read.table('~/Dropbox/Biology_projects/dros_X_A_mut/data/table.txt',
                      header = T)

# levels(as.factor(as.character(myTable[myTable$dmelname == names.Dmel.SI.A[1], 1:2]$melcontig)))
# max(myTable[myTable$dmelname == names.Dmel.SI.A[1], 1:2]$melsite)
# min(myTable[myTable$dmelname == names.Dmel.SI.A[1], 1:2]$melsite)
# median(myTable[myTable$dmelname == names.Dmel.SI.A[1], 1:2]$melsite)

system.time(
linesToWrite <- sapply(names.Dmel.SI.A, FUN = function(nm){
  mySubset <- myTable[myTable$dmelname == nm, 1:2]
  if(nlevels(as.factor(as.character(mySubset$melcontig))) != 1){
    print('more than one contig')
  }
  if(max(mySubset$melsite) - min(mySubset$melsite) > 23){
    print('range greater than 23')
  }
  contig <- levels(as.factor(as.character(mySubset$melcontig)))
  a <- min(mySubset$melsite)
  b <- max(mySubset$melsite)
  return(paste0(contig, ':', a, '..', b, collapse = ''))
})
)

fileConn <- file('~/Dropbox/Biology_projects/dros_X_A_mut/RRC/RRC_input.txt')
writeLines(text=paste0(linesToWrite, collapse = '\n'),
           con = fileConn)
close(fileConn)  




#