# parse the recombination rate calculate output and sanity check it

RRCtemp <- read.table('~/Dropbox/Biology_projects/dros_X_A_mut/RRC/RRC_input.txt.rrc',
                      stringsAsFactors = F)
names(RRCtemp) <- c('info', 'RRCstart', 'RRCmid', 'RRCstop', 'Comstart', 'Commid', 'Comstop')

contig <- sapply(RRCtemp[,1], FUN = function(x){
  strsplit(x, ':')[[1]][1]
})

strt <- sapply(RRCtemp[,1], FUN = function(x){
  as.numeric(strsplit(strsplit(x, ':')[[1]][2], '\\.\\.')[[1]][1])
})

stp <- sapply(RRCtemp[,1], FUN = function(x){
  as.numeric(strsplit(strsplit(x, ':')[[1]][2], '\\.\\.')[[1]][2])
})

RRC <- cbind (contig, strt, stp, RRCtemp[,2:7])

# and the smoothed data:
RRCsmoothedtemp <- read.table('~/Dropbox/Biology_projects/dros_X_A_mut/RRC/RRC_input_smoothed.txt.rrc',
                      stringsAsFactors = F)
names(RRCsmoothedtemp) <- c('info', 'RRCstart', 'RRCmid', 'RRCstop', 'Comstart', 'Commid', 'Comstop')

contigsmoothed <- sapply(RRCsmoothedtemp[,1], FUN = function(x){
  strsplit(x, ':')[[1]][1]
})

strtsmoothed <- sapply(RRCsmoothedtemp[,1], FUN = function(x){
  as.numeric(strsplit(strsplit(x, ':')[[1]][2], '\\.\\.')[[1]][1])
})

stpsmoothed <- sapply(RRCsmoothedtemp[,1], FUN = function(x){
  as.numeric(strsplit(strsplit(x, ':')[[1]][2], '\\.\\.')[[1]][2])
})

RRCsmoothed <- cbind (contigsmoothed, strtsmoothed, stpsmoothed, RRCsmoothedtemp[,2:7])


# ##################################################
# # now sanity check the order of the introns in the table is the same as the order
# # of alignments in my list
# 
# load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_alignments.RData')
# 
# names.Dmel.SI.A <- sapply(files.intersection, FUN = function(x){
#   strsplit(strsplit(x, '/')[[1]][8], '&')[[1]][3]
# })
# names(names.Dmel.SI.A) <- NULL
# 
# myTable <- read.table('~/Dropbox/Biology_projects/dros_X_A_mut/data/table.txt',
#                       header = T)
# 
# for(i in seq_along(names.Dmel.SI.A)){
#   nm <- names.Dmel.SI.A[i]
#   mySubset <- myTable[myTable$dmelname == nm, 1:2]
#   mylower <- min(mySubset$melsite)
#   myupper <- max(mySubset$melsite)
#   rrc <- RRC[i,]
#   if(rrc$strt != mylower){print('nope')}
#   if(rrc$stp != myupper){print('nope')}
# }
# 
# for(i in seq_along(names.Dmel.SI.A)){
#   nm <- names.Dmel.SI.A[i]
#   mySubset <- myTable[myTable$dmelname == nm, 1:2]
#   mylower <- min(mySubset$melsite)
#   myupper <- max(mySubset$melsite)
#   rrc <- RRCsmoothed[i,]
#   if(rrc$strt != mylower){print('nope')}
#   if(rrc$stp != myupper){print('nope')}
# }
# 
# # Everything works
# 
# ###############################################################################

recomb.comeron.dmel.SI.A <- RRC$Commid
recomb.comeronSmoothed.dmel.SI.A <- RRCsmoothed$Commid

# save stuff
save(RRC, RRCsmoothed,
     recomb.comeron.dmel.SI.A, recomb.comeronSmoothed.dmel.SI.A,
     file = '~/Dropbox/Biology_projects/dros_X_A_mut/data/recombination_Dmel.RData')








#
