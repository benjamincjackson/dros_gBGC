setwd(dirname(rstudioapi::getSourceEditorContext()$path))

load('../data/est-sfs_tables.RData')
load('../data/est-sfs_matrices.RData')
load('../data/SI_binindices_GCcontent.RData')
load('../data/SI_alignments.RData')
source('../scripts/R_helper_functions.R')
library('seqinr')

# concatenate the sequences according to GC content bin:
temp.concat.alignments.GCbin.mean.MD.SI.A <- get_concat_seqs(alignments.intersection.MD21.SI.A,
                                                             GCbinindex.5bins.mean.SI.A)
temp.concat.alignments.GCbin.mean.ZI69.SI.A <- get_concat_seqs(alignments.intersection.ZI69.SI.A,
                                                               GCbinindex.5bins.mean.SI.A)
temp.concat.alignments.GCbin.mean.ZI21.SI.A <- get_concat_seqs(alignments.intersection.ZI21.SI.A,
                                                               GCbinindex.5bins.mean.SI.A)
temp.concat.alignments.GCbin.mean.RG21.SI.A <- get_concat_seqs(alignments.intersection.RG21.SI.A,
                                                               GCbinindex.5bins.mean.SI.A)

temp.concat.alignments.GCbin.sim.MD.SI.A <- get_concat_seqs(alignments.intersection.MD21.SI.A,
                                                            GCbinindex.5bins.Dsim.SI.A)

temp.concat.alignments.GCbin.mel.ZI69.SI.A <- get_concat_seqs(alignments.intersection.ZI69.SI.A,
                                                              GCbinindex.5bins.Dmel.SI.A)
temp.concat.alignments.GCbin.mel.ZI21.SI.A <- get_concat_seqs(alignments.intersection.ZI21.SI.A,
                                                              GCbinindex.5bins.Dmel.SI.A)
temp.concat.alignments.GCbin.mel.RG21.SI.A <- get_concat_seqs(alignments.intersection.RG21.SI.A,
                                                              GCbinindex.5bins.Dmel.SI.A)

get_indx <- function(alignment){
  # make a matrix of the alignment
  mat <- makeMatrix(alignment)
  # maybe the first three rows should be sim, mel yak references
  # then the next x (21) lines should be the polymophism sample:
  linesToKeep <- 'MD|ZI|RG'
  newMat <- rbind(mat[grep('dsim', alignment$nam),],
                  mat[grep('dmel', alignment$nam),],
                  mat[grep('dyak', alignment$nam),],
                  mat[grep(linesToKeep, alignment$nam),])
  # so the sim ref is row 1. the mel and yak refs are rows 2 & 3.
  # and these are the polymophism rows:
  polyRows <- 4:(length(grep(linesToKeep, alignment$nam)) + 3)
  # now retain the sites where there are at most two alleles in the polymorphism dataset AND
  # no missind data:
  indx <- apply(newMat, 2, FUN = function(cl){
    all(nlevels(as.factor(cl[polyRows])) <= 2,
        all(levels(as.factor(cl)) != 'n'))
  })
  return(indx)
}

length.alignments.intersection <- sapply(alignments.intersection, FUN = function(x){
  nchar(x$seq[1])
})

files <- sapply(files.intersection, FUN = function(x){
  strsplit(x, split = '/')[[1]][4]
})
names(files) <- NULL

files.sites <- cbind(files, length.alignments.intersection)

for(i in 1:5){
  temp <- files.sites[GCbinindex.5bins.Dsim.SI.A == i,]
  colnames(temp) <- NULL
  temp2 <- unlist(apply(temp, 1,
                        FUN = function(rw){
                          rep(rw[1], as.integer(rw[2]))
                        }))
  indx <- get_indx(temp.concat.alignments.GCbin.sim.MD.SI.A[[i]])
  temp3 <- temp2[indx]
  assign(x = paste0('files.sites.MD21.sim.bin0', i, collapse = ''),
         value = temp3
  )

  temp <- files.sites[GCbinindex.5bins.Dmel.SI.A == i,]
  colnames(temp) <- NULL
  temp2 <- unlist(apply(temp, 1,
                        FUN = function(rw){
                          rep(rw[1], as.integer(rw[2]))
                        }))
  indx <- get_indx(temp.concat.alignments.GCbin.mel.RG21.SI.A[[i]])
  temp3 <- temp2[indx]
  assign(x = paste0('files.sites.RG21.mel.bin0', i, collapse = ''),
         value = temp3
  )

  temp <- files.sites[GCbinindex.5bins.Dmel.SI.A == i,]
  colnames(temp) <- NULL
  temp2 <- unlist(apply(temp, 1,
                        FUN = function(rw){
                          rep(rw[1], as.integer(rw[2]))
                        }))
  indx <- get_indx(temp.concat.alignments.GCbin.mel.ZI21.SI.A[[i]])
  temp3 <- temp2[indx]
  assign(x = paste0('files.sites.ZI21.mel.bin0', i, collapse = ''),
         value = temp3
  )

  temp <- files.sites[GCbinindex.5bins.Dmel.SI.A == i,]
  colnames(temp) <- NULL
  temp2 <- unlist(apply(temp, 1,
                        FUN = function(rw){
                          rep(rw[1], as.integer(rw[2]))
                        }))
  indx <- get_indx(temp.concat.alignments.GCbin.mel.ZI69.SI.A[[i]])
  temp3 <- temp2[indx]
  assign(x = paste0('files.sites.ZI69.mel.bin0', i, collapse = ''),
         value = temp3
  )

  temp <- files.sites[GCbinindex.5bins.mean.SI.A == i,]
  colnames(temp) <- NULL
  temp2 <- unlist(apply(temp, 1,
                        FUN = function(rw){
                          rep(rw[1], as.integer(rw[2]))
                        }))
  indx <- get_indx(temp.concat.alignments.GCbin.mean.MD.SI.A[[i]])
  temp3 <- temp2[indx]
  assign(x = paste0('files.sites.MD21.mean.bin0', i, collapse = ''),
         value = temp3
  )

  temp <- files.sites[GCbinindex.5bins.mean.SI.A == i,]
  colnames(temp) <- NULL
  temp2 <- unlist(apply(temp, 1,
                        FUN = function(rw){
                          rep(rw[1], as.integer(rw[2]))
                        }))
  indx <- get_indx(temp.concat.alignments.GCbin.mean.RG21.SI.A[[i]])
  temp3 <- temp2[indx]
  assign(x = paste0('files.sites.RG21.mean.bin0', i, collapse = ''),
         value = temp3
  )

  temp <- files.sites[GCbinindex.5bins.mean.SI.A == i,]
  colnames(temp) <- NULL
  temp2 <- unlist(apply(temp, 1,
                        FUN = function(rw){
                          rep(rw[1], as.integer(rw[2]))
                        }))
  indx <- get_indx(temp.concat.alignments.GCbin.mean.ZI21.SI.A[[i]])
  temp3 <- temp2[indx]
  assign(x = paste0('files.sites.ZI21.mean.bin0', i, collapse = ''),
         value = temp3
  )

  temp <- files.sites[GCbinindex.5bins.mean.SI.A == i,]
  colnames(temp) <- NULL
  temp2 <- unlist(apply(temp, 1,
                        FUN = function(rw){
                          rep(rw[1], as.integer(rw[2]))
                        }))
  indx <- get_indx(temp.concat.alignments.GCbin.mean.ZI69.SI.A[[i]])
  temp3 <- temp2[indx]
  assign(x = paste0('files.sites.ZI69.mean.bin0', i, collapse = ''),
         value = temp3
  )

  rm(temp, temp2, temp3, indx)
}

popSFSs <- function(rw, n){
  newrw <- unname(unlist(rw))
  SFS_AT_2_AT <- vector('numeric', length = n + 1)
  SFS_AT_2_GC <- vector('numeric', length = n + 1)
  SFS_GC_2_GC <- vector('numeric', length = n + 1)
  SFS_GC_2_AT <- vector('numeric', length = n + 1)

  # counts of alleles in the poly data:
  cA <- newrw[8]
  cC <-  newrw[9]
  cG <-  newrw[10]
  cT <-  newrw[11]

  countVector <- c('A' = cA, 'C' = cC, 'G' = cG, 'T' = cT)

  # probability of the major allele being ancestral:
  pMajAnc <- newrw[3]
  # the identity of the major allele:
  majAl <- names(which.max(countVector))
  minAl <- names(countVector)[order(countVector, decreasing = T)[2]]

  # populate the SFSs:
  if(majAl == 'A'){
    # the sfs c olumn to populate if A is ancestral:
    indx_anc <- min(cA + 1, n + 1 - cA)
    # the sfs column to populate if A is derived:
    indx_dev <- max(cA + 1, n + 1 - cA)

    if (minAl == 'C'){
      # either this is C -> A or A -> C
      SFS_AT_2_GC[indx_anc] <- SFS_AT_2_GC[indx_anc] + 1 * pMajAnc
      SFS_GC_2_AT[indx_dev] <- SFS_GC_2_AT[indx_dev] + 1 * (1 - pMajAnc)
    }
    if (minAl == 'G'){
      SFS_AT_2_GC[indx_anc] <- SFS_AT_2_GC[indx_anc] + 1 * pMajAnc
      SFS_GC_2_AT[indx_dev] <- SFS_GC_2_AT[indx_dev] + 1 * (1 - pMajAnc)
    }
    if (minAl == 'T'){
      SFS_AT_2_AT[indx_anc] <- SFS_AT_2_AT[indx_anc] + 1 * pMajAnc
      SFS_AT_2_AT[indx_dev] <- SFS_AT_2_AT[indx_dev] + 1 * (1 - pMajAnc)
    }
  }

  if(majAl == 'C'){
    indx_anc <- min(cC + 1, n + 1 - cC)
    indx_dev <- max(cC + 1, n + 1 - cC)

    if (minAl == 'A'){
      SFS_GC_2_AT[indx_anc] <- SFS_GC_2_AT[indx_anc] + 1 * pMajAnc
      SFS_AT_2_GC[indx_dev] <- SFS_AT_2_GC[indx_dev] + 1 * (1 - pMajAnc)
    }
    if (minAl == 'G'){
      SFS_GC_2_GC[indx_anc] <- SFS_GC_2_GC[indx_anc] + 1 * pMajAnc
      SFS_GC_2_GC[indx_dev] <- SFS_GC_2_GC[indx_dev] + 1 * (1 - pMajAnc)
    }
    if (minAl == 'T'){
      SFS_GC_2_AT[indx_anc] <- SFS_GC_2_AT[indx_anc] + 1 * pMajAnc
      SFS_AT_2_GC[indx_dev] <- SFS_AT_2_GC[indx_dev] + 1 * (1 - pMajAnc)
    }
  }

  if(majAl == 'G'){
    indx_anc <- min(cG + 1, n + 1 - cG)
    indx_dev <- max(cG + 1, n + 1 - cG)

    if (minAl == 'A'){
      SFS_GC_2_AT[indx_anc] <- SFS_GC_2_AT[indx_anc] + 1 * pMajAnc
      SFS_AT_2_GC[indx_dev] <- SFS_AT_2_GC[indx_dev] + 1 * (1 - pMajAnc)
    }
    if (minAl == 'C'){
      SFS_GC_2_GC[indx_anc] <- SFS_GC_2_GC[indx_anc] + 1 * pMajAnc
      SFS_GC_2_GC[indx_dev] <- SFS_GC_2_GC[indx_dev] + 1 * (1 - pMajAnc)
    }
    if (minAl == 'T'){
      SFS_GC_2_AT[indx_anc] <- SFS_GC_2_AT[indx_anc] + 1 * pMajAnc
      SFS_AT_2_GC[indx_dev] <- SFS_AT_2_GC[indx_dev] + 1 * (1 - pMajAnc)
    }
  }

  if(majAl == 'T'){
    indx_anc <- min(cT + 1, n + 1 - cT)
    indx_dev <- max(cT + 1, n + 1 - cT)

    if (minAl == 'A'){
      SFS_AT_2_AT[indx_anc] <- SFS_AT_2_AT[indx_anc] + 1 * pMajAnc
      SFS_AT_2_AT[indx_dev] <- SFS_AT_2_AT[indx_dev] + 1 * (1 - pMajAnc)
    }
    if (minAl == 'C'){
      SFS_AT_2_GC[indx_anc] <- SFS_AT_2_GC[indx_anc] + 1 * pMajAnc
      SFS_GC_2_AT[indx_dev] <- SFS_GC_2_AT[indx_dev] + 1 * (1 - pMajAnc)
    }
    if (minAl == 'G'){
      SFS_AT_2_GC[indx_anc] <- SFS_AT_2_GC[indx_anc] + 1 * pMajAnc
      SFS_GC_2_AT[indx_dev] <- SFS_GC_2_AT[indx_dev] + 1 * (1 - pMajAnc)
    }
  }

  return(list('SFS_AT_2_AT' = SFS_AT_2_AT,
              'SFS_AT_2_GC' = SFS_AT_2_GC,
              'SFS_GC_2_AT' = SFS_GC_2_AT,
              'SFS_GC_2_GC' = SFS_GC_2_GC))
}


getDAFs <- function(estsfsTable){
  n <- sum(estsfsTable[1,8:11])

  SFS_AT_2_AT <- vector('numeric', length = n + 1)
  SFS_AT_2_GC <- vector('numeric', length = n + 1)
  SFS_GC_2_GC <- vector('numeric', length = n + 1)
  SFS_GC_2_AT <- vector('numeric', length = n + 1)

  for(i in seq_len(nrow(estsfsTable))){
    results <- NULL
    rw <- estsfsTable[i,]
    newrw <- unname(unlist(rw))
    countVector <- c(newrw[8:11])
    # if this site isn't segregating, then skip it:
    if(identical(sum(countVector != 0), 1L)){
      next
    }
    else{
      results <- popSFSs(rw, n)
    }

    SFS_AT_2_AT <- SFS_AT_2_AT + results$SFS_AT_2_AT
    SFS_AT_2_GC <- SFS_AT_2_GC + results$SFS_AT_2_GC
    SFS_GC_2_AT <- SFS_GC_2_AT + results$SFS_GC_2_AT
    SFS_GC_2_GC <- SFS_GC_2_GC + results$SFS_GC_2_GC
  }

  SFS_N_2_N <-  SFS_AT_2_AT + SFS_GC_2_GC

  # get the derived allele frequencies from the unfolded SFSs:
  DAF_AT_2_GC <- sum(SFS_AT_2_GC[2:n] * 1:(n-1)/n) /  sum(SFS_AT_2_GC[2:n])
  DAF_GC_2_AT <- sum(SFS_GC_2_AT[2:n] * 1:(n-1)/n) /  sum(SFS_GC_2_AT[2:n])
  DAF_N_2_N <-  sum(SFS_N_2_N[2:n] * 1:(n-1)/n) /  sum(SFS_N_2_N[2:n])

  return(data.frame('DAF_N_2_N' = DAF_N_2_N,
                    'DAF_AT_2_GC' = DAF_AT_2_GC,
                    'DAF_GC_2_AT' = DAF_GC_2_AT))
}

# a function to get bootstrapped DAFs by resampling introns from an est-sfs
# output table
getDAFCIs <- function(estSFStable, filenameVec, nbootstraps){
  # for each bootstrap, make a new table by sampling introns
  splittable <- split(estSFStable, filenameVec)
  DAFs <- NULL
  for(i in seq_len(nbootstraps)){
    newtable <- do.call(rbind,
                        sample(splittable, size = length(splittable), replace = T))
    DAFs <- rbind(DAFs, getDAFs(newtable))

  }
  return(DAFs)
}


for(i in 1:5){
  assign(paste0('DAFCIs.MD21.A.sim.GCbin0', i, collapse = ''),
         getDAFCIs(get(paste0('estsfsTable.MD.A.sim.GCbin0', i, collapse = NULL)),
                   get(paste0('files.sites.MD21.sim.bin0', i, collapse = NULL)),
                   1000))
}
for(i in 1:5){
  assign(paste0('DAFCIs.MD21.A.mean.GCbin0', i, collapse = ''),
         getDAFCIs(get(paste0('estsfsTable.MD.A.mean.GCbin0', i, collapse = NULL)),
                   get(paste0('files.sites.MD21.mean.bin0', i, collapse = NULL)),
                   1000))
}
# for(i in 1:5){
#   assign(paste0('DAFCIs.RG21.A.mel.GCbin0', i, collapse = ''),
#          getDAFCIs(get(paste0('estsfsTable.RG21.A.mel.GCbin0', i, collapse = NULL)),
#                    get(paste0('files.sites.RG21.mel.bin0', i, collapse = NULL)),
#                    1000))
# }
# for(i in 1:5){
#     assign(paste0('DAFCIs.RG21.A.mean.GCbin0', i, collapse = ''),
#            getDAFCIs(get(paste0('estsfsTable.RG21.A.mean.GCbin0', i, collapse = NULL)),
#                      get(paste0('files.sites.RG21.mean.bin0', i, collapse = NULL)),
#                      1000))
  }
# for(i in 1:5){
#   assign(paste0('DAFCIs.ZI21.A.mel.GCbin0', i, collapse = ''),
#          getDAFCIs(get(paste0('estsfsTable.ZI21.A.mel.GCbin0', i, collapse = NULL)),
#                    get(paste0('files.sites.ZI21.mel.bin0', i, collapse = NULL)),
#                    1000))
# }
# for(i in 1:5){
#   assign(paste0('DAFCIs.ZI21.A.mean.GCbin0', i, collapse = ''),
#          getDAFCIs(get(paste0('estsfsTable.ZI21.A.mean.GCbin0', i, collapse = NULL)),
#                    get(paste0('files.sites.ZI21.mean.bin0', i, collapse = NULL)),
#                    1000))
# }
for(i in 1:5){
  assign(paste0('DAFCIs.ZI69.A.mel.GCbin0', i, collapse = ''),
         getDAFCIs(get(paste0('estsfsTable.ZI69.A.mel.GCbin0', i, collapse = NULL)),
                   get(paste0('files.sites.ZI69.mel.bin0', i, collapse = NULL)),
                   1000))
}
for(i in 1:5){
  assign(paste0('DAFCIs.ZI69.A.mean.GCbin0', i, collapse = ''),
         getDAFCIs(get(paste0('estsfsTable.ZI69.A.mean.GCbin0', i, collapse = NULL)),
                   get(paste0('files.sites.ZI69.mean.bin0', i, collapse = NULL)),
                   1000))
}

save(list = ls(pattern = 'DAFCIs.'),
     file = '../data/DAF_CI_bootstraps.RData')
