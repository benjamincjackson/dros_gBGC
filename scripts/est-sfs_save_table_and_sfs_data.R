setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load the matrices that we used to write the est-sfs input files with:
load('../data/est-sfs_matrices.RData')
load('../data/est-sfs_matrices_DIFFS.RData')

# load the bin indices and GC contents
load('../data/SI_binindices_GCcontent.RData')

# read in the persite output files from est-sfs
for(i in 1:5){
  assign(paste0('estsfsTable.MD.mean.GCbin', i),
         read.table(paste0('../est-sfs/mean/MD_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}
for(i in 1:5){
  assign(paste0('estsfsTable.MD.sim.GCbin', i),
         read.table(paste0('../est-sfs/sim/MD_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}
estsfsTable.MD.ALL <- read.table('../est-sfs/all/MD_SI_A/output_files/ALL.pvalues',
                                 skip = 8)

for(i in 1:5){
  assign(paste0('estsfsTable.ZI69.mean.GCbin', i),
         read.table(paste0('../est-sfs/mean/ZI69_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}
for(i in 1:5){
  assign(paste0('estsfsTable.ZI69.mel.GCbin', i),
         read.table(paste0('../est-sfs/mel/ZI69_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}
estsfsTable.ZI69.ALL <- read.table('../est-sfs/all/ZI69_SI_A/output_files/ALL.pvalues',
                                  skip = 8)

for(i in 1:5){
  assign(paste0('estsfsTable.ZI21.mean.GCbin', i),
         read.table(paste0('../est-sfs/mean/ZI21_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}
for(i in 1:5){
  assign(paste0('estsfsTable.ZI21.mel.GCbin', i),
         read.table(paste0('../est-sfs/mel/ZI21_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}
estsfsTable.ZI21.ALL <- read.table('../est-sfs/all/ZI21_SI_A/output_files/ALL.pvalues',
                                   skip = 8)

for(i in 1:5){
  assign(paste0('estsfsTable.RG21.mean.GCbin', i),
         read.table(paste0('../est-sfs/mean/RG21_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}
for(i in 1:5){
  assign(paste0('estsfsTable.RG21.mel.GCbin', i),
         read.table(paste0('../est-sfs/mel/RG21_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}
estsfsTable.RG21.ALL <- read.table('../est-sfs/all/RG21_SI_A/output_files/ALL.pvalues',
                                   skip = 8)

for(i in 1:5){
  assign(paste0('estsfsTable.MD.diff.GCbin', i),
         read.table(paste0('../est-sfs/diff/MD_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}
for(i in 1:5){
  assign(paste0('estsfsTable.ZI69.diff.GCbin', i),
         read.table(paste0('../est-sfs/diff/ZI69_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}

###############################################
 
countNucs <- function(cl){
  LUTindex <- c('A' = 1L, 'C' = 2L, 'G' = 3L, 'T' = 4L)
  counts <- rep(0, 4)
  for(i in seq_along(cl)){
    counts[LUTindex[toupper(cl[i])]] <- counts[LUTindex[toupper(cl[i])]] + 1
  }
  return(counts)
}

# a function to return the allelic counts of the polymorphism sample
getPolyCounts <- function(mat){
  # just get the polymorphism rows, if there are three outgroup rows first
  polyMat <- mat[4:nrow(mat),]
  counts <- apply(polyMat, 2, countNucs)
  return(t(counts))
}


# get tables of counts of the alleles for each bin using the function above
for(i in 1:5){
  assign(paste0('polyCounts.MD.SI.A.mean.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.mean.GCbin0', i, '.MD.A'))))
}
for(i in 1:5){
  assign(paste0('polyCounts.MD.SI.A.sim.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.sim.GCbin0', i, '.MD.A'))))
}

for(i in 1:5){
  assign(paste0('polyCounts.ZI69.SI.A.mean.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.mean.GCbin0', i, '.ZI69.A'))))
}
for(i in 1:5){
  assign(paste0('polyCounts.ZI69.SI.A.mel.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.mel.GCbin0', i, '.ZI69.A'))))
}

for(i in 1:5){
  assign(paste0('polyCounts.ZI21.SI.A.mean.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.mean.GCbin0', i, '.ZI21.A'))))
}
for(i in 1:5){
  assign(paste0('polyCounts.ZI21.SI.A.mel.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.mel.GCbin0', i, '.ZI21.A'))))
}

for(i in 1:5){
  assign(paste0('polyCounts.RG21.SI.A.mean.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.mean.GCbin0', i, '.RG21.A'))))
}
for(i in 1:5){
  assign(paste0('polyCounts.RG21.SI.A.mel.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.mel.GCbin0', i, '.RG21.A'))))
}
for(i in 1:5){
  assign(paste0('polyCounts.MD.SI.A.diff.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.diff.GCbin0', i, '.MD.A'))))
}
for(i in 1:5){
  assign(paste0('polyCounts.ZI69.SI.A.diff.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.diff.GCbin0', i, '.ZI69.A'))))
}

polyCounts.MD.SI.A.ALL <- getPolyCounts(estsfsMat.ALL.MD.A)
polyCounts.ZI69.SI.A.ALL <- getPolyCounts(estsfsMat.ALL.ZI69.A)
polyCounts.ZI21.SI.A.ALL <- getPolyCounts(estsfsMat.ALL.ZI21.A)
polyCounts.RG21.SI.A.ALL <- getPolyCounts(estsfsMat.ALL.RG21.A)


# stick these tables together with the probability tables, which we can
# work with to calculate DAF

for(i in 1:5){
  assign(paste0('estsfsTable.MD.A.mean.GCbin0', i),
         cbind(get(paste0('estsfsTable.MD.mean.GCbin', i)), get(paste0('polyCounts.MD.SI.A.mean.GCbin0', i))))
  assign(paste0('estsfsTable.MD.A.sim.GCbin0', i),
         cbind(get(paste0('estsfsTable.MD.sim.GCbin', i)), get(paste0('polyCounts.MD.SI.A.sim.GCbin0', i))))
  
  assign(paste0('estsfsTable.ZI69.A.mean.GCbin0', i),
         cbind(get(paste0('estsfsTable.ZI69.mean.GCbin', i)), get(paste0('polyCounts.ZI69.SI.A.mean.GCbin0', i))))
  assign(paste0('estsfsTable.ZI69.A.mel.GCbin0', i),
         cbind(get(paste0('estsfsTable.ZI69.mel.GCbin', i)), get(paste0('polyCounts.ZI69.SI.A.mel.GCbin0', i))))
  
  assign(paste0('estsfsTable.ZI21.A.mean.GCbin0', i),
         cbind(get(paste0('estsfsTable.ZI21.mean.GCbin', i)), get(paste0('polyCounts.ZI21.SI.A.mean.GCbin0', i))))
  assign(paste0('estsfsTable.ZI21.A.mel.GCbin0', i),
         cbind(get(paste0('estsfsTable.ZI21.mel.GCbin', i)), get(paste0('polyCounts.ZI21.SI.A.mel.GCbin0', i))))
  
  assign(paste0('estsfsTable.RG21.A.mean.GCbin0', i),
         cbind(get(paste0('estsfsTable.RG21.mean.GCbin', i)), get(paste0('polyCounts.RG21.SI.A.mean.GCbin0', i))))
  assign(paste0('estsfsTable.RG21.A.mel.GCbin0', i),
         cbind(get(paste0('estsfsTable.RG21.mel.GCbin', i)), get(paste0('polyCounts.RG21.SI.A.mel.GCbin0', i))))
  
  assign(paste0('estsfsTable.MD.A.diff.GCbin0', i),
         cbind(get(paste0('estsfsTable.MD.diff.GCbin', i)), get(paste0('polyCounts.MD.SI.A.diff.GCbin0', i))))
  assign(paste0('estsfsTable.ZI69.A.diff.GCbin0', i),
         cbind(get(paste0('estsfsTable.ZI69.diff.GCbin', i)), get(paste0('polyCounts.ZI69.SI.A.diff.GCbin0', i))))
}

estsfsTable.MD.A.ALL <- cbind(estsfsTable.MD.ALL, polyCounts.MD.SI.A.ALL)
estsfsTable.ZI69.A.ALL <- cbind(estsfsTable.ZI69.ALL, polyCounts.ZI69.SI.A.ALL)
estsfsTable.ZI21.A.ALL <- cbind(estsfsTable.ZI21.ALL, polyCounts.ZI21.SI.A.ALL)
estsfsTable.RG21.A.ALL <- cbind(estsfsTable.RG21.ALL, polyCounts.RG21.SI.A.ALL)

save(file = '../data/est-sfs_tables.RData', 
     list = ls(pattern = 'estsfsTable'))

# a function to take a row of one of the tables created above and
# populate the four relavant SFSs
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


# And now need to populate the SFSs using these tables
getSFSs <- function(estsfsTable){
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
  return(list('SFS_AT_2_AT' = SFS_AT_2_AT,
              'SFS_AT_2_GC' = SFS_AT_2_GC,
              'SFS_GC_2_AT' = SFS_GC_2_AT,
              'SFS_GC_2_GC' = SFS_GC_2_GC))
}


SFS_and_num_sites.MD.mean <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.MD.mean[[i]] <- getSFSs(get(paste0('estsfsTable.MD.A.mean.GCbin0', i)))
  }
)

SFS_and_num_sites.MD.sim <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.MD.sim[[i]] <- getSFSs(get(paste0('estsfsTable.MD.A.sim.GCbin0', i)))
  }
)


SFS_and_num_sites.ZI69.mean <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.ZI69.mean[[i]] <- getSFSs(get(paste0('estsfsTable.ZI69.A.mean.GCbin0', i)))
  }
)

SFS_and_num_sites.ZI69.mel <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.ZI69.mel[[i]] <- getSFSs(get(paste0('estsfsTable.ZI69.A.mel.GCbin0', i)))
  }
)


SFS_and_num_sites.ZI21.mean <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.ZI21.mean[[i]] <- getSFSs(get(paste0('estsfsTable.ZI21.A.mean.GCbin0', i)))
  }
)

SFS_and_num_sites.ZI21.mel <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.ZI21.mel[[i]] <- getSFSs(get(paste0('estsfsTable.ZI21.A.mel.GCbin0', i)))
  }
)


SFS_and_num_sites.RG21.mean <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.RG21.mean[[i]] <- getSFSs(get(paste0('estsfsTable.RG21.A.mean.GCbin0', i)))
  }
)

SFS_and_num_sites.RG21.mel <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.RG21.mel[[i]] <- getSFSs(get(paste0('estsfsTable.RG21.A.mel.GCbin0', i)))
  }
)


SFS_and_num_sites.MD.diff <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.MD.diff[[i]] <- getSFSs(get(paste0('estsfsTable.MD.A.diff.GCbin0', i)))
  }
)


SFS_and_num_sites.ZI69.diff <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.ZI69.diff[[i]] <- getSFSs(get(paste0('estsfsTable.ZI69.A.diff.GCbin0', i)))
  }
)

# sum the AT > AT and GC > GC SFSes to make a N > N SFS, for each bin
N_to_N_ALL_SFSes.MD.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.MD.mean, FUN = function(x){ x[[1]] + x[[4]] })
AT_to_GC_ALL_SFSes.MD.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.MD.mean, FUN = function(x){ x[[2]] })
GC_to_AT_ALL_SFSes.MD.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.MD.mean, FUN = function(x){ x[[3]] })

N_to_N_ALL_SFSes.MD.SI.A.sim.5GCbins <- lapply(SFS_and_num_sites.MD.sim, FUN = function(x){ x[[1]] + x[[4]] })
AT_to_GC_ALL_SFSes.MD.SI.A.sim.5GCbins <- lapply(SFS_and_num_sites.MD.sim, FUN = function(x){ x[[2]] })
GC_to_AT_ALL_SFSes.MD.SI.A.sim.5GCbins <- lapply(SFS_and_num_sites.MD.sim, FUN = function(x){ x[[3]] })


N_to_N_ALL_SFSes.ZI69.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.ZI69.mean, FUN = function(x){ x[[1]] + x[[4]] })
AT_to_GC_ALL_SFSes.ZI69.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.ZI69.mean, FUN = function(x){ x[[2]] })
GC_to_AT_ALL_SFSes.ZI69.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.ZI69.mean, FUN = function(x){ x[[3]] })

N_to_N_ALL_SFSes.ZI69.SI.A.mel.5GCbins <- lapply(SFS_and_num_sites.ZI69.mel, FUN = function(x){ x[[1]] + x[[4]] })
AT_to_GC_ALL_SFSes.ZI69.SI.A.mel.5GCbins <- lapply(SFS_and_num_sites.ZI69.mel, FUN = function(x){ x[[2]] })
GC_to_AT_ALL_SFSes.ZI69.SI.A.mel.5GCbins <- lapply(SFS_and_num_sites.ZI69.mel, FUN = function(x){ x[[3]] })


N_to_N_ALL_SFSes.ZI21.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.ZI21.mean, FUN = function(x){ x[[1]] + x[[4]] })
AT_to_GC_ALL_SFSes.ZI21.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.ZI21.mean, FUN = function(x){ x[[2]] })
GC_to_AT_ALL_SFSes.ZI21.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.ZI21.mean, FUN = function(x){ x[[3]] })

N_to_N_ALL_SFSes.ZI21.SI.A.mel.5GCbins <- lapply(SFS_and_num_sites.ZI21.mel, FUN = function(x){ x[[1]] + x[[4]] })
AT_to_GC_ALL_SFSes.ZI21.SI.A.mel.5GCbins <- lapply(SFS_and_num_sites.ZI21.mel, FUN = function(x){ x[[2]] })
GC_to_AT_ALL_SFSes.ZI21.SI.A.mel.5GCbins <- lapply(SFS_and_num_sites.ZI21.mel, FUN = function(x){ x[[3]] })


N_to_N_ALL_SFSes.RG21.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.RG21.mean, FUN = function(x){ x[[1]] + x[[4]] })
AT_to_GC_ALL_SFSes.RG21.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.RG21.mean, FUN = function(x){ x[[2]] })
GC_to_AT_ALL_SFSes.RG21.SI.A.mean.5GCbins <- lapply(SFS_and_num_sites.RG21.mean, FUN = function(x){ x[[3]] })

N_to_N_ALL_SFSes.RG21.SI.A.mel.5GCbins <- lapply(SFS_and_num_sites.RG21.mel, FUN = function(x){ x[[1]] + x[[4]] })
AT_to_GC_ALL_SFSes.RG21.SI.A.mel.5GCbins <- lapply(SFS_and_num_sites.RG21.mel, FUN = function(x){ x[[2]] })
GC_to_AT_ALL_SFSes.RG21.SI.A.mel.5GCbins <- lapply(SFS_and_num_sites.RG21.mel, FUN = function(x){ x[[3]] })

N_to_N_ALL_SFSes.MD.SI.A.diff.5GCbins <- lapply(SFS_and_num_sites.MD.diff, FUN = function(x){ x[[1]] + x[[4]] })
AT_to_GC_ALL_SFSes.MD.SI.A.diff.5GCbins <- lapply(SFS_and_num_sites.MD.diff, FUN = function(x){ x[[2]] })
GC_to_AT_ALL_SFSes.MD.SI.A.diff.5GCbins <- lapply(SFS_and_num_sites.MD.diff, FUN = function(x){ x[[3]] })

N_to_N_ALL_SFSes.ZI69.SI.A.diff.5GCbins <- lapply(SFS_and_num_sites.ZI69.diff, FUN = function(x){ x[[1]] + x[[4]] })
AT_to_GC_ALL_SFSes.ZI69.SI.A.diff.5GCbins <- lapply(SFS_and_num_sites.ZI69.diff, FUN = function(x){ x[[2]] })
GC_to_AT_ALL_SFSes.ZI69.SI.A.diff.5GCbins <- lapply(SFS_and_num_sites.ZI69.diff, FUN = function(x){ x[[3]] })

# save the SFSs, probably?
save(file = '../data/est-sfs_SFSs.RData', 
     list = ls(pattern = '_SFSes'))
