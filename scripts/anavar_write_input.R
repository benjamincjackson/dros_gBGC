# write input files to run glemin's models using Kai's implementation (anavar)
# THIS IS JUST A TEST AT THE MOMENT - RUN THE MOST EXTREME diff BINS

# load the matrices that we used to write the est-sfs input files with:
load('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/data/est-sfs_matrices_DIFFS.RData')
load('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/data/est-sfs_matrices.RData')

# load the bin indices and GC contents
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_binindices_GCcontent_DIFFERENCES.RData')
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_binindices_GCcontent.RData')

# read in the persite output files from est-sfs
for(i in 1:5){
  assign(paste0('estsfsTable.MD.diff.GCbin', i),
         read.table(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/est-sfs/diff/MD_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}

for(i in 1:5){
  assign(paste0('estsfsTable.ZI69.diff.GCbin', i),
         read.table(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/est-sfs/diff/ZI69_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}

for(i in 1:5){
  assign(paste0('estsfsTable.MD.mean.GCbin', i),
         read.table(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/est-sfs/mean/MD_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}

for(i in 1:5){
  assign(paste0('estsfsTable.ZI69.mean.GCbin', i),
         read.table(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/est-sfs/mean/ZI69_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}

for(i in 1:5){
  assign(paste0('estsfsTable.MD.sim.GCbin', i),
         read.table(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/est-sfs/sim/MD_SI_A/output_files/',
                           'GC_bin_0', i, '.pvalues'),
                    skip = 8))
}

for(i in 1:5){
  assign(paste0('estsfsTable.ZI69.mel.GCbin', i),
         read.table(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/est-sfs/mel/ZI69_SI_A/output_files/',
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
  assign(paste0('polyCounts.MD.SI.A.diff.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.diff.GCbin0', i, '.MD.A'))))
}

for(i in 1:5){
  assign(paste0('polyCounts.ZI69.SI.A.diff.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.diff.GCbin0', i, '.ZI69.A'))))
}


for(i in 1:5){
  assign(paste0('polyCounts.MD.SI.A.mean.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.mean.GCbin0', i, '.MD.A'))))
}

for(i in 1:5){
  assign(paste0('polyCounts.ZI69.SI.A.mean.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.mean.GCbin0', i, '.ZI69.A'))))
}


for(i in 1:5){
  assign(paste0('polyCounts.MD.SI.A.sim.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.sim.GCbin0', i, '.MD.A'))))
}

for(i in 1:5){
  assign(paste0('polyCounts.ZI69.SI.A.mel.GCbin0', i),
         getPolyCounts(get(paste0('estsfsMat.mel.GCbin0', i, '.ZI69.A'))))
}

# stick these tables together with the probability tables, which we can
# work with to calculate DAF

for(i in 1:5){
  assign(paste0('estsfsTable.MD.A.diff.GCbin0', i),
         cbind(get(paste0('estsfsTable.MD.diff.GCbin', i)), get(paste0('polyCounts.MD.SI.A.diff.GCbin0', i))))
  
  assign(paste0('estsfsTable.ZI69.A.diff.GCbin0', i),
         cbind(get(paste0('estsfsTable.ZI69.diff.GCbin', i)), get(paste0('polyCounts.ZI69.SI.A.diff.GCbin0', i))))
}

for(i in 1:5){
  assign(paste0('estsfsTable.MD.A.mean.GCbin0', i),
         cbind(get(paste0('estsfsTable.MD.mean.GCbin', i)), get(paste0('polyCounts.MD.SI.A.mean.GCbin0', i))))
  
  assign(paste0('estsfsTable.ZI69.A.mean.GCbin0', i),
         cbind(get(paste0('estsfsTable.ZI69.mean.GCbin', i)), get(paste0('polyCounts.ZI69.SI.A.mean.GCbin0', i))))
}

for(i in 1:5){
  assign(paste0('estsfsTable.MD.A.sim.GCbin0', i),
         cbind(get(paste0('estsfsTable.MD.sim.GCbin', i)), get(paste0('polyCounts.MD.SI.A.sim.GCbin0', i))))
  
  assign(paste0('estsfsTable.ZI69.A.mel.GCbin0', i),
         cbind(get(paste0('estsfsTable.ZI69.mel.GCbin', i)), get(paste0('polyCounts.ZI69.SI.A.mel.GCbin0', i))))
}

# a function to take a row of one of the tables created above and
# populate the four relevant SFSs when the site is fixed
getFixedCounts <- function(rw, n){
  newrw <- unname(unlist(rw))

  AT_count <- 0L
  GC_count <- 0L
  
  # counts of alleles in the poly data:
  cA <- newrw[8]
  cC <-  newrw[9]
  cG <-  newrw[10]
  cT <-  newrw[11]
  
  countVector <- c('A' = cA, 'C' = cC, 'G' = cG, 'T' = cT)
  
  # which allele is fixed?
  fixedAl <- names(countVector)[which.max(countVector)]
  
  # populate the SFSs:
  if(fixedAl == 'A'){
    AT_count <- AT_count + 1
  }

  if(fixedAl == 'C'){
    GC_count <- GC_count + 1
  }

  if(fixedAl == 'G'){
    GC_count <- GC_count + 1
  }
  
  if(fixedAl == 'T'){
    AT_count <- AT_count + 1
  }
  
  return(list('AT_count' = AT_count,
              'GC_count' = GC_count))
}

# a function to take a row of one of the tables created above and
# populate the four relevant SFSs when the site is segregating
popSFSs_poly <- function(rw, n){
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
    # the sfs column to populate if A is ancestral:
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


# And now need to populate the SFSs using these tables and the two functions
# above:
getSFSs_andFixedCounts <- function(estsfsTable){
  n <- sum(estsfsTable[1,8:11])
  
  SFS_AT_2_AT <- vector('numeric', length = n + 1)
  SFS_AT_2_GC <- vector('numeric', length = n + 1)
  SFS_GC_2_GC <- vector('numeric', length = n + 1)
  SFS_GC_2_AT <- vector('numeric', length = n + 1)
  
  GC_count <- 0L
  AT_count <- 0L
  
  for(i in seq_len(nrow(estsfsTable))){
    results <- NULL
    results_fixed <- NULL
    rw <- estsfsTable[i,]
    newrw <- unname(unlist(rw))
    countVector <- c(newrw[8:11])
    
    # if this site isn't segregating, then apply the appropriate function:
    if(identical(sum(countVector != 0), 1L)){
      results_fixed <- getFixedCounts(rw, n)
      
      GC_count <- GC_count + results_fixed$GC_count
      AT_count <- AT_count + results_fixed$AT_count
    }
    # if this site is segregating, then apply the other function:
    else{
      results <- popSFSs_poly(rw, n)
      
      SFS_AT_2_AT <- SFS_AT_2_AT + results$SFS_AT_2_AT
      SFS_AT_2_GC <- SFS_AT_2_GC + results$SFS_AT_2_GC
      SFS_GC_2_AT <- SFS_GC_2_AT + results$SFS_GC_2_AT
      SFS_GC_2_GC <- SFS_GC_2_GC + results$SFS_GC_2_GC
    }
  }
  
  return(list('SFS_AT_2_AT' = SFS_AT_2_AT,
              'SFS_AT_2_GC' = SFS_AT_2_GC,
              'SFS_GC_2_AT' = SFS_GC_2_AT,
              'SFS_GC_2_GC' = SFS_GC_2_GC,
              'GC_fixed' = GC_count,
              'AT_fixed' = AT_count))
}

SFS_and_num_sites.MD.diff <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.MD.diff[[i]] <- getSFSs_andFixedCounts(get(paste0('estsfsTable.MD.A.diff.GCbin0', i)))
  }
)

SFS_and_num_sites.ZI69.diff <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.ZI69.diff[[i]] <- getSFSs_andFixedCounts(get(paste0('estsfsTable.ZI69.A.diff.GCbin0', i)))
  }
)


SFS_and_num_sites.MD.mean <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.MD.mean[[i]] <- getSFSs_andFixedCounts(get(paste0('estsfsTable.MD.A.mean.GCbin0', i)))
  }
)

SFS_and_num_sites.ZI69.mean <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.ZI69.mean[[i]] <- getSFSs_andFixedCounts(get(paste0('estsfsTable.ZI69.A.mean.GCbin0', i)))
  }
)


SFS_and_num_sites.MD.species <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.MD.species[[i]] <- getSFSs_andFixedCounts(get(paste0('estsfsTable.MD.A.sim.GCbin0', i)))
  }
)

SFS_and_num_sites.ZI69.species <- vector('list', 5)
system.time(
  for(i in 1:5){
    SFS_and_num_sites.ZI69.species[[i]] <- getSFSs_andFixedCounts(get(paste0('estsfsTable.ZI69.A.mel.GCbin0', i)))
  }
)


# save the SFSs, probably?
save(file = '/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/data/anavar_SFSs.RData', 
     list = ls(pattern = '_SFSes'))


# a function to write an anavar input file given some things:
# 1) an object returned by getSFSs_andFixedCounts
# 2) a file path (and file name?)
# 3) optimisation algorithm to be used by the file
# 4) the model subset to use (M1star, M1, M0star, M0 etc.)

writeAnavarCtrlFile <- function(fileNameAndPath,
                                SFSinfo,
                                algorithmToUse = c('LBFGS', 'SLSQP'),
                                modelSubset = c('none', 'M0',  'M1', 'M0*', 'M1*')){
  
  fileConn <- file(fileNameAndPath)
  
  modelSubset <- match.arg(modelSubset)
  algorithmToUse <- match.arg(algorithmToUse)
  
  SFSneu <- SFSinfo$SFS_AT_2_AT + SFSinfo$SFS_GC_2_GC
  L_all <- sum(unlist(SFSinfo[1:4])) + SFSinfo$GC_fixed + SFSinfo$AT_fixed
  L_AT <- sum(unlist(SFSinfo$SFS_AT_2_AT)) + sum(unlist(SFSinfo$SFS_AT_2_GC)) + SFSinfo$AT_fixed
  L_GC <- sum(unlist(SFSinfo$SFS_GC_2_GC)) + sum(unlist(SFSinfo$SFS_GC_2_AT)) + SFSinfo$GC_fixed
  n <- length(SFSinfo$SFS_AT_2_AT) - 1
  
  writeLines(text = c(
    '[algorithm_commands]',
    paste0('search_algorithm: NLOPT_LD_', algorithmToUse, collapse = ''),
    'maxeval: 100000',
    'maxtime: 600',
    'num_searches: 1000',
    'nnoimp: 3',
    'maximp: 10',
    'optional: true',
    'size: 5000',
    'key: 5',
    'epsabs: 1e-40',
    'epsrel: 1e-8',
    'rftol: 1e-9',
    'init: random',
    '',
    '[model_commands]',
    'model: gBGC_EXTENDED_M1*',
    '',
    paste0('n: ', n, collapse = ''),
    'r_range: 0.01, 10',
    '',
    '[neutral_SNPs]',
    paste0('m: ', round(L_all, 2), collapse = ''),
    paste0('sfs: ', paste0(round(SFSneu[2:n], 2), collapse = ', ')),
    'theta_range: 1e-10, 0.1',
    'gamma_range: -100, 50',
    'e_range: 0, 0.5',
    '\n',
    '[ws_SNPs]',
    paste0('m: ', round(L_AT, 2), collapse = ''),
    paste0('sfs: ', paste0(round(SFSinfo$SFS_AT_2_GC[2:n], 2), collapse = ', ')),
    'theta_range: 1e-10, 0.1',
    'gamma_range: -100, 50',
    'e_range: 0, 0.5',
    '',
    '[sw_SNPs]',
    paste0('m: ', round(L_GC, 2), collapse = ''),
    paste0('sfs: ', paste0(round(SFSinfo$SFS_GC_2_AT[2:n], 2), collapse = ', ')),
    'theta_range: 1e-10, 0.1',
    'gamma_range: -100, 50',
    'e_range: 0, 0.5',
    '',
    paste0('constraint: ', modelSubset, collapse = '')
  ),
  con = fileConn)
  close(fileConn)
}

# TEST
# writeAnavarCtrlFile(fileNameAndPath = '~/Desktop/test.anavar',
#                     SFSinfo = SFS_and_num_sites.MD.diff[[1]],
#                     algorithmToUse = 'LBFGS')


commandsFile <- '/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/anavar/commands.txt'
# And then let's make folders and populate then with control files?
for(bintype in c('mean', 'diff', 'species')){
  for(binnumber in 1:5){
    for(pop in c('MD', 'ZI69')){
      for(model in c('none', 'M0',  'M1', 'M0*', 'M1*')){
          for(algorithm in c('LBFGS', 'SLSQP')){
          
          if(model == 'M0*') {
            wd <- paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/anavar', '/',
                         bintype, '/',
                         pop, '/',
                         binnumber, '/',
                         'M0star/',
                         collapse = '')
            
          } else if(model == 'M1*'){
            wd <- paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/anavar', '/',
                         bintype, '/',
                         pop, '/',
                         binnumber, '/',
                         'M1star/',
                         collapse = '')       
          } else {
            wd <- paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/anavar', '/',
                         bintype, '/',
                         pop, '/',
                         binnumber, '/',
                         model, '/',
                         collapse = '')
          }

          if(!file.exists(wd)){
            dir.create(wd,
                       recursive = T)
          }
          
          writeAnavarCtrlFile(fileNameAndPath = paste0(wd, '/',
                                                       paste0('gbgc_extended_m1s_',
                                                              algorithm,
                                                              '.txt',
                                                              collapse = ''),
                                                       collapse = ''),
                              SFSinfo = get(paste0('SFS_and_num_sites.', pop, '.', bintype, collapse = ''))[[binnumber]],
                              algorithmToUse = algorithm,
                              modelSubset = model)
          
          line1 <- paste0('/Users/ben/programs/anavar1.4/anavar_src/anavar ',
                          wd, 'gbgc_extended_m1s_LBFGS.txt ',
                          wd, 'results_LBFGS.txt ',
                          wd, 'log.txt $RANDOM; rm ',
                          wd, 'log.txt', collapse = '')
          line2 <- paste0('/Users/ben/programs/anavar1.4/anavar_src/anavar ',
                          wd, 'gbgc_extended_m1s_SLSQP.txt ',
                          wd, 'results_SLSQP.txt ',
                          wd, 'log.txt $RANDOM; rm ',
                          wd, 'log.txt', collapse = '')
          
          write(x = line1, file = commandsFile, append = TRUE)
          write(x = line2, file = commandsFile, append = TRUE)
          
        }
      }
    }
  }
}



