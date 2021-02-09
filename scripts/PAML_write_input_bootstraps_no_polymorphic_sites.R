setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# A script to write input files for PAML by bootstrapping introns per bin
# this is to get CIs

library(seqinr)
library(ape)
source('../scripts/R_helper_functions.R')

# read in the sequences
load('../data/SI_alignments.RData')
# and the bin indices
load('../data/SI_binindices_GCcontent.RData')
load('../data/SI_binindices_GCcontent_DIFFERENCES.RData')

get_index <- function(mat){
  indx <- apply(mat, 2, FUN=function(cl){
    mel <- cl[1]
    sim <- cl[2]
    ZI <- cl[4:72]
    MD <- cl[73:93]
    cl_sim <- c(sim, MD)[c(sim, MD) != 'n']
    nalleles_sim <- nlevels(as.factor(cl_sim))    
    cl_mel <- c(mel, ZI)[c(mel, ZI) != 'n']
    nalleles_mel <- nlevels(as.factor(cl_mel))
    
    if(nalleles_mel > 1 || nalleles_sim > 1){
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
  return(indx)
}

get_filtered_alignment <- function(alignment){
  al <- editAlignment(alignment, 'mel|sim|yak|ZI|MD')
  mat <- makeMatrix(al)
  index <- get_index(mat)
  newmat <- mat[,index, drop = F]
  newseqs <- reconstructSequences(newmat)
  temp <- seqinr::as.alignment(nb = al$nb,
                               nam = al$nam,
                               seq = newseqs,
                               com = NULL)
  newal <- editAlignment(temp, 'mel|sim|yak')
  return(newal)
}

temp.refs.alignments.nopoly <- lapply(alignments.intersection, get_filtered_alignment)


# now concatenate the sequences according to GC bin index:
# a function to bootstrap and concatenate sequences by bin:
get_bootstrapped_concat_seqs <- function(sequences, binIndexVector, binIndex){

  sequences_temp <- sample(sequences[binIndexVector == binIndex],
                             size = length(sequences[binIndexVector == binIndex]),
                             replace = TRUE)
  # sequences_temp <- sequences_temp[!sapply(sequences_temp, is.null)] # because some of the alignments may be empty, which screws everything up
  concat_seqs <- makeGiantAlignmentCDS(sequences_temp)
}


for(i in 1:5){
  assign(paste0('seqs_concat_bootstrap.meanbins.SI.A.GCbin0', i),
         replicate(n = 1000,
                   expr = get_bootstrapped_concat_seqs(temp.refs.alignments.nopoly,
                                                       GCbinindex.5bins.mean.SI.A,
                                                       i)))
}

# test for the length of all the sequences for the SIs:
apply(seqs_concat_bootstrap.meanbins.SI.A.GCbin01, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.meanbins.SI.A.GCbin02, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.meanbins.SI.A.GCbin03, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.meanbins.SI.A.GCbin04, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.meanbins.SI.A.GCbin05, 2, FUN = function(x){
  nchar(x$seq[[1]])
})


# the order of the names (sequences) in the alignments is: mel, sim, yak:
seqs_concat_bootstrap.meanbins.SI.A.GCbin05[2,1]

for(i in 1:5){
  wd <- paste0('../PAML/no_poly/bootstraps/mean/inputfiles/GCbin0', i)
  # make the folders if they don't already exist
  if(!file.exists(wd)){
    dir.create(wd,
               recursive = T)
  }
  
  for(j in seq_len(ncol(get(paste0('seqs_concat_bootstrap.meanbins.SI.A.GCbin0', i))))){
    seqs <- as.list(get(paste0('seqs_concat_bootstrap.meanbins.SI.A.GCbin0', i)))[,j]$seq
    names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
    fileConn <- file(paste(wd,
                           '/bootstrap',
                           sprintf("%04d", j),
                           '.PHYLIP.seq', sep = ''))
    writeLines(text=c(paste(length(seqs), nchar(seqs[[1]])), # number and length of sequences
                      names(seqs)[2],
                      toupper(seqs[[2]]),
                      names(seqs)[1],
                      toupper(seqs[[1]]),
                      names(seqs)[3],
                      toupper(seqs[[3]])),
               con = fileConn)
    close(fileConn)
  }
}

###

for(i in 1:5){
  assign(paste0('seqs_concat_bootstrap.simbins.SI.A.GCbin0', i),
         replicate(n = 1000,
                   expr = get_bootstrapped_concat_seqs(temp.refs.alignments.nopoly,
                                                       GCbinindex.5bins.Dsim.SI.A,
                                                       i)))
}


# test for the length of all the sequences for the SIs:
apply(seqs_concat_bootstrap.simbins.SI.A.GCbin01, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.simbins.SI.A.GCbin02, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.simbins.SI.A.GCbin03, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.simbins.SI.A.GCbin04, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.simbins.SI.A.GCbin05, 2, FUN = function(x){
  nchar(x$seq[[1]])
})


# the order of the names (sequences) in the alignments is: mel, sim, yak:
seqs_concat_bootstrap.simbins.SI.A.GCbin01[2,1]

# now write the files
# according to the PAML manual, you can put the whole sequence on one line
# nb. want to write the variable number 'i' to two digits, to make ordering things
# easier (write '01' instead of '1', etc.)

for(i in 1:5){
  wd <- paste0('../PAML/no_poly/bootstraps/sim/inputfiles/GCbin0', i)
  # make the folders if they don't already exist
  if(!file.exists(wd)){
    dir.create(wd,
               recursive = T)
  }
  
  for(j in seq_len(ncol(get(paste0('seqs_concat_bootstrap.simbins.SI.A.GCbin0', i))))){
    seqs <- as.list(get(paste0('seqs_concat_bootstrap.simbins.SI.A.GCbin0', i)))[,j]$seq
    names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
    fileConn <- file(paste(wd,
                           '/bootstrap',
                           sprintf("%04d", j),
                           '.PHYLIP.seq', sep = ''))
    writeLines(text=c(paste(length(seqs), nchar(seqs[[1]])), # number and length of sequences
                      names(seqs)[2],
                      toupper(seqs[[2]]),
                      names(seqs)[1],
                      toupper(seqs[[1]]),
                      names(seqs)[3],
                      toupper(seqs[[3]])),
               con = fileConn)
    close(fileConn)
  }
}



for(i in 1:5){
  assign(paste0('seqs_concat_bootstrap.melbins.SI.A.GCbin0', i),
         replicate(n = 1000,
                   expr = get_bootstrapped_concat_seqs(temp.refs.alignments.nopoly,
                                                       GCbinindex.5bins.Dmel.SI.A,
                                                       i)))
}

# test for the length of all the sequences for the SIs:
apply(seqs_concat_bootstrap.melbins.SI.A.GCbin01, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.melbins.SI.A.GCbin02, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.melbins.SI.A.GCbin03, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.melbins.SI.A.GCbin04, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.melbins.SI.A.GCbin05, 2, FUN = function(x){
  nchar(x$seq[[1]])
})


# the order of the names (sequences) in the alignments is: mel, sim, yak:
seqs_concat_bootstrap.melbins.SI.A.GCbin01[2,1]

for(i in 1:5){
  wd <- paste0('../PAML/no_poly/bootstraps/mel/inputfiles/GCbin0', i)
  # make the folders if they don't already exist
  if(!file.exists(wd)){
    dir.create(wd,
               recursive = T)
  }
  
  for(j in seq_len(ncol(get(paste0('seqs_concat_bootstrap.melbins.SI.A.GCbin0', i))))){
    seqs <- as.list(get(paste0('seqs_concat_bootstrap.melbins.SI.A.GCbin0', i)))[,j]$seq
    names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
    fileConn <- file(paste(wd,
                           '/bootstrap',
                           sprintf("%04d", j),
                           '.PHYLIP.seq', sep = ''))
    writeLines(text=c(paste(length(seqs), nchar(seqs[[1]])), # number and length of sequences
                      names(seqs)[2],
                      toupper(seqs[[2]]),
                      names(seqs)[1],
                      toupper(seqs[[1]]),
                      names(seqs)[3],
                      toupper(seqs[[3]])),
               con = fileConn)
    close(fileConn)
  }
}

###

for(i in 1:5){
  assign(paste0('seqs_concat_bootstrap.diffbins.SI.A.GCbin0', i),
         replicate(n = 1000,
                   expr = get_bootstrapped_concat_seqs(temp.refs.alignments.nopoly,
                                                       GCbinindex.5bins.diff.SI.A,
                                                       i)))
}

# test for the length of all the sequences for the SIs:
apply(seqs_concat_bootstrap.diffbins.SI.A.GCbin01, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.diffbins.SI.A.GCbin02, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.diffbins.SI.A.GCbin03, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.diffbins.SI.A.GCbin04, 2, FUN = function(x){
  nchar(x$seq[[1]])
})
apply(seqs_concat_bootstrap.diffbins.SI.A.GCbin05, 2, FUN = function(x){
  nchar(x$seq[[1]])
})


# the order of the names (sequences) in the alignments is: mel, sim, yak:
seqs_concat_bootstrap.diffbins.SI.A.GCbin01[2,1]

for(i in 1:5){
  wd <- paste0('../PAML/no_poly/bootstraps/diff/inputfiles/GCbin0', i)
  # make the folders if they don't already exist
  if(!file.exists(wd)){
    dir.create(wd,
               recursive = T)
  }
  
  for(j in seq_len(ncol(get(paste0('seqs_concat_bootstrap.diffbins.SI.A.GCbin0', i))))){
    seqs <- as.list(get(paste0('seqs_concat_bootstrap.diffbins.SI.A.GCbin0', i)))[,j]$seq
    names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
    fileConn <- file(paste(wd,
                           '/bootstrap',
                           sprintf("%04d", j),
                           '.PHYLIP.seq', sep = ''))
    writeLines(text=c(paste(length(seqs), nchar(seqs[[1]])), # number and length of sequences
                      names(seqs)[2],
                      toupper(seqs[[2]]),
                      names(seqs)[1],
                      toupper(seqs[[1]]),
                      names(seqs)[3],
                      toupper(seqs[[3]])),
               con = fileConn)
    close(fileConn)
  }
}

###

seqs_concat_bootstrap.meanbins.SI.A.ALL <-
         replicate(n = 1000,
                   expr = get_bootstrapped_concat_seqs(temp.refs.alignments.nopoly,
                                                       rep(1, length(temp.refs.alignments.nopoly)),
                                                       1))


# test for the length of all the sequences for the SIs:
apply(seqs_concat_bootstrap.meanbins.SI.A.ALL, 2, FUN = function(x){
  nchar(x$seq[[1]])
})

# the order of the names (sequences) in the alignments is: mel, sim, yak:
seqs_concat_bootstrap.meanbins.SI.A.ALL[2,1]

wd <- '../PAML/no_poly/bootstraps/all/inputfiles'
# make the folders if they don't already exist
if(!file.exists(wd)){
  dir.create(wd,
             recursive = T)
}

for(j in seq_len(ncol(seqs_concat_bootstrap.meanbins.SI.A.ALL))){
  seqs <- as.list(seqs_concat_bootstrap.meanbins.SI.A.ALL)[,j]$seq
  names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
  fileConn <- file(paste(wd,
                         '/bootstrap',
                         sprintf("%04d", j),
                         '.PHYLIP.seq', sep = ''))
  writeLines(text=c(paste(length(seqs), nchar(seqs[[1]])), # number and length of sequences
                    names(seqs)[2],
                    toupper(seqs[[2]]),
                    names(seqs)[1],
                    toupper(seqs[[1]]),
                    names(seqs)[3],
                    toupper(seqs[[3]])),
             con = fileConn)
  close(fileConn)
}









#
