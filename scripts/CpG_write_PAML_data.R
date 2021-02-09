setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# write PAML input from the (non)CpG/prone -split data

library(seqinr)
library(pegas)
library(ape)
library(parallel)
source('R_helper_functions.R')

as.alignment.seqinr <- seqinr::as.alignment

files.CpG <- grep(list.files('../CpG/fasta/CpG/', full.names = T),
                           pattern = 'X',
                           invert = T,
                           value = T)

files.nonCpG <- grep(list.files('../CpG/fasta/nonCpG/', full.names = T),
                  pattern = 'X',
                  invert = T,
                  value = T)

files.CpGprone <- grep(list.files('../CpG/fasta/CpGprone/', full.names = T),
                  pattern = 'X',
                  invert = T,
                  value = T)

files.nonCpGprone <- grep(list.files('../CpG/fasta/nonCpGprone/', full.names = T),
                     pattern = 'X',
                     invert = T,
                     value = T)

alignments.CpG <- lapply(files.CpG,
                         read.alignment,
                         format = 'fasta')

alignments.nonCpG <- lapply(files.nonCpG,
                            read.alignment,
                            format = 'fasta')

alignments.CpGprone <- lapply(files.CpGprone,
                              read.alignment,
                              format = 'fasta')

alignments.nonCpGprone <- lapply(files.nonCpGprone,
                                 read.alignment,
                                 format = 'fasta')

######

alignments.CpG.refs <- lapply(alignments.CpG,
                              editAlignment,
                              'sim|mel|yak')
alignments.CpG.refs.concat <- get_concat_seqs(alignments.CpG.refs,
                                            rep(1, length(alignments.CpG.refs)))

alignments.nonCpG.refs <- lapply(alignments.nonCpG,
                              editAlignment,
                              'sim|mel|yak')
alignments.nonCpG.refs.concat <- get_concat_seqs(alignments.nonCpG.refs,
                                                rep(1, length(alignments.nonCpG.refs)))

alignments.CpGprone.refs <- lapply(alignments.CpGprone,
                              editAlignment,
                              'sim|mel|yak')
alignments.CpGprone.refs.concat <- get_concat_seqs(alignments.CpGprone.refs,
                                                   rep(1, length(alignments.CpGprone.refs)))

alignments.nonCpGprone.refs <- lapply(alignments.nonCpGprone,
                                 editAlignment,
                                 'sim|mel|yak')
alignments.nonCpGprone.refs.concat <- get_concat_seqs(alignments.nonCpGprone.refs,
                                                     rep(1, length(alignments.nonCpGprone.refs)))

# the order of the names (sequences) in the sim alignments is mel, sim, yak:
alignments.CpG.refs.concat[[1]]$nam
alignments.nonCpG.refs.concat[[1]]$nam
alignments.CpGprone.refs.concat[[1]]$nam
alignments.nonCpGprone.refs.concat[[1]]$nam


## write PAML input

input_file_folder_CpG <- '../CpG/PAML/CpG/pointestimate'
# make the folders if they don't already exist
if(!file.exists(input_file_folder_CpG)){
  dir.create(input_file_folder_CpG,
             recursive = T)
}
wd <- input_file_folder_CpG

for(i in seq_along(alignments.CpG.refs.concat)){
  seqs <- as.list(alignments.CpG.refs.concat[[i]]$seq)
  names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
  fileConn <- file(paste(wd,
                         '/CpG.PHYLIP.seq', sep = ''))
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

input_file_folder_nonCpG <- '../CpG/PAML/nonCpG/pointestimate'
# make the folders if they don't already exist
if(!file.exists(input_file_folder_nonCpG)){
  dir.create(input_file_folder_nonCpG,
             recursive = T)
}
wd <- input_file_folder_nonCpG

for(i in seq_along(alignments.nonCpG.refs.concat)){
  seqs <- as.list(alignments.nonCpG.refs.concat[[i]]$seq)
  names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
  fileConn <- file(paste(wd,
                         '/nonCpG.PHYLIP.seq', sep = ''))
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

input_file_folder_CpGprone <- '../CpG/PAML/CpGprone/pointestimate'
# make the folders if they don't already exist
if(!file.exists(input_file_folder_CpGprone)){
  dir.create(input_file_folder_CpGprone,
             recursive = T)
}
wd <- input_file_folder_CpGprone

for(i in seq_along(alignments.CpGprone.refs.concat)){
  seqs <- as.list(alignments.CpGprone.refs.concat[[i]]$seq)
  names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
  fileConn <- file(paste(wd,
                         '/CpGprone.PHYLIP.seq', sep = ''))
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

input_file_folder_nonCpGprone <- '../CpG/PAML/nonCpGprone/pointestimate'
# make the folders if they don't already exist
if(!file.exists(input_file_folder_nonCpGprone)){
  dir.create(input_file_folder_nonCpGprone,
             recursive = T)
}
wd <- input_file_folder_nonCpGprone

for(i in seq_along(alignments.nonCpGprone.refs.concat)){
  seqs <- as.list(alignments.nonCpGprone.refs.concat[[i]]$seq)
  names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
  fileConn <- file(paste(wd,
                         '/nonCpGprone.PHYLIP.seq', sep = ''))
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


#### write the bootstraps for CIs

# a function to bootstrap and concatenate sequences by bin:
get_bootstrapped_concat_seqs <- function(sequences, binIndexVector, binIndex){
  
  sequences_temp <- sample(sequences[binIndexVector == binIndex],
                           size = length(sequences[binIndexVector == binIndex]),
                           replace = TRUE)
  # sequences_temp <- sequences_temp[!sapply(sequences_temp, is.null)] # because some of the alignments may be empty, which screws everything up
  concat_seqs <- makeGiantAlignmentCDS(sequences_temp)
}


alignments.CpG.refs.concat.bootstrap <-
  replicate(n = 1000,
            expr = get_bootstrapped_concat_seqs(alignments.CpG.refs,
                                                rep(1, length(alignments.CpG.refs.concat)),
                                                1))

alignments.nonCpG.refs.concat.bootstrap <-
  replicate(n = 1000,
            expr = get_bootstrapped_concat_seqs(alignments.nonCpG.refs,
                                                rep(1, length(alignments.nonCpG.refs.concat)),
                                                1))

alignments.CpGprone.refs.concat.bootstrap <-
  replicate(n = 1000,
            expr = get_bootstrapped_concat_seqs(alignments.CpGprone.refs,
                                                rep(1, length(alignments.CpGprone.refs.concat)),
                                                1))

alignments.nonCpGprone.refs.concat.bootstrap <-
  replicate(n = 1000,
            expr = get_bootstrapped_concat_seqs(alignments.nonCpGprone.refs,
                                                rep(1, length(alignments.nonCpGprone.refs.concat)),
                                                1))


# test for the length of all the boostrapped samples:
nchar(alignments.CpG.refs.concat[[1]]$seq[[1]])
apply(alignments.CpG.refs.concat.bootstrap, 2, FUN = function(x){
  nchar(x$seq[[1]])
})

nchar(alignments.nonCpG.refs.concat[[1]]$seq[[1]])
apply(alignments.nonCpG.refs.concat.bootstrap, 2, FUN = function(x){
  nchar(x$seq[[1]])
})

nchar(alignments.CpG.refs.concat[[1]]$seq[[1]]) + nchar(alignments.nonCpG.refs.concat[[1]]$seq[[1]])


nchar(alignments.CpGprone.refs.concat[[1]]$seq[[1]])
apply(alignments.CpGprone.refs.concat.bootstrap, 2, FUN = function(x){
  nchar(x$seq[[1]])
})

nchar(alignments.nonCpGprone.refs.concat[[1]]$seq[[1]])
apply(alignments.nonCpGprone.refs.concat.bootstrap, 2, FUN = function(x){
  nchar(x$seq[[1]])
})

nchar(alignments.CpGprone.refs.concat[[1]]$seq[[1]]) + nchar(alignments.nonCpGprone.refs.concat[[1]]$seq[[1]])


# write the bootstraps

wd <- '../CpG/PAML/CpG/bootstraps/inputfiles'
# make the folders if they don't already exist
if(!file.exists(wd)){
  dir.create(wd,
             recursive = T)
}

for(j in seq_len(ncol(alignments.CpG.refs.concat.bootstrap))){
  seqs <- as.list(alignments.CpG.refs.concat.bootstrap)[,j]$seq
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

wd <- '../CpG/PAML/nonCpG/bootstraps/inputfiles'
# make the folders if they don't already exist
if(!file.exists(wd)){
  dir.create(wd,
             recursive = T)
}

for(j in seq_len(ncol(alignments.nonCpG.refs.concat.bootstrap))){
  seqs <- as.list(alignments.nonCpG.refs.concat.bootstrap)[,j]$seq
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

wd <- '../CpG/PAML/CpGprone/bootstraps/inputfiles'
# make the folders if they don't already exist
if(!file.exists(wd)){
  dir.create(wd,
             recursive = T)
}

for(j in seq_len(ncol(alignments.CpGprone.refs.concat.bootstrap))){
  seqs <- as.list(alignments.CpGprone.refs.concat.bootstrap)[,j]$seq
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

wd <- '../CpG/PAML/nonCpGprone/bootstraps/inputfiles'
# make the folders if they don't already exist
if(!file.exists(wd)){
  dir.create(wd,
             recursive = T)
}

for(j in seq_len(ncol(alignments.nonCpGprone.refs.concat.bootstrap))){
  seqs <- as.list(alignments.nonCpGprone.refs.concat.bootstrap)[,j]$seq
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

