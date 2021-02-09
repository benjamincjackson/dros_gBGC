setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(seqinr)
library(pegas)
library(ape)
library(parallel)
as.alignment.seqinr <- seqinr::as.alignment
source('../scripts/R_helper_functions.R')

load('../data/SI_alignments.RData')
load('../data/SI_binindices_GCcontent.RData')
load('../data/SI_binindices_GCcontent_DIFFERENCES.RData')


temp.refs.alignments.intersection <- lapply(alignments.intersection, editAlignment, 'sim|mel|yak')

temp.seqs.concat.refs.alignments.intersection.Dsim.5bins <- get_concat_seqs(temp.refs.alignments.intersection,
                                                                      GCbinindex.5bins.Dsim.SI.A)

temp.seqs.concat.refs.alignments.intersection.Dmel.5bins <- get_concat_seqs(temp.refs.alignments.intersection,
                                                                            GCbinindex.5bins.Dmel.SI.A)

temp.seqs.concat.refs.alignments.intersection.mean.5bins <- get_concat_seqs(temp.refs.alignments.intersection,
                                                                            GCbinindex.5bins.mean.SI.A)

temp.seqs.concat.refs.alignments.intersection.diff.5bins <- get_concat_seqs(temp.refs.alignments.intersection,
                                                                            GCbinindex.5bins.diff.SI.A)




# test for the length of all the sequences:
sapply(temp.seqs.concat.refs.alignments.intersection.Dsim.5bins, FUN = function(x){
  nchar(x$seq[[1]])
}) # ... looks about right

# and for mel:
sapply(temp.seqs.concat.refs.alignments.intersection.Dmel.5bins, FUN = function(x){
  nchar(x$seq[[1]])
})

# and for the mean bins
sapply(temp.seqs.concat.refs.alignments.intersection.mean.5bins, FUN = function(x){
  nchar(x$seq[[1]])
})

# and for the diff bins
sapply(temp.seqs.concat.refs.alignments.intersection.diff.5bins, FUN = function(x){
  nchar(x$seq[[1]])
})


# the order of the names (sequences) in the sim alignments is mel, sim, yak:
temp.seqs.concat.refs.alignments.intersection.Dsim.5bins[[1]]$nam
temp.seqs.concat.refs.alignments.intersection.Dmel.5bins[[1]]$nam
temp.seqs.concat.refs.alignments.intersection.mean.5bins[[1]]$nam
temp.seqs.concat.refs.alignments.intersection.diff.5bins[[1]]$nam


input_file_folder_sim_SI_A <- '../PAML/sim_SI_A_GCbins/input_files_GC_bins/'
# make the folders if they don't already exist
if(!file.exists(input_file_folder_sim_SI_A)){
  dir.create(input_file_folder_sim_SI_A,
             recursive = T)
}
wd <- input_file_folder_sim_SI_A

# now write the files
# according to the PAML manual, you can put the whole sequence on one line
# nb. want to write the variable number 'i' to two digits, to make ordering things
# easier (write '01' instead of '1', etc.)
for(i in seq_along(temp.seqs.concat.refs.alignments.intersection.Dsim.5bins)){
  seqs <- as.list(temp.seqs.concat.refs.alignments.intersection.Dsim.5bins[[i]]$seq)
  names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
  fileConn <- file(paste(wd,
                         '/GC_bin_',
                         sprintf("%02d", i),
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

input_file_folder_mel_SI_A <- '../PAML/mel_SI_A_GCbins/input_files_GC_bins/'
# make the folders if they don't already exist
if(!file.exists(input_file_folder_mel_SI_A)){
  dir.create(input_file_folder_mel_SI_A,
             recursive = T)
}
wd <- input_file_folder_mel_SI_A

# now write the files
# according to the PAML manual, you can put the whole sequence on one line
# nb. want to write the variable number 'i' to two digits, to make ordering things
# easier (write '01' instead of '1', etc.)
for(i in seq_along(temp.seqs.concat.refs.alignments.intersection.Dmel.5bins)){
  seqs <- as.list(temp.seqs.concat.refs.alignments.intersection.Dmel.5bins[[i]]$seq)
  names(seqs) <- c('mel', 'sim', 'yak') # # CHECK THIS IS THE RIGHT ORDER (above)
  fileConn <- file(paste(wd,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.PHYLIP.seq', sep = ''))
  writeLines(text=c(paste(length(seqs), nchar(seqs[[1]])), # write in the order SIM MEL YAK for consistency with above
                    names(seqs)[2],
                    toupper(seqs[[2]]),
                    names(seqs)[1],
                    toupper(seqs[[1]]),
                    names(seqs)[3],
                    toupper(seqs[[3]])),
             con = fileConn)
  close(fileConn)
}

# here's also a folder where the sequences are binned by mean GC content
input_file_folder_mean_SI_A <- '../PAML/mean_SI_A_GCbins/input_files_GC_bins/'
# make the folders if they don't already exist
if(!file.exists(input_file_folder_mean_SI_A)){
  dir.create(input_file_folder_mean_SI_A,
             recursive = T)
}
wd <- input_file_folder_mean_SI_A

# now write the files
# according to the PAML manual, you can put the whole sequence on one line
# nb. want to write the variable number 'i' to two digits, to make ordering things
# easier (write '01' instead of '1', etc.)
for(i in seq_along(temp.seqs.concat.refs.alignments.intersection.mean.5bins)){
  seqs <- as.list(temp.seqs.concat.refs.alignments.intersection.mean.5bins[[i]]$seq)
  names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
  fileConn <- file(paste(wd,
                         '/GC_bin_',
                         sprintf("%02d", i),
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


# here's a folder where the sequences are binned by mean GC content
input_file_folder_diff_SI_A <- '../PAML/diff_SI_A_GCbins/input_files_GC_bins/'
# make the folders if they don't already exist
if(!file.exists(input_file_folder_diff_SI_A)){
  dir.create(input_file_folder_diff_SI_A,
             recursive = T)
}
wd <- input_file_folder_diff_SI_A

# now write the files
# according to the PAML manual, you can put the whole sequence on one line
# nb. want to write the variable number 'i' to two digits, to make ordering things
# easier (write '01' instead of '1', etc.)
for(i in seq_along(temp.seqs.concat.refs.alignments.intersection.diff.5bins)){
  seqs <- as.list(temp.seqs.concat.refs.alignments.intersection.diff.5bins[[i]]$seq)
  names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
  fileConn <- file(paste(wd,
                         '/GC_bin_',
                         sprintf("%02d", i),
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


# also write an input file for ALL the SI data, to check the overall branch lengths
temp.seqs.concat.refs.alignments.intersection.ALL <- get_concat_seqs(temp.refs.alignments.intersection,
                                                                      rep(1, length(temp.refs.alignments.intersection)))


input_file_folder_SI_A <- '../PAML/SI_A_ALL/'
# make the folders if they don't already exist
if(!file.exists(input_file_folder_SI_A)){
  dir.create(input_file_folder_SI_A,
             recursive = T)
}
wd <- input_file_folder_SI_A

# now write the files
# according to the PAML manual, you can put the whole sequence on one line
# nb. want to write the variable number 'i' to two digits, to make ordering things
# easier (write '01' instead of '1', etc.)
for(i in seq_along(temp.seqs.concat.refs.alignments.intersection.ALL)){
  seqs <- as.list(temp.seqs.concat.refs.alignments.intersection.ALL[[i]]$seq)
  names(seqs) <- c('mel', 'sim', 'yak') # CHECK THIS IS THE RIGHT ORDER (above)
  fileConn <- file(paste(wd,
                         '/ALL.PHYLIP.seq', sep = ''))
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

