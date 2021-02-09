setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(seqinr)
library(pegas)
library(ape)
library(parallel)
source('../scripts/R_helper_functions.R')

load('../data/SI_alignments.RData')
load('../data/SI_binindices_GCcontent.RData')
load('../data/SI_binindices_GCcontent_DIFFERENCES.RData')


temp.seqs.noRG.alignments.intersection <- lapply(alignments.intersection, editAlignment, 'mel|sim|yak|ZI|MD')

temp.seqs.concat.alignments.intersection.Dsim.5bins <- get_concat_seqs(temp.seqs.noRG.alignments.intersection,
                                                                      GCbinindex.5bins.Dsim.SI.A)

temp.seqs.concat.alignments.intersection.Dmel.5bins <- get_concat_seqs(temp.seqs.noRG.alignments.intersection,
                                                                            GCbinindex.5bins.Dmel.SI.A)

temp.seqs.concat.alignments.intersection.mean.5bins <- get_concat_seqs(temp.seqs.noRG.alignments.intersection,
                                                                            GCbinindex.5bins.mean.SI.A)

temp.seqs.concat.alignments.intersection.diff.5bins <- get_concat_seqs(temp.seqs.noRG.alignments.intersection,
                                                                            GCbinindex.5bins.diff.SI.A)


# test for the length of all the sequences:
sapply(temp.seqs.concat.alignments.intersection.Dsim.5bins, FUN = function(x){
  nchar(x$seq[[1]])
}) # ... looks about right

# and for mel:
sapply(temp.seqs.concat.alignments.intersection.Dmel.5bins, FUN = function(x){
  nchar(x$seq[[1]])
})

# and for the mean bins
sapply(temp.seqs.concat.alignments.intersection.mean.5bins, FUN = function(x){
  nchar(x$seq[[1]])
})

sapply(temp.seqs.concat.alignments.intersection.diff.5bins, FUN = function(x){
  nchar(x$seq[[1]])
})


# the order of the names (sequences) in the sim alignments is mel, sim, yak, ZI, MD:
temp.seqs.concat.alignments.intersection.Dsim.5bins[[1]]$nam
temp.seqs.concat.alignments.intersection.Dmel.5bins[[1]]$nam
temp.seqs.concat.alignments.intersection.mean.5bins[[1]]$nam
temp.seqs.concat.alignments.intersection.diff.5bins[[1]]$nam


# now I need to filter out sites that are not monomorphic in either polymorphism
# sample

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

seqs.mean.refs <- lapply(temp.seqs.concat.alignments.intersection.mean.5bins, FUN = function(al){
  mat <- makeMatrix(al)
  index <- get_index(mat)
  newmat <- mat[,index]
  newseqs <- reconstructSequences(newmat)
  temp <- seqinr::as.alignment(nb = al$nb,
                                  nam = al$nam,
                                  seq = newseqs,
                                  com = NULL)
  newal <- editAlignment(temp, 'mel|sim|yak')
  return(newal)
})

seqs.diff.refs <- lapply(temp.seqs.concat.alignments.intersection.diff.5bins, FUN = function(al){
  mat <- makeMatrix(al)
  index <- get_index(mat)
  newmat <- mat[,index]
  newseqs <- reconstructSequences(newmat)
  temp <- seqinr::as.alignment(nb = al$nb,
                               nam = al$nam,
                               seq = newseqs,
                               com = NULL)
  newal <- editAlignment(temp, 'mel|sim|yak')
  return(newal)
})

seqs.dsim.refs <- lapply(temp.seqs.concat.alignments.intersection.Dsim.5bins, FUN = function(al){
  mat <- makeMatrix(al)
  index <- get_index(mat)
  newmat <- mat[,index]
  newseqs <- reconstructSequences(newmat)
  temp <- seqinr::as.alignment(nb = al$nb,
                               nam = al$nam,
                               seq = newseqs,
                               com = NULL)
  newal <- editAlignment(temp, 'mel|sim|yak')
  return(newal)
})
  
seqs.dmel.refs <- lapply(temp.seqs.concat.alignments.intersection.Dmel.5bins, FUN = function(al){
  mat <- makeMatrix(al)
  index <- get_index(mat)
  newmat <- mat[,index]
  newseqs <- reconstructSequences(newmat)
  temp <- seqinr::as.alignment(nb = al$nb,
                               nam = al$nam,
                               seq = newseqs,
                               com = NULL)
  newal <- editAlignment(temp, 'mel|sim|yak')
  return(newal)
})

# the order of the names (sequences) in the alignments is mel, sim, yak:
seqs.dmel.refs[[1]]$nam
seqs.dsim.refs[[1]]$nam
seqs.mean.refs[[1]]$nam
seqs.diff.refs[[1]]$nam


# checking the length of things:
sapply(temp.seqs.concat.alignments.intersection.mean.5bins, FUN = function(x){
  nchar(x$seq[[1]])
})
sapply(seqs.mean.refs, FUN = function(x){
  nchar(x$seq[[1]])
})

sapply(temp.seqs.concat.alignments.intersection.diff.5bins, FUN = function(x){
  nchar(x$seq[[1]])
})
sapply(seqs.diff.refs, FUN = function(x){
  nchar(x$seq[[1]])
})

sapply(temp.seqs.concat.alignments.intersection.Dmel.5bins, FUN = function(x){
  nchar(x$seq[[1]])
})
sapply(seqs.dmel.refs, FUN = function(x){
  nchar(x$seq[[1]])
})

sapply(temp.seqs.concat.alignments.intersection.Dsim.5bins, FUN = function(x){
  nchar(x$seq[[1]])
})
sapply(seqs.dsim.refs, FUN = function(x){
  nchar(x$seq[[1]])
})


input_file_folder_sim_SI_A <- '../PAML/no_poly/sim_SI_A_GCbins/input_files_GC_bins/'
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
for(i in seq_along(seqs.dsim.refs)){
  seqs <- as.list(seqs.dsim.refs[[i]]$seq)
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

input_file_folder_mel_SI_A <- '../PAML/no_poly/mel_SI_A_GCbins/input_files_GC_bins/'
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
for(i in seq_along(seqs.dmel.refs)){
  seqs <- as.list(seqs.dmel.refs[[i]]$seq)
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

# here's also a fodler where the sequences are binned by mean GC content
input_file_folder_mean_SI_A <- '../PAML/no_poly/mean_SI_A_GCbins/input_files_GC_bins/'
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
for(i in seq_along(seqs.mean.refs)){
  seqs <- as.list(seqs.mean.refs[[i]]$seq)
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


# here's a folder where the sequences are binned by diff GC content
input_file_folder_diff_SI_A <- '../PAML/no_poly/diff_SI_A_GCbins/input_files_GC_bins/'
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
for(i in seq_along(seqs.diff.refs)){
  seqs <- as.list(seqs.diff.refs[[i]]$seq)
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
temp.seqs.concat.all <- get_concat_seqs(seqs.mean.refs, rep(1, length(5)))


input_file_folder_SI_A <- '../PAML/no_poly/SI_A_ALL/'
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
for(i in seq_along(temp.seqs.concat.all)){
  seqs <- as.list(temp.seqs.concat.all[[i]]$seq)
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

