setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# write est-sfs input for GC-difference bins - we want to know if gamma is
# particularly strong in regions that have divergence heavily in GC content
# between mel and sim

library(seqinr)
library(pegas)
library(ape)
library(parallel)
source('../scripts/R_helper_functions.R')

# load the bin index and GC content data
load('../data/SI_binindices_GCcontent.RData')
load('../data/SI_binindices_GCcontent_DIFFERENCES.RData')

# load the alignments
load('../data/SI_alignments.RData')


temp.concat.alignments.GCbin.diff.MD.SI.A <- get_concat_seqs(alignments.intersection.MD21.SI.A,
                                                             GCbinindex.5bins.diff.SI.A)

temp.concat.alignments.GCbin.diff.ZI69.SI.A <- get_concat_seqs(alignments.intersection.ZI69.SI.A,
                                                               GCbinindex.5bins.diff.SI.A)

# a look up table to add counts to the relavent column for the est-sfs input format

# function to return a line of est-sfs input file from a column in a matrix of nucleotides
# apply this function over a matrix created from a fasta-format alignment
get_inputfile_line <- function(cl, ingroup = c('dsim', 'dmel')){
  
  LUTindex <- c('A' = 1L, 'C' = 2L, 'G' = 3L, 'T' = 4L)
  
  simAl <- toupper(cl[1])
  melAl <- toupper(cl[2])
  yakAl <- toupper(cl[3])
  polyAls <- toupper(cl[4:length(cl)])
  
  outcol1 <- c(0,0,0,0)
  outcol2 <- c(0,0,0,0)
  
  polycol <- c(0,0,0,0)
  
  if(identical(ingroup, 'dsim')) {
    outcol1[LUTindex[melAl]] <- 1
  }
  
  else {
    outcol1[LUTindex[simAl]] <- 1
  }
  
  outcol2[LUTindex[yakAl]] <- 1
  
  for(i in seq_along(polyAls)) {
    al <- polyAls[i]
    polycol[LUTindex[al]] <- polycol[LUTindex[al]] + 1
  }
  
  inputstring <- paste(paste0(polycol, collapse = ','), paste0(outcol1, collapse = ','), paste0(outcol2, collapse = ','), sep = '\t')
  return(inputstring)
}

# function to return a matrix of  sites, with some filters:
# max two alleles in the polymorphism dataset
# no missing data
get_mat <- function(alignment){
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
  # so this is the matrix we want to do everything to:
  newMatBiallelic <- newMat[,indx]
  return(newMatBiallelic)
}

input_file_folder_MD_SI_A_diff <- '../est-sfs/diff/MD_SI_A/input_files/'
input_file_folder_ZI69_SI_A_diff <- '../est-sfs/diff/ZI69_SI_A/input_files/'

# make the folders if they don't already exist
for(i in ls(pattern = 'input_file')){
  j <- get(i)
  if(!file.exists(j)){
    dir.create(j,
               recursive = T)
  }
}

for(i in 1:5){
  assign(paste0('estsfsMat.diff.GCbin0', i, '.MD.A'), get_mat(temp.concat.alignments.GCbin.diff.MD.SI.A[[i]]))
  lineList <- apply(get(paste0('estsfsMat.diff.GCbin0', i, '.MD.A')), 2, get_inputfile_line, 'dsim')
  fileConn <- file(paste(input_file_folder_MD_SI_A_diff,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.txt', sep = ''))
  writeLines(text = lineList,
             con = fileConn)
  close(fileConn)
}

for(i in 1:5){
  assign(paste0('estsfsMat.diff.GCbin0', i, '.ZI69.A'), get_mat(temp.concat.alignments.GCbin.diff.ZI69.SI.A[[i]]))
  lineList <- apply(get(paste0('estsfsMat.diff.GCbin0', i, '.ZI69.A')), 2, get_inputfile_line, 'dmel')
  fileConn <- file(paste(input_file_folder_ZI69_SI_A_diff,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.txt', sep = ''))
  writeLines(text = lineList,
             con = fileConn)
  close(fileConn)
}

# save the matrices to load later when plotting
save(file = '../data/est-sfs_matrices_DIFFS.RData', 
     list = ls(pattern = 'estsfsMat'))












#