setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Write input for est-sfs for SI data, in bins and overall
# This is the format of the input file for ml-est-sfs:

library('seqinr')

# 3,0,0,17    0,0,0,1 0,0,0,1 0,0,1,0

# where the order of alleles is A, C, G, T, columns are white-space separated (!?),
# the first column is the ingroup data, and the following columns are the outgroup data.

# load the gc content and bin index stuff:
load('../data/SI_binindices_GCcontent.RData')
# load the alignments:
load('../data/SI_alignments.RData')
# load helper functions:
source('../scripts/R_helper_functions.R')

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

# and all of them:
temp.concat.alignments.all.Dsim.v202.MD.SI.A <- get_concat_seqs(alignments.intersection.MD21.SI.A,
                                                                rep(1, length(alignments.intersection.MD21.SI.A)))
temp.concat.alignments.all.Dmel.v557.ZI69.SI.A <- get_concat_seqs(alignments.intersection.ZI69.SI.A,
                                                                rep(1, length(alignments.intersection.ZI69.SI.A)))
temp.concat.alignments.all.Dmel.v557.ZI21.SI.A <- get_concat_seqs(alignments.intersection.ZI21.SI.A,
                                                                  rep(1, length(alignments.intersection.ZI21.SI.A)))
temp.concat.alignments.all.Dmel.v557.RG21.SI.A <- get_concat_seqs(alignments.intersection.RG21.SI.A,
                                                                  rep(1, length(alignments.intersection.RG21.SI.A)))

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
  # no missing data:
  indx <- apply(newMat, 2, FUN = function(cl){
    all(nlevels(as.factor(cl[polyRows])) <= 2,
        all(levels(as.factor(cl)) != 'n'))
  })
  # so this is the matrix we want to do everything to:
  newMatBiallelic <- newMat[,indx]
  return(newMatBiallelic)
}


input_file_folder_MD_SI_A_mean <- '../est-sfs/mean/MD_SI_A/input_files/'
input_file_folder_ZI69_SI_A_mean <- '../est-sfs/mean/ZI69_SI_A/input_files/'
input_file_folder_ZI21_SI_A_mean <- '../est-sfs/mean/ZI21_SI_A/input_files/'
input_file_folder_RG21_SI_A_mean <- '../est-sfs/mean/RG21_SI_A/input_files/'

input_file_folder_MD_SI_A_sim <- '../est-sfs/sim/MD_SI_A/input_files/'

input_file_folder_ZI69_SI_A_mel <- '../est-sfs/mel/ZI69_SI_A/input_files/'
input_file_folder_ZI21_SI_A_mel <- '../est-sfs/mel/ZI21_SI_A/input_files/'
input_file_folder_RG21_SI_A_mel <- '../est-sfs/mel/RG21_SI_A/input_files/'

input_file_folder_MD_SI_A_all <- '../est-sfs/all/MD_SI_A/input_files/'
input_file_folder_ZI69_SI_A_all <- '../est-sfs/all/ZI69_SI_A/input_files/'
input_file_folder_ZI21_SI_A_all <- '../est-sfs/all/ZI21_SI_A/input_files/'
input_file_folder_RG21_SI_A_all <- '../est-sfs/all/RG21_SI_A/input_files/'

# make the folders if they don't already exist
for(i in ls(pattern = 'input_file')){
  j <- get(i)
  if(!file.exists(j)){
    dir.create(j,
               recursive = T)
  }
}

### mean bins
for(i in 1:5){
  assign(paste0('estsfsMat.mean.GCbin0', i, '.MD.A'), get_mat(temp.concat.alignments.GCbin.mean.MD.SI.A[[i]]))
  lineList <- apply(get(paste0('estsfsMat.mean.GCbin0', i, '.MD.A')), 2, get_inputfile_line, 'dsim')
  fileConn <- file(paste(input_file_folder_MD_SI_A_mean,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.txt', sep = ''))
  writeLines(text = lineList,
             con = fileConn)
  close(fileConn)
}

for(i in 1:5){
  assign(paste0('estsfsMat.mean.GCbin0', i, '.ZI69.A'), get_mat(temp.concat.alignments.GCbin.mean.ZI69.SI.A[[i]]))
  lineList <- apply(get(paste0('estsfsMat.mean.GCbin0', i, '.ZI69.A')), 2, get_inputfile_line, 'dmel')
  fileConn <- file(paste(input_file_folder_ZI69_SI_A_mean,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.txt', sep = ''))
  writeLines(text = lineList,
             con = fileConn)
  close(fileConn)
}

for(i in 1:5){
  assign(paste0('estsfsMat.mean.GCbin0', i, '.ZI21.A'), get_mat(temp.concat.alignments.GCbin.mean.ZI21.SI.A[[i]]))
  lineList <- apply(get(paste0('estsfsMat.mean.GCbin0', i, '.ZI21.A')), 2, get_inputfile_line, 'dmel')
  fileConn <- file(paste(input_file_folder_ZI21_SI_A_mean,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.txt', sep = ''))
  writeLines(text = lineList,
             con = fileConn)
  close(fileConn)
}

for(i in 1:5){
  assign(paste0('estsfsMat.mean.GCbin0', i, '.RG21.A'), get_mat(temp.concat.alignments.GCbin.mean.RG21.SI.A[[i]]))
  lineList <- apply(get(paste0('estsfsMat.mean.GCbin0', i, '.RG21.A')), 2, get_inputfile_line, 'dmel')
  fileConn <- file(paste(input_file_folder_RG21_SI_A_mean,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.txt', sep = ''))
  writeLines(text = lineList,
             con = fileConn)
  close(fileConn)
}

### species-specific bins
for(i in 1:5){
  assign(paste0('estsfsMat.sim.GCbin0', i, '.MD.A'), get_mat(temp.concat.alignments.GCbin.sim.MD.SI.A[[i]]))
  lineList <- apply(get(paste0('estsfsMat.sim.GCbin0', i, '.MD.A')), 2, get_inputfile_line, 'dsim')
  fileConn <- file(paste(input_file_folder_MD_SI_A_sim,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.txt', sep = ''))
  writeLines(text = lineList,
             con = fileConn)
  close(fileConn)
}

for(i in 1:5){
  assign(paste0('estsfsMat.mel.GCbin0', i, '.ZI69.A'), get_mat(temp.concat.alignments.GCbin.mel.ZI69.SI.A[[i]]))
  lineList <- apply(get(paste0('estsfsMat.mel.GCbin0', i, '.ZI69.A')), 2, get_inputfile_line, 'dmel')
  fileConn <- file(paste(input_file_folder_ZI69_SI_A_mel,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.txt', sep = ''))
  writeLines(text = lineList,
             con = fileConn)
  close(fileConn)
}

for(i in 1:5){
  assign(paste0('estsfsMat.mel.GCbin0', i, '.ZI21.A'), get_mat(temp.concat.alignments.GCbin.mel.ZI21.SI.A[[i]]))
  lineList <- apply(get(paste0('estsfsMat.mel.GCbin0', i, '.ZI21.A')), 2, get_inputfile_line, 'dmel')
  fileConn <- file(paste(input_file_folder_ZI21_SI_A_mel,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.txt', sep = ''))
  writeLines(text = lineList,
             con = fileConn)
  close(fileConn)
}

for(i in 1:5){
  assign(paste0('estsfsMat.mel.GCbin0', i, '.RG21.A'), get_mat(temp.concat.alignments.GCbin.mel.RG21.SI.A[[i]]))
  lineList <- apply(get(paste0('estsfsMat.mel.GCbin0', i, '.RG21.A')), 2, get_inputfile_line, 'dmel')
  fileConn <- file(paste(input_file_folder_RG21_SI_A_mel,
                         '/GC_bin_',
                         sprintf("%02d", i),
                         '.txt', sep = ''))
  writeLines(text = lineList,
             con = fileConn)
  close(fileConn)
}

### all sites
estsfsMat.ALL.MD.A <- get_mat(temp.concat.alignments.all.Dsim.v202.MD.SI.A[[1]])
lineList <- apply(estsfsMat.ALL.MD.A, 2, get_inputfile_line, 'dsim')
fileConn <- file(paste(input_file_folder_MD_SI_A_all,
                       '/ALL.txt', sep = ''))
writeLines(text = lineList,
           con = fileConn)
close(fileConn)

estsfsMat.ALL.ZI69.A <- get_mat(temp.concat.alignments.all.Dmel.v557.ZI69.SI.A[[1]])
lineList <- apply(estsfsMat.ALL.ZI69.A, 2, get_inputfile_line, 'dmel')
fileConn <- file(paste(input_file_folder_ZI69_SI_A_all,
                       '/ALL.txt', sep = ''))
writeLines(text = lineList,
           con = fileConn)
close(fileConn)

estsfsMat.ALL.ZI21.A <- get_mat(temp.concat.alignments.all.Dmel.v557.ZI21.SI.A[[1]])
lineList <- apply(estsfsMat.ALL.ZI21.A, 2, get_inputfile_line, 'dmel')
fileConn <- file(paste(input_file_folder_ZI21_SI_A_all,
                       '/ALL.txt', sep = ''))
writeLines(text = lineList,
           con = fileConn)
close(fileConn)

estsfsMat.ALL.RG21.A <- get_mat(temp.concat.alignments.all.Dmel.v557.RG21.SI.A[[1]])
lineList <- apply(estsfsMat.ALL.RG21.A, 2, get_inputfile_line, 'dmel')
fileConn <- file(paste(input_file_folder_RG21_SI_A_all,
                       '/ALL.txt', sep = ''))
writeLines(text = lineList,
           con = fileConn)
close(fileConn)


# save the matrices to load later when plotting
save(file = '../data/est-sfs_matrices.RData', 
     list = ls(pattern = 'estsfsMat'))


###
test <- get_mat(temp.concat.alignments.GCbin.Dsim.v202.MD.SI.A[[1]])
get_inputfile_line(test[,1], 'dsim')
testlist <- apply(test, 2, get_inputfile_line, 'dsim')

fileConn <- file('~/Desktop/test.txt')
writeLines(text = testlist,
           con = fileConn)
close(fileConn)


#