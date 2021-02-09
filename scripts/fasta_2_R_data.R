setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(seqinr)
library(pegas)
library(ape)
library(parallel)
source('R_helper_functions.R')

as.alignment.seqinr <- seqinr::as.alignment

files.intersection <- grep(list.files('../fasta_SIs/intersection', full.names = T),
                           pattern = 'X',
                           invert = T,
                           value = T)

alignments.intersection <- lapply(files.intersection,
                                       read.alignment,
                                       format = 'fasta')


# downsample the alignments where required (69 non-admixed ZI lines, 21 of
# these subsampled, 21 least admixed RG lines)

ZI69 <- c('ZI114N', 'ZI117',  'ZI161',  'ZI181',  'ZI184',  'ZI194',  'ZI207',
          'ZI210',  'ZI211',  'ZI213',  'ZI214',  'ZI219',  'ZI232',  'ZI233',
          'ZI235',  'ZI239',  'ZI250',  'ZI252',  'ZI253',  'ZI254N', 'ZI255',
          'ZI264',  'ZI265',  'ZI267', 'ZI268',  'ZI27',   'ZI271',  'ZI292',
          'ZI296',  'ZI303',  'ZI311N', 'ZI320',  'ZI321',  'ZI324',  'ZI332',
          'ZI333',  'ZI339',  'ZI341',  'ZI344',  'ZI348',  'ZI358',  'ZI364',
          'ZI365',  'ZI368',  'ZI378',  'ZI379',  'ZI384',  'ZI386',  'ZI388',
          'ZI398',  'ZI400',  'ZI402',  'ZI418N', 'ZI420',  'ZI437',  'ZI443',
          'ZI447',  'ZI455N', 'ZI456',  'ZI457',  'ZI460',  'ZI476',  'ZI477',
          'ZI486',  'ZI517',  'ZI523',  'ZI527',  'ZI85',  'ZI90')

# I ran "sample(ZI69, 21)" to get 21 random samples from the above vector:
ZI21 <- c("ZI232", "ZI437", "ZI268", "ZI364", "ZI476", "ZI233", "ZI339",
          "ZI400", "ZI210", "ZI398", "ZI456", "ZI321", "ZI368", "ZI265",
          "ZI420", "ZI447", "ZI303", "ZI320", "ZI477", "ZI253", "ZI344")

RG21 <- c('RG34', 'RG19', 'RG7', 'RG24', 'RG36', 'RG37N', 'RG39', 'RG4N', 'RG2',
          'RG32N', 'RG9', 'RG13N', 'RG18N', 'RG22', 'RG25', 'RG28', 'RG3', 'RG33',
          'RG38N', 'RG5', 'RG6N')

alignments.intersection.ZI69.SI.A <- lapply(alignments.intersection,
                                         editAlignment,
                                         paste0('sim|mel|yak|^', paste0(ZI69, collapse = '$|^'), '$'))

alignments.intersection.ZI21.SI.A <- lapply(alignments.intersection,
                                         editAlignment,
                                         paste0('sim|mel|yak|^', paste0(ZI21, collapse = '$|^'), '$'))

alignments.intersection.RG21.SI.A <- lapply(alignments.intersection,
                                         editAlignment,
                                         paste0('sim|mel|yak|^', paste0(RG21, collapse = '$|^'), '$'))

alignments.intersection.MD21.SI.A <- lapply(alignments.intersection,
                                            editAlignment,
                                            'sim|mel|yak|MD')

# save the alignments:
save(list = c('files.intersection',
              'alignments.intersection',
              'alignments.intersection.MD21.SI.A',
              'alignments.intersection.RG21.SI.A',
              'alignments.intersection.ZI21.SI.A',
              'alignments.intersection.ZI69.SI.A'),
     file = '../data/SI_alignments.RData')

load('../data/SI_alignments.RData')

GCcontent.Dmel.SI.A <- sapply(files.intersection, FUN = function(x){
  as.numeric(strsplit(strsplit(x, '/')[[1]][4], '&')[[1]][5])
})
names(GCcontent.Dmel.SI.A) <- NULL

GCcontent.Dsim.SI.A <- sapply(files.intersection, FUN = function(x){
  as.numeric(strsplit(strsplit(strsplit(x, '/')[[1]][4], '&')[[1]][6], '.fa')[[1]][1])
})
names(GCcontent.Dsim.SI.A) <- NULL

GCcontent.mean.SI.A <- apply(cbind(GCcontent.Dmel.SI.A, GCcontent.Dsim.SI.A), 1, FUN = function(rw){
  (rw[1] + rw[2]) / 2
})

# make a function to bin by GC content, with a programmable number of bins:
# This returns the same number of genes per bin
get_GC_bin_index <- function(GC_vector, num_bins){
  
  gc <- GC_vector
  
  gc.rank <- rank(gc, ties.method = 'first')
  
  bin_size <- floor(length(gc) / num_bins) # size of the first n - 1 bins
  final_bin_size <- length(gc) - ((num_bins - 1) * bin_size) # size of the last bin
  bin_sizes <- c(rep(bin_size, num_bins - 1), final_bin_size)
  
  # a loop to population a new vector with the correct bin allocation?
  gc.bin <- vector('integer', length = length(gc)) # vector entirely populated with 0s
  for(i in seq_len(num_bins)){
    # can use the fact that things which haven't been allocated a bin yet are 0s
    # to population the next bin:
    indx <- (gc.rank <= sum(bin_sizes[1:i]) & gc.bin == 0L)
    gc.bin[indx] <- i
  }
  
  return(gc.bin)
}

GCbinindex.5bins.Dsim.SI.A <- get_GC_bin_index(GCcontent.Dsim.SI.A, 5L)
GCbinindex.5bins.Dmel.SI.A <- get_GC_bin_index(GCcontent.Dmel.SI.A, 5L)
GCbinindex.5bins.mean.SI.A <- get_GC_bin_index(GCcontent.mean.SI.A, 5L)


temp.refs.alignments.intersection <- lapply(alignments.intersection, editAlignment, 'sim|mel|yak')

temp.seqs.concat.refs.alignments.intersection.Dsim.5bins <- get_concat_seqs(temp.refs.alignments.intersection,
                                                                      GCbinindex.5bins.Dsim.SI.A)

temp.seqs.concat.refs.alignments.intersection.Dmel.5bins <- get_concat_seqs(temp.refs.alignments.intersection,
                                                                            GCbinindex.5bins.Dmel.SI.A)

temp.seqs.concat.refs.alignments.intersection.mean.5bins <- get_concat_seqs(temp.refs.alignments.intersection,
                                                                            GCbinindex.5bins.mean.SI.A)

temp.seqs.concat.refs.alignments.intersection.diff.5bins <- get_concat_seqs(temp.refs.alignments.intersection,
                                                                            GCbinindex.5bins.diff.SI.A)

GCcontent.concat.5bins.sim.simbins.SI.A <- sapply(temp.seqs.concat.refs.alignments.intersection.Dsim.5bins, FUN = function(x){
  GC.content(as.DNAbin(editAlignment(x, 'sim')))
})

GCcontent.concat.5bins.sim.meanbins.SI.A <- sapply(temp.seqs.concat.refs.alignments.intersection.mean.5bins, FUN = function(x){
  GC.content(as.DNAbin(editAlignment(x, 'sim')))
})

GCcontent.concat.5bins.mel.meanbins.SI.A <- sapply(temp.seqs.concat.refs.alignments.intersection.mean.5bins, FUN = function(x){
  GC.content(as.DNAbin(editAlignment(x, 'mel')))
})

GCcontent.concat.5bins.mel.melbins.SI.A <- sapply(temp.seqs.concat.refs.alignments.intersection.Dmel.5bins, FUN = function(x){
  GC.content(as.DNAbin(editAlignment(x, 'mel')))
})


# save the gc content and bin index stuff:
save(list = c('GCbinindex.5bins.Dsim.SI.A',
              'GCbinindex.5bins.Dmel.SI.A',
              'GCbinindex.5bins.mean.SI.A',
              'GCcontent.Dmel.SI.A',
              'GCcontent.Dsim.SI.A',
              'GCcontent.mean.SI.A',
              'GCcontent.concat.5bins.sim.simbins.SI.A',
              'GCcontent.concat.5bins.sim.meanbins.SI.A',
              'GCcontent.concat.5bins.mel.meanbins.SI.A',
              'GCcontent.concat.5bins.mel.melbins.SI.A'),
     file = '../data/SI_binindices_GCcontent.RData')


# and make a difference GC bin index
GCcontent.diff <- GCcontent.Dsim.SI.A - GCcontent.Dmel.SI.A

GCbinindex.5bins.diff.SI.A <- get_GC_bin_index(GCcontent.diff, 5L)

save(GCbinindex.5bins.diff.SI.A, GCcontent.diff,
     file = '../data/SI_binindices_GCcontent_DIFFERENCES.RData')

