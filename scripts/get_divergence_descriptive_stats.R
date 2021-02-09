setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# calculate descriptive population genetics statistics for within-population
# diversity and between-species/population divergence from the SI sites
# that we use for all the analyses

source('../scripts/PAML_parse_output_file.R')
myOutputFile <- list.files(path = '../PAML/SI_A_ALL/run_1/',
                                 pattern = 'output', full.names = T)

point_estimate <- readOutput(myOutputFile)
point_estimate$sim_blength
point_estimate$mel_blength

point_estimate$sim_blength + point_estimate$mel_blength

source('../scripts/PAML_parse_output_file.R')

myOutputFiles.all.bs <- list.files(path = '../PAML/bootstraps/all/outputfiles/',
                                   pattern = 'output', full.names = T)

results.all.bs <- lapply(myOutputFiles.all.bs, readOutput)

ms.all.CIs <- quantile(
  sapply(results.all.bs, FUN = function(x){
  sim <- x$sim_blength
  mel <- x$mel_blength
  vec <- sim + mel
  return(vec)
}), probs = c(0.025, 0.975))

m.all.CIs <- quantile(
  sapply(results.all.bs, FUN = function(x){
    mel <- x$mel_blength
    return(mel)
  }), probs = c(0.025, 0.975))

s.all.CIs <- quantile(
  sapply(results.all.bs, FUN = function(x){
    sim <- x$sim_blength
    return(sim)
  }), probs = c(0.025, 0.975))


source('../scripts/R_helper_functions.R')
load('../data/SI_alignments.RData')
load('../data/SI_binindices_GCcontent.RData')

library(ape)

raw.divergence.simbins <- sapply(lapply(get_concat_seqs(alignments.intersection, GCbinindex.5bins.Dsim.SI.A),
                                 editAlignment, 'sim|mel'), FUN = function(x){
                                   dist.dna(as.DNAbin(x), model = 'raw')
                                 })

raw.divergence.melbins <- sapply(lapply(get_concat_seqs(alignments.intersection, GCbinindex.5bins.Dmel.SI.A),
                                 editAlignment, 'sim|mel'), FUN = function(x){
                                   dist.dna(as.DNAbin(x), model = 'raw')
                                 })

raw.divergence.meanbins <- sapply(lapply(get_concat_seqs(alignments.intersection, GCbinindex.5bins.mean.SI.A),
                                 editAlignment, 'sim|mel'), FUN = function(x){
                                   dist.dna(as.DNAbin(x), model = 'raw')
                                 })

raw.divergence.all.concat <- sapply(lapply(get_concat_seqs(alignments.intersection, rep(1, length(alignments.intersection))),
                                editAlignment, 'sim|mel'), FUN = function(x){
                                  dist.dna(as.DNAbin(x), model = 'raw')
                                })

#