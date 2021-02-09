setwd(dirname(rstudioapi::getSourceEditorContext()$path))

load('../data/SI_binindices_GCcontent.RData')
source('../scripts/PAML_parse_output_file.R')


# myOutputFiles.mean <- list.files(path = '../PAML/mean_SI_A_GCbins/output_files_GC_bins',
#                                  pattern = 'output_GC_bin', full.names = T)

# results.mean <- lapply(myOutputFiles.mean, readOutput)
# 
# mel.a.1 <- data.frame(S2W = sum(sapply(results.mean, FUN = function(x){x$mel_GC_2_AT})),
#                       W2S = sum(sapply(results.mean, FUN = function(x){x$mel_AT_2_CG})),
#                       row.names = c('A'))
# 
# chisq.test(mel.a.1)
# 
# sim.a.1 <- data.frame(S2W = sum(sapply(results.mean, FUN = function(x){x$sim_GC_2_AT})),
#                       W2S = sum(sapply(results.mean, FUN = function(x){x$sim_AT_2_CG})),
#                       row.names = c('A'))
# 
# chisq.test(sim.a.1)


all <- readOutput('../PAML/SI_A_ALL/run_1/output')
all.table.sim <- data.frame(S2W = all$sim_GC_2_AT,
                            W2S = all$sim_AT_2_CG,
                            row.names = c('A'))
chisq.test(all.table.sim)

all.table.mel <- data.frame(S2W = all$mel_GC_2_AT,
                            W2S = all$mel_AT_2_CG,
                            row.names = c('A'))
chisq.test(all.table.mel)
