setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# write csv tables of sumstats that can be exported into Word for the supplementary
# information

# load the bin indices and GC contents
load('../data/SI_binindices_GCcontent.RData')
source('../scripts/PAML_parse_output_file.R')

# load the bootstrap data:
load('../data/PAML_bootstraps.RData')

table_write_dir <- '../manuscript/tables/'
# make the folders if they don't already exist
if(!file.exists(table_write_dir)){
  dir.create(table_write_dir,
             recursive = T)
}

myOutputFiles.mean <- list.files(path = '../PAML/mean_SI_A_GCbins/output_files_GC_bins',
                                 pattern = 'output_GC_bin', full.names = T)
myOutputFiles.diff <- list.files(path = '../PAML/diff_SI_A_GCbins/output_files_GC_bins',
                                 pattern = 'output_GC_bin', full.names = T)
myOutputFiles.mel <- list.files(path = '../PAML/mel_SI_A_GCbins/output_files_GC_bins',
                                pattern = 'output_GC_bin', full.names = T)
myOutputFiles.sim <- list.files(path = '../PAML/sim_SI_A_GCbins/output_files_GC_bins',
                                pattern = 'output_GC_bin', full.names = T)

results.mean <- lapply(myOutputFiles.mean, readOutput)
results.diff <- lapply(myOutputFiles.diff, readOutput)
results.mel <- lapply(myOutputFiles.mel, readOutput)
results.sim <- lapply(myOutputFiles.sim, readOutput)

# now write the tables

# mean

mean.df <- data.frame("bin" = rep(1:5, each = 2),
                      
                      "species" = rep(c("sim", "mel"), 5),
                      
                      "num_sites" = as.character(rep(sapply(results.mean, FUN = function(x){x$length}), each = 2)),
                      
                      "N[S>W]" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.mean, FUN = function(x){x$sim_GC_2_AT}),
                                                                        sapply(results.mean, FUN = function(x){x$mel_GC_2_AT})))),
                                        " (",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                        sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                        ", ",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                        sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                        ")"),
                      
                      "N[W>S]" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.mean, FUN = function(x){x$sim_AT_2_CG}),
                                                                        sapply(results.mean, FUN = function(x){x$mel_AT_2_CG})))),
                                        " (",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                        sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                        ", ",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                        sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                        ")"),
                      
                      "N[neu]" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.mean, FUN = function(x){x$sim_N_2_N}),
                                                                        sapply(results.mean, FUN = function(x){x$mel_N_2_N})))),
                                        " (",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                        sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                        ", ",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                        sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                        ")"),
                      
                      "Total changes" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.mean, FUN = function(x){x$sim_GC_2_AT + x$sim_AT_2_CG + x$sim_N_2_N}),
                                                                               sapply(results.mean, FUN = function(x){x$mel_GC_2_AT + x$mel_AT_2_CG + x$mel_N_2_N})))),
                                               " (",
                                               sprintf("%.1f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT + y$sim_AT_2_CG + y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                               sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT + y$mel_AT_2_CG + y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                               ", ",
                                               sprintf("%.1f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT + y$sim_AT_2_CG + y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                               sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT + y$mel_AT_2_CG + y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                               ")"),
                      
                      "GC anc" = paste0(sprintf("%.3f", rep(sapply(results.mean, FUN = function(x){x$GC_ms}), each = 2)),
                                        " (",
                                        sprintf("%.3f", rep(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_ms}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,], each = 2)),
                                        ", ",
                                        sprintf("%.3f", rep(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_ms}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,], each = 2)),
                                        ")"),
                      
                      "GC spp" = paste0(sprintf("%.3f", as.vector(rbind(sapply(results.mean, FUN = function(x){x$GC_sim}),
                                                                        sapply(results.mean, FUN = function(x){x$GC_mel})))),
                                        " (",
                                        sprintf("%.3f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                        sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                        ", ",
                                        sprintf("%.3f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                        sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                        ")"),
                      
                      "GC alt spp" = paste0(sprintf("%.3f", as.vector(rbind(sapply(results.mean, FUN = function(x){x$GC_mel}),
                                                                            sapply(results.mean, FUN = function(x){x$GC_sim})))),
                                            " (",
                                            sprintf("%.3f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                            sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                            ", ",
                                            sprintf("%.3f", as.vector(rbind(sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                            sapply(results.mean.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                            ")"))

write.table(t(mean.df), file = paste0(table_write_dir, "tableS1.mean.txt"), sep = "\t", quote = FALSE, row.names = T, col.names = F)


# species

spp.df <- data.frame("bin" = rep(1:5, each = 2),
                    
                    "species" = rep(c("sim", "mel"), 5),
                    
                    "num_sites" = as.character(as.vector(rbind(sapply(results.sim, FUN = function(x){x$length}),
                                                               sapply(results.mel, FUN = function(x){x$length})))),
                    
                    "N[S>W]" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.sim, FUN = function(x){x$sim_GC_2_AT}),
                                                                      sapply(results.mel, FUN = function(x){x$mel_GC_2_AT})))),
                                      " (",
                                      sprintf("%.1f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                      sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                      ", ",
                                      sprintf("%.1f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                      sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                      ")"),
                    
                    "N[W>S]" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.sim, FUN = function(x){x$sim_AT_2_CG}),
                                                                      sapply(results.mel, FUN = function(x){x$mel_AT_2_CG})))),
                                      " (",
                                      sprintf("%.1f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                      sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                      ", ",
                                      sprintf("%.1f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                      sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                      ")"),
                    
                    "N[neu]" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.sim, FUN = function(x){x$sim_N_2_N}),
                                                                      sapply(results.mel, FUN = function(x){x$mel_N_2_N})))),
                                      " (",
                                      sprintf("%.1f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                      sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                      ", ",
                                      sprintf("%.1f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                      sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                      ")"),
                    
                    "Total changes" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.sim, FUN = function(x){x$sim_GC_2_AT + x$sim_AT_2_CG + x$sim_N_2_N}),
                                                                             sapply(results.mel, FUN = function(x){x$mel_GC_2_AT + x$mel_AT_2_CG + x$mel_N_2_N})))),
                                             " (",
                                             sprintf("%.1f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT + y$sim_AT_2_CG + y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                             sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT + y$mel_AT_2_CG + y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                             ", ",
                                             sprintf("%.1f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT + y$sim_AT_2_CG + y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                             sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT + y$mel_AT_2_CG + y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                             ")"),
                    
                    "GC anc" = paste0(sprintf("%.3f", as.vector(rbind(sapply(results.sim, FUN = function(x){x$GC_ms}),
                                                                      sapply(results.mel, FUN = function(x){x$GC_ms})))),
                                      " (",
                                      sprintf("%.3f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_ms}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                      sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_ms}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                      ", ",
                                      sprintf("%.3f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_ms}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                      sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_ms}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                      ")"),
                    
                    "GC spp" = paste0(sprintf("%.3f", as.vector(rbind(sapply(results.sim, FUN = function(x){x$GC_sim}),
                                                                      sapply(results.mel, FUN = function(x){x$GC_mel})))),
                                      " (",
                                      sprintf("%.3f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                      sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                      ", ",
                                      sprintf("%.3f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                      sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                      ")"),
                    
                    "GC alt spp" = paste0(sprintf("%.3f", as.vector(rbind(sapply(results.sim, FUN = function(x){x$GC_mel}),
                                                                          sapply(results.mel, FUN = function(x){x$GC_sim})))),
                                          " (",
                                          sprintf("%.3f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                          sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                          ", ",
                                          sprintf("%.3f", as.vector(rbind(sapply(results.sim.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                          sapply(results.mel.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                          ")"))

write.table(t(spp.df), file = paste0(table_write_dir, "tableS1.spp.txt"), sep = "\t", quote = FALSE, row.names = T, col.names = F)


# diff

diff.df <- data.frame("bin" = rep(1:5, each = 2),
                      
                      "species" = rep(c("sim", "mel"), 5),
                      
                      "num_sites" = as.character(rep(sapply(results.diff, FUN = function(x){x$length}), each = 2)),
                      
                      "N[S>W]" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.diff, FUN = function(x){x$sim_GC_2_AT}),
                                                                        sapply(results.diff, FUN = function(x){x$mel_GC_2_AT})))),
                                        " (",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                        sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                        ", ",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                        sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                        ")"),
                      
                      "N[W>S]" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.diff, FUN = function(x){x$sim_AT_2_CG}),
                                                                        sapply(results.diff, FUN = function(x){x$mel_AT_2_CG})))),
                                        " (",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                        sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                        ", ",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                        sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_AT_2_CG}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                        ")"),
                      
                      "N[neu]" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.diff, FUN = function(x){x$sim_N_2_N}),
                                                                        sapply(results.diff, FUN = function(x){x$mel_N_2_N})))),
                                        " (",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                        sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                        ", ",
                                        sprintf("%.1f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                        sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                        ")"),
                      
                      "Total changes" = paste0(sprintf("%.1f", as.vector(rbind(sapply(results.diff, FUN = function(x){x$sim_GC_2_AT + x$sim_AT_2_CG + x$sim_N_2_N}),
                                                                               sapply(results.diff, FUN = function(x){x$mel_GC_2_AT + x$mel_AT_2_CG + x$mel_N_2_N})))),
                                               " (",
                                               sprintf("%.1f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT + y$sim_AT_2_CG + y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                               sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT + y$mel_AT_2_CG + y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                               ", ",
                                               sprintf("%.1f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$sim_GC_2_AT + y$sim_AT_2_CG + y$sim_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                               sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$mel_GC_2_AT + y$mel_AT_2_CG + y$mel_N_2_N}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                               ")"),
                      
                      "GC anc" = paste0(sprintf("%.3f", rep(sapply(results.diff, FUN = function(x){x$GC_ms}), each = 2)),
                                        " (",
                                        sprintf("%.3f", rep(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_ms}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,], each = 2)),
                                        ", ",
                                        sprintf("%.3f", rep(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_ms}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,], each = 2)),
                                        ")"),
                      
                      "GC spp" = paste0(sprintf("%.3f", as.vector(rbind(sapply(results.diff, FUN = function(x){x$GC_sim}),
                                                                        sapply(results.diff, FUN = function(x){x$GC_mel})))),
                                        " (",
                                        sprintf("%.3f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                        sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                        ", ",
                                        sprintf("%.3f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                        sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                        ")"),
                      
                      "GC alt spp" = paste0(sprintf("%.3f", as.vector(rbind(sapply(results.diff, FUN = function(x){x$GC_mel}),
                                                                            sapply(results.diff, FUN = function(x){x$GC_sim})))),
                                            " (",
                                            sprintf("%.3f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,],
                                                                            sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[1,]))),
                                            ", ",
                                            sprintf("%.3f", as.vector(rbind(sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_mel}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,],
                                                                            sapply(results.diff.bs, FUN = function(x){vec <- sapply(x, FUN = function(y){y$GC_sim}); CIs <- quantile(vec, probs = c(0.025, 0.975))})[2,]))),
                                            ")"))

write.table(t(diff.df), file = paste0(table_write_dir, "tableS1.diff.txt"), sep = "\t", quote = FALSE, row.names = T, col.names = F)


