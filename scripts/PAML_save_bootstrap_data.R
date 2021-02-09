setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# read the PAML bootstrap data and write it to an RData object, for use later
# in plotting

source('../scripts/PAML_parse_output_file.R')

myOutputFiles.mel.bs <- lapply(list.dirs(path = '../PAML/bootstraps/mel/outputfiles')[2:6],
                               FUN = function(mydir){
                                 list.files(path = mydir,
                                            pattern = 'output', full.names = T)
                               })

myOutputFiles.sim.bs <- lapply(list.dirs(path = '../PAML/bootstraps/sim/outputfiles')[2:6],
                               FUN = function(mydir){
                                 list.files(path = mydir,
                                            pattern = 'output', full.names = T)
                               })

myOutputFiles.mean.bs <- lapply(list.dirs(path = '../PAML/bootstraps/mean/outputfiles')[2:6],
                                FUN = function(mydir){
                                  list.files(path = mydir,
                                             pattern = 'output', full.names = T)
                                })

myOutputFiles.diff.bs <- lapply(list.dirs(path = '../PAML/bootstraps/diff/outputfiles')[2:6],
                                FUN = function(mydir){
                                  list.files(path = mydir,
                                             pattern = 'output', full.names = T)
                                })


results.mel.bs <- lapply(myOutputFiles.mel.bs,
                         FUN = function(x){
                           stats <- lapply(x, readOutput)
                           return(stats)
                         })

results.sim.bs <- lapply(myOutputFiles.sim.bs,
                         FUN = function(x){
                           stats <- lapply(x, readOutput)
                           return(stats)
                         })

results.mean.bs <- lapply(myOutputFiles.mean.bs,
                          FUN = function(x){
                            stats <- lapply(x, readOutput)
                            return(stats)
                          })

results.diff.bs <- lapply(myOutputFiles.diff.bs,
                          FUN = function(x){
                            stats <- lapply(x, readOutput)
                            return(stats)
                          })

save(list = ls(pattern = 'results'),
     file = '../data/PAML_bootstraps.RData')
