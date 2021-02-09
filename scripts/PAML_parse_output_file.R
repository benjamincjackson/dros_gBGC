# write a function to parse the 'output' file of PAML

readOutput <- function(filepath){
  myCon <- file(filepath)
  myFile <- readLines(con = myCon)
  close(myCon)
  
  L <- as.numeric(strsplit(myFile[[1]], split = '\\s+')[[1]][3])
  
  GCstart <- grep('(frequency parameters for branches)', myFile)
  fGC <- myFile[GCstart:(GCstart + 7)]
  msGCline <- fGC[grep('Node # 5 ', fGC)]
  msC <- as.numeric(strsplit(msGCline, split = '\\s+')[[1]][12])
  msG <- as.numeric(strsplit(msGCline, split = '\\s+')[[1]][14])

  mGCline <- fGC[grep('Node # 2 ', fGC)]
  mC <- as.numeric(strsplit(mGCline, split = '\\s+')[[1]][12])
  mG <- as.numeric(strsplit(mGCline, split = '\\s+')[[1]][14])
  
  sGCline <- fGC[grep('Node # 1 ', fGC)]
  sC <- as.numeric(strsplit(sGCline, split = '\\s+')[[1]][12])
  sG <- as.numeric(strsplit(sGCline, split = '\\s+')[[1]][14])
  
  mainStart <- grep('Expected numbers of nucleotide changes on branches', myFile)
  f <- myFile[mainStart:length(myFile)]
  simStart <- grep('sim', f)
  melStart <- grep('mel', f)
  
  simBranchLench <- as.numeric(strsplit(strsplit(f[simStart], split = '\\s+')[[1]][7], ',')[[1]][1])
  melBranchLench <- as.numeric(strsplit(strsplit(f[melStart], split = '\\s+')[[1]][7], ',')[[1]][1])
  
  simLines <- f[(simStart + 4):(simStart + 7)]
  simMat <- matrix(
    sapply(simLines, FUN = function(x){
      as.numeric(strsplit(x, split = '\\s+')[[1]][6:9])
    }), nrow = 4, ncol = 4, byrow = T,
    dimnames = list(c('T', 'C', 'A', 'G'), c('T', 'C', 'A', 'G')))
  
  sim_N_2_N <- simMat['A', 'T'] + simMat['T', 'A'] + simMat['C', 'G'] + simMat['G', 'C']
  sim_AT_2_CG <- simMat['A', 'C'] + simMat['A', 'G'] + simMat['T', 'C'] + simMat['T', 'G']
  sim_GC_2_AT <- simMat['C', 'A'] + simMat['C', 'T'] + simMat['G', 'A'] + simMat['G', 'T']
  
  melLines <- f[(melStart + 4):(melStart + 7)]
  melMat <- matrix(
    sapply(melLines, FUN = function(x){
      as.numeric(strsplit(x, split = '\\s+')[[1]][6:9])
    }), nrow = 4, ncol = 4, byrow = T,
    dimnames = list(c('T', 'C', 'A', 'G'), c('T', 'C', 'A', 'G')))
  
  mel_N_2_N <- melMat['A', 'T'] + melMat['T', 'A'] + melMat['C', 'G'] + melMat['G', 'C']
  mel_AT_2_CG <- melMat['A', 'C'] + melMat['A', 'G'] + melMat['T', 'C'] + melMat['T', 'G']
  mel_GC_2_AT <- melMat['C', 'A'] + melMat['C', 'T'] + melMat['G', 'A'] + melMat['G', 'T']
 
  return(list(
    'length' = L,
    'sim_blength' = simBranchLench,
    'mel_blength' = melBranchLench,
    'GC_ms' = msG + msC,
    'GC_mel' = mG + mC,
    'GC_sim' = sG + sC,
    'sim_N_2_N' = sim_N_2_N,
    'sim_AT_2_CG' = sim_AT_2_CG,
    'sim_GC_2_AT' = sim_GC_2_AT,
    'mel_N_2_N' = mel_N_2_N,
    'mel_AT_2_CG' = mel_AT_2_CG,
    'mel_GC_2_AT' = mel_GC_2_AT,
    'sim_AT_2_CG_rate' = sim_AT_2_CG / (L - (L * (msG + msC))),
    'sim_GC_2_AT_rate' = sim_GC_2_AT / (L * (msG + msC)),
    'mel_AT_2_CG_rate' = mel_AT_2_CG / (L - (L * (msG + msC))),
    'mel_GC_2_AT_rate' = mel_GC_2_AT / (L * (msG + msC)),
    'sim_kappa' = (sim_GC_2_AT / (L * (msG + msC))) / (sim_AT_2_CG / (L - (L * (msG + msC)))),
    'mel_kappa' = (mel_GC_2_AT / (L * (msG + msC))) / (mel_AT_2_CG / (L - (L * (msG + msC)))),
    'sim_ratio_counts' = sim_AT_2_CG / sim_GC_2_AT,
    'mel_ratio_counts' = mel_AT_2_CG / mel_GC_2_AT
  ))
}