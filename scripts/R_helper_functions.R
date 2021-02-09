as.alignment.seqinr <- seqinr::as.alignment

compareToRef <- function(targ, ref){
  identical(nchar(ref), nchar(targ))
}

compareLength <- function(alignment){
  reference <- alignment$seq[1]
  targets <- alignment$seq[-1]
  target_equals_ref <- sapply(targets, compareToRef, ref = reference)
  return(all(target_equals_ref))
}

compareLengthVec <- function(vectorOfSequences){
  reference <- vectorOfSequences[1]
  targets <- vectorOfSequences[-1]
  target_equals_ref <- sapply(targets, compareToRef, ref = reference)
  return(all(target_equals_ref))
}


getSequenceList <- function(alignment, popstring){
  index <- grep(popstring, alignment$nam, value = FALSE)
  return(index)
}

editAlignment <- function(oldal, popstr){
  newalignment <- as.alignment.seqinr(nb = length(getSequenceList(alignment = oldal, popstring = popstr)),
                                      nam = oldal$nam[getSequenceList(alignment = oldal, popstring = popstr)],
                                      seq = oldal$seq[getSequenceList(alignment = oldal, popstring = popstr)],
                                      com = oldal$com[getSequenceList(alignment = oldal, popstring = popstr)])
  return(newalignment)
}

concatenateSequences <- function(vectorOfSequences){
  require(seqinr)
  combinedSequences <- c2s(unlist(lapply(vectorOfSequences, s2c)))
  return(combinedSequences)
}

stripSequences <- function(sequences){
  justSequences <- lapply(sequences, FUN = function(al){al$seq})
  return(justSequences)
}

transposeSequences <- function(listOfListOfSequences){
  #the below code transposes a list of list of sequences
  #see http://stackoverflow.com/questions/16179197/transpose-a-list-of-lists
  #assumes there are the same number of elements in [[1]] as the rest of the list
  newList <- lapply(seq_along(listOfListOfSequences[[1]]), FUN = function(i){
    unlist(lapply(listOfListOfSequences, FUN = '[[', i))})
  return(newList)
}

makeGiantAlignmentCDS <- function(sequences){

  listOfTransposedSeqs <- transposeSequences(stripSequences(sequences))
  concatSeq <- lapply(listOfTransposedSeqs, concatenateSequences)

  nb  <- length(concatSeq)
  nam <- sequences[[1]]$nam

  newGiantAlignment <- seqinr::as.alignment(nb = nb,
                                            nam = nam,
                                            seq = concatSeq,
                                            com = NULL)
  return(newGiantAlignment)
}


get_concat_seqs <- function(sequences, binIndex){
  #get the bin index:
  bin_index <- binIndex
  #make an empty vector to fill with the new concatenated sequences
  concat_seqs <- vector('list', length = nlevels(as.factor(bin_index)))
  #fill the vector:
  for(i in seq_len(nlevels(as.factor(bin_index)))){
    sequences_temp <- sequences[bin_index == i]
    sequences_temp <- sequences_temp[!sapply(sequences_temp, is.null)] # because some of the alignments may be empty, which screws everything up
    concat_seqs[[i]] <- makeGiantAlignmentCDS(sequences_temp)
  }
  return(concat_seqs)
}

# this makes a wide (landscape) matrix of the sequences
makeMatrix <- function(alignment){ ##NB THIS ONLY WORKS IF LENGTHS OF ALL SEQS ARE EQUAL
  require(seqinr)
  charstrlist <- lapply(alignment$seq, s2c)
  charstrs <- unlist(charstrlist)
  numofseqs <- length(alignment$seq) # no. of sequences
  m <- matrix(charstrs, nrow = numofseqs, byrow = T)
  return(m)
}

# makes sequences from rows of a matrix
reconstructSequences <- function(mat){
  newSeqs <- apply(mat, 1, FUN = function(rw){ paste0(rw, collapse = '') })
  return(newSeqs)
}

getEverythingPAML <- function(probTable, extantNode){
  # First we want to get the state at the EXTANT (mel or sim) node, which is the first 
  # nucleotide of the triplet for sim, and the second for mel
  extantNuc <- sapply(rownames(probTable), FUN = function(x){
    strsplit(x, '')[[1]][extantNode]
  })
  names(extantNuc) <- NULL
  # and cbind it to the other data to make a dataframe
  probTable <- data.frame(probTable, extantNuc, stringsAsFactors = F)
  # now we want to get the changes out of this table
  
  # subset the whole table to get just the sites where the sim nucleotide is C, and the
  # the mel-sim site was A (with a given probability)
  A_2_C_table <- subset(probTable[,c(1:2, 6)], extantNuc == 'C')
  # then you can multiply the freq and the prob columns together to get the (probabilistic)
  # total number of sites
  A_2_C_total <- sum(A_2_C_table$freq * A_2_C_table$A)
  
  # repeat for A to G changes, T to C changes, and T to G changes
  A_2_G_table <- subset(probTable[,c(1:2, 6)], extantNuc == 'G')
  A_2_G_total <- sum(A_2_G_table$freq * A_2_G_table$A)
  
  T_2_C_table <- subset(probTable[,c(1, 3, 6)], extantNuc == 'C')
  T_2_C_total <- sum(T_2_C_table$freq * T_2_C_table$T)
  
  T_2_G_table <- subset(probTable[,c(1, 3, 6)], extantNuc == 'G')
  T_2_G_total <- sum(T_2_G_table$freq * T_2_G_table$T)
  
  # then the total number of AT_2_GC changes is these four totals added together:
  AT_2_GC_total <- A_2_C_total + A_2_G_total + T_2_C_total + T_2_G_total
  
  # This is the total number of AT sites in the mel-sim ancestor, which we get by
  # multiplying the frequencies of the sites by the probability of having an A or
  # a T in the ancestor at each of those sites
  AT_total_sites <- sum(probTable$A * probTable$freq) + sum(probTable$T * probTable$freq)
  
  # Then the rate is the total number of changes multiplied by the total number of sites
  AT_2_GC_rate <- AT_2_GC_total / AT_total_sites
  
  ### THEN EVERYTHING AGAIN FOR THE GC > AT CHANGES ###
  
  C_2_A_table <- subset(probTable[,c(1, 5, 6)], extantNuc == 'A')
  C_2_A_total <- sum(C_2_A_table$freq * C_2_A_table$C)
  
  C_2_T_table <- subset(probTable[,c(1, 5, 6)], extantNuc == 'T')
  C_2_T_total <- sum(C_2_T_table$freq * C_2_T_table$C)
  
  G_2_A_table <- subset(probTable[,c(1, 4, 6)], extantNuc == 'A')
  G_2_A_total <- sum(G_2_A_table$freq * G_2_A_table$G)
  
  G_2_T_table <- subset(probTable[,c(1, 4, 6)], extantNuc == 'T')
  G_2_T_total <- sum(G_2_T_table$freq * G_2_T_table$G)
  
  # then the total number of AT_2_GC changes is these four totals added together:
  GC_2_AT_total <- C_2_A_total + C_2_T_total + G_2_A_total + G_2_T_total
  
  # This is the total number of GC sites in the mel-sim ancestor, which we get by
  # multiplying the frequencies of the sites by the probability of having an G or
  # a C in the ancestor at each of those sites
  GC_total_sites <- sum(probTable$G * probTable$freq) + sum(probTable$C * probTable$freq)
  
  # Then the rate is the total number of changes multiplied by the total number of sites
  GC_2_AT_rate <- GC_2_AT_total / GC_total_sites
  
  mydf <- data.frame('GC_2_AT_total' = GC_2_AT_total,
                     'AT_2_GC_total' = AT_2_GC_total,
                     'GC_total_sites' = GC_total_sites,
                     'AT_total_sites' = AT_total_sites,
                     'GCcontent' = GC_total_sites / (GC_total_sites + AT_total_sites),
                     'GC_2_AT_rate' = GC_2_AT_rate,
                     'AT_2_GC_rate' = AT_2_GC_rate,
                     'WS_SW_ratio' = AT_2_GC_total / GC_2_AT_total)
  
  return(mydf)
}
