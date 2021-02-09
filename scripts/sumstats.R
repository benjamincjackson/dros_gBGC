pi_from_uSFS <- function(SFS, per_site = TRUE) {
  
  if(sum(SFS) == 0) {
    return(-99)
  }
  
  N = length(SFS) - 1
  binom = (N * (N - 1)) / 2
  
  temp <- vector('numeric', length = length(SFS) - 1)
  for(j in seq_along(SFS)){
    if(j == 1){ next }
    i <- j - 1
    temp[i] <- (i * (N - i) * SFS[j]) / binom
  }
  
  pi = sum(temp)
  
  if(per_site == TRUE) {
    return(pi/sum(SFS))
  }
  else {
    return(pi)
  }
}

theta_W_from_uSFS <- function(SFS, per_site = TRUE) {
  # if(identical(length(SFS), 0)) {
  #   return(-99)
  # }
  N <- length(SFS) - 1
  S <- sum(SFS[2:N]) ## takes the interior of the SFS, gets S
  
  harmonic <- sum( 1 / seq_len(N - 1))
  
  if(per_site == TRUE) {
    return(S / (harmonic * sum(SFS)))
  }
  else {
    return(S / (harmonic))
  }
}

tajimas_D_from_uSFS <- function(SFS) {
  th_pi <- pi_from_uSFS(SFS, per_site=FALSE)
  n <- length(SFS) - 1
  S <- sum(SFS[2:n]) ## takes the interior of the SFS, gets S
  
  tmp <- 1:(n - 1)
  a1 <- sum(1/tmp)
  a2 <- sum(1/tmp^2)
  b1 <- (n + 1)/(3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
  e1 <- c1/a1
  e2 <- c2/(a1^2 + a2)
  D <- (th_pi - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
  return(D)
}

delta_pi_from_uSFS <- function(SFS){
  th_pi = pi_from_uSFS(SFS, per_site=FALSE)
  N = length(SFS)-1
  S = sum(SFS[2:(length(SFS) - 1)])
  
  tmp <- 1:(N - 1)
  a1 <- sum(1/tmp)
  
  delta_pi = (th_pi / S) - (1 / a1)
  return(delta_pi)
}

GCcontent_from_freqs <-function(input_file){
  # input file should look like:
  # 16:4:0:0
  # 0:20:0:0
  # 0:0:0:20
  # ...
  # A:T:G:C
  
  con <- file(input_file, open="r")
  lin <-readLines(con)
  close(con)
  
  lin_split <- lapply(lin, FUN = function(x){
    as.integer(strsplit(x, split = ':')[[1]])
  })
  
  myMat <- do.call(rbind, lin_split)
  
  GC_count <- sum(myMat[,3:4])
  total_count <- sum(myMat)
  
  GC_content <- GC_count / total_count
  
  return(GC_content)
}

prop_singletons_from_uSFS <- function(SFS){
  prop <- SFS[2] / sum(SFS[2:(length(SFS) - 1)])
  return(prop)
}


bootstrap_sumstats <- function(list_of_SFSs, n = 1000){
  max_SFS_length <- max(sapply(list_of_SFSs, length))
  SFSs_samelength <- list_of_SFSs[sapply(list_of_SFSs, FUN = function(x){ length(x) == max_SFS_length })]
  summed_sfs <- apply(do.call(rbind, SFSs_samelength), 2, sum)
  
  num_seg_sites <- sum(summed_sfs[2:length(summed_sfs)] - 1)
  
  index <- seq_along(SFSs_samelength)
  
  pi <- vector('numeric', length = n)
  theta <- vector('numeric', length = n)
  tajD <- vector('numeric', length = n)
  delta_pi <- vector('numeric', length = n)
  prop_singletons  <- vector('numeric', length = n)
  
  for(j in seq_len(n)){
    i <- 0
    sfs <- rep(0, max_SFS_length)
    while(i < num_seg_sites){
      random_index <- sample(index, 1)
      random_sfs <- SFSs_samelength[[random_index]]
      n_seg_sites_random_sfs <- sum(random_sfs[2:(length(random_sfs) - 1)])
      sfs <- sfs + random_sfs
      i <- i + n_seg_sites_random_sfs
    }
    pi[j] <- pi_from_uSFS(sfs)
    theta[j] <- theta_W_from_uSFS(sfs)
    tajD[j] <- tajimas_D_from_uSFS(sfs)
    delta_pi[j] <- delta_pi_from_uSFS(sfs)
    prop_singletons[j] <- prop_singletons_from_uSFS(sfs)
  }
  return(data.frame('pi' = pi, 'theta' = theta, 'D' = tajD, 'delta_pi' = delta_pi, 'prop_single' = prop_singletons))
}

get_SFS_from_file <- function(fileName, splt = ','){
  con <- file(fileName, open="r")
  lin <-readLines(con)
  sfs <- lin[1]
  sfs <- as.numeric(strsplit(sfs, split = splt)[[1]])
  close(con)
  return(sfs)
}

sum_SFS <- function(listofSFSs){
  lengths <- sapply(listofSFSs, length)
  maxLength <- max(lengths)
  newList <- listofSFSs[lengths == maxLength]
  temp <- do.call(rbind, newList)
  summedSFS <- apply(temp, 2, sum)
  return(summedSFS)
}

get_SFS_from_alignment <- function(al, popstring){
  newAl <- editAlignment(al, popstring)
  mymat <- makeMatrix(newAl)
  n <- newAl$nb
  sfs <- vector('numeric', length = n + 1)

  indx <- apply(mymat, 2, FUN = function(x){
    all(all(levels(as.factor(x)) != 'n'),
        nlevels(as.factor(x)) <= 2)
  })
  
  filtmat <- mymat[,indx, drop = F]
  
  for(i in seq_len(ncol(filtmat))){
    cl <- filtmat[,i]
    if(nlevels(as.factor(cl)) == 1){
      sfs[1] <- sfs[1] + 1
    }
    else {
      ac1 <- summary(as.factor(cl))[[1]]
      ac2 <- summary(as.factor(cl))[[2]]
      MAC <- min(ac1, ac2)
      sfs[(MAC + 1)] <- sfs[(MAC + 1)] + 1
    }
  }
  
  return(sfs)
}









#