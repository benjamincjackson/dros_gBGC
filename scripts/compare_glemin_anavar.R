load('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/data/anavar_results.Rdata')
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_binindices_GCcontent_DIFFERENCES.RData')
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_binindices_GCcontent.RData')
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/glemin_results_DIFFS.RData')
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/glemin_results.RData')


# pop/bin - M0, M1, M0star, M1star. 
GClabelsSimMean <- sprintf("%.2f", signif(GCcontent.concat.5bins.sim.meanbins.SI.A, 2))
GClabelsSimSim <- sprintf("%.2f", signif(GCcontent.concat.5bins.sim.simbins.SI.A, 2))
GClabelsDiff <- 1:5
GClabelsMelMean <- sprintf("%.2f", signif(GCcontent.concat.5bins.mel.meanbins.SI.A, 2))
GClabelsMelMel <- sprintf("%.2f", signif(GCcontent.concat.5bins.mel.melbins.SI.A, 2))




popn <- 'MD'
bin_glemin <- 'mean'
bin_anavar <- 'mean'
pdf(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/plots/glemin_anavar_comparison_', popn, '_', bin_anavar, '.pdf'),
    height = 5, width = 10)
par(lwd = 2,
    cex = 1.6,
    mfrow = c(1,2))
anavar_data <- lapply(c("M1", "none"), FUN = function(x){
  subset(anavar_results, 
         subset = anavar_results$bintype == bin_anavar &
           anavar_results$pop == popn &
           anavar_results$algorithm == 'LBFGS' &
           anavar_results$model == x)
})

glemin_data <- vector('list', 2)
glemin_data[[1]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_no_error')), 
                           FUN = function(x){
                             x$without_error$model1$B
                           })
glemin_data[[2]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_error')), 
                           FUN = function(x){
                             x$with_error$model1$B
                           })
for(i in 1:2){
  if(i == 1){
    model <- 'm1'
  }
  if(i == 2){
    model <- 'm1*'
  }
  plot(x = 1:5,
       y = anavar_data[[i]]$ws_gamma_1,
       main = paste0(popn, ', ', bin_anavar, ' bins, ', model),
       xlab = 'GC content',
       ylab = expression(gamma),
       axes = F,
       pch = 1,
       ylim = c(-3, 6))
  
  points(x = 1:5,
         y = anavar_data[[i]]$sw_gamma_1,
         pch = 2)
  
  points(x = 1:5,
         y = glemin_data[[i]],
         pch = 3)
  
  axis(1, at = 1:5,
       labels = GClabelsSimMean,
       las = 2)
  
  legend('topleft',
         bty = 'n',
         legend = c('anavar_gamma_WS', 'anavar_gamma_SW', 'glemin_gamma'),
         pch = c(1, 2, 3),
         cex = 0.6)
  axis(2)
  box()
}
dev.off()


popn <- 'MD'
bin_glemin <- 'sim'
bin_anavar <- 'species'
pdf(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/plots/glemin_anavar_comparison_', popn, '_', bin_anavar, '.pdf'),
    height = 5, width = 10)
par(lwd = 2,
    cex = 1.6,
    mfrow = c(1,2))
anavar_data <- lapply(c("M1", "none"), FUN = function(x){
  subset(anavar_results, 
         subset = anavar_results$bintype == bin_anavar &
           anavar_results$pop == popn &
           anavar_results$algorithm == 'LBFGS' &
           anavar_results$model == x)
})

glemin_data <- vector('list', 2)
glemin_data[[1]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_no_error')), 
                           FUN = function(x){
                             x$without_error$model1$B
                           })
glemin_data[[2]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_error')), 
                           FUN = function(x){
                             x$with_error$model1$B
                           })
for(i in 1:2){
  if(i == 1){
    model <- 'm1'
  }
  if(i == 2){
    model <- 'm1*'
  }
  plot(x = 1:5,
       y = anavar_data[[i]]$ws_gamma_1,
       main = paste0(popn, ', ', bin_anavar, ' bins, ', model),
       xlab = 'GC content',
       ylab = expression(gamma),
       axes = F,
       pch = 1,
       ylim = c(-3, 6))
  
  points(x = 1:5,
         y = anavar_data[[i]]$sw_gamma_1,
         pch = 2)
  
  points(x = 1:5,
         y = glemin_data[[i]],
         pch = 3)
  
  axis(1, at = 1:5,
       labels = GClabelsSimSim,
       las = 2)
  
  legend('topleft',
         bty = 'n',
         legend = c('anavar_gamma_WS', 'anavar_gamma_SW', 'glemin_gamma'),
         pch = c(1, 2, 3),
         cex = 0.6)
  axis(2)
  box()
}
dev.off()

popn <- 'MD'
bin_glemin <- 'diff'
bin_anavar <- 'diff'
pdf(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/plots/glemin_anavar_comparison_', popn, '_', bin_anavar, '.pdf'),
    height = 5, width = 10)
par(lwd = 2,
    cex = 1.6,
    mfrow = c(1,2))
anavar_data <- lapply(c("M1", "none"), FUN = function(x){
  subset(anavar_results, 
         subset = anavar_results$bintype == bin_anavar &
           anavar_results$pop == popn &
           anavar_results$algorithm == 'LBFGS' &
           anavar_results$model == x)
})

glemin_data <- vector('list', 2)
glemin_data[[1]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_no_error')), 
                           FUN = function(x){
                             x$without_error$model1$B
                           })
glemin_data[[2]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_error')), 
                           FUN = function(x){
                             x$with_error$model1$B
                           })
for(i in 1:2){
  if(i == 1){
    model <- 'm1'
  }
  if(i == 2){
    model <- 'm1*'
  }
  plot(x = 1:5,
       y = anavar_data[[i]]$ws_gamma_1,
       main = paste0(popn, ', ', bin_anavar, ' bins, ', model),
       xlab = 'GC content',
       ylab = expression(gamma),
       axes = F,
       pch = 1,
       ylim = c(-3, 6))
  
  points(x = 1:5,
         y = anavar_data[[i]]$sw_gamma_1,
         pch = 2)
  
  points(x = 1:5,
         y = glemin_data[[i]],
         pch = 3)
  
  axis(1, at = 1:5,
       labels = GClabelsDiff,
       las = 2)
  
  legend('topleft',
         bty = 'n',
         legend = c('anavar_gamma_WS', 'anavar_gamma_SW', 'glemin_gamma'),
         pch = c(1, 2, 3),
         cex = 0.6)
  axis(2)
  box()
}
dev.off()

popn <- 'ZI69'
bin_glemin <- 'mean'
bin_anavar <- 'mean'
pdf(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/plots/glemin_anavar_comparison_', popn, '_', bin_anavar, '.pdf'),
    height = 5, width = 10)
par(lwd = 2,
    cex = 1.6,
    mfrow = c(1,2))
anavar_data <- lapply(c("M1", "none"), FUN = function(x){
  subset(anavar_results, 
         subset = anavar_results$bintype == bin_anavar &
           anavar_results$pop == popn &
           anavar_results$algorithm == 'LBFGS' &
           anavar_results$model == x)
})

glemin_data <- vector('list', 2)
glemin_data[[1]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_no_error')), 
                           FUN = function(x){
                             x$without_error$model1$B
                           })
glemin_data[[2]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_error')), 
                           FUN = function(x){
                             x$with_error$model1$B
                           })
for(i in 1:2){
  if(i == 1){
    model <- 'm1'
  }
  if(i == 2){
    model <- 'm1*'
  }
  plot(x = 1:5,
       y = anavar_data[[i]]$ws_gamma_1,
       main = paste0(popn, ', ', bin_anavar, ' bins, ', model),
       xlab = 'GC content',
       ylab = expression(gamma),
       axes = F,
       pch = 1,
       ylim = c(-3, 6))
  
  points(x = 1:5,
         y = anavar_data[[i]]$sw_gamma_1,
         pch = 2)
  
  points(x = 1:5,
         y = glemin_data[[i]],
         pch = 3)
  
  axis(1, at = 1:5,
       labels = GClabelsMelMean,
       las = 2)
  
  legend('topleft',
         bty = 'n',
         legend = c('anavar_gamma_WS', 'anavar_gamma_SW', 'glemin_gamma'),
         pch = c(1, 2, 3),
         cex = 0.6)
  axis(2)
  box()
}
dev.off()


popn <- 'ZI69'
bin_glemin <- 'mel'
bin_anavar <- 'species'
pdf(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/plots/glemin_anavar_comparison_', popn, '_', bin_anavar, '.pdf'),
    height = 5, width = 10)
par(lwd = 2,
    cex = 1.6,
    mfrow = c(1,2))
anavar_data <- lapply(c("M1", "none"), FUN = function(x){
  subset(anavar_results, 
         subset = anavar_results$bintype == bin_anavar &
           anavar_results$pop == popn &
           anavar_results$algorithm == 'LBFGS' &
           anavar_results$model == x)
})

glemin_data <- vector('list', 2)
glemin_data[[1]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_no_error')), 
                           FUN = function(x){
                             x$without_error$model1$B
                           })
glemin_data[[2]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_error')), 
                           FUN = function(x){
                             x$with_error$model1$B
                           })
for(i in 1:2){
  if(i == 1){
    model <- 'm1'
  }
  if(i == 2){
    model <- 'm1*'
  }
  plot(x = 1:5,
       y = anavar_data[[i]]$ws_gamma_1,
       main = paste0(popn, ', ', bin_anavar, ' bins, ', model),
       xlab = 'GC content',
       ylab = expression(gamma),
       axes = F,
       pch = 1,
       ylim = c(-3, 6))
  
  points(x = 1:5,
         y = anavar_data[[i]]$sw_gamma_1,
         pch = 2)
  
  points(x = 1:5,
         y = glemin_data[[i]],
         pch = 3)
  
  axis(1, at = 1:5,
       labels = GClabelsMelMel,
       las = 2)
  
  legend('topleft',
         bty = 'n',
         legend = c('anavar_gamma_WS', 'anavar_gamma_SW', 'glemin_gamma'),
         pch = c(1, 2, 3),
         cex = 0.6)
  axis(2)
  box()
}
dev.off()

popn <- 'ZI69'
bin_glemin <- 'diff'
bin_anavar <- 'diff'
pdf(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/plots/glemin_anavar_comparison_', popn, '_', bin_anavar, '.pdf'),
    height = 5, width = 10)
par(lwd = 2,
    cex = 1.6,
    mfrow = c(1,2))
anavar_data <- lapply(c("M1", "none"), FUN = function(x){
  subset(anavar_results, 
         subset = anavar_results$bintype == bin_anavar &
           anavar_results$pop == popn &
           anavar_results$algorithm == 'LBFGS' &
           anavar_results$model == x)
})

glemin_data <- vector('list', 2)
glemin_data[[1]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_no_error')), 
                           FUN = function(x){
                             x$without_error$model1$B
                           })
glemin_data[[2]] <- sapply(get(paste0('results_', popn, '.SI.A.5GCbins.', bin_glemin, '_error')), 
                           FUN = function(x){
                             x$with_error$model1$B
                           })
for(i in 1:2){
  if(i == 1){
    model <- 'm1'
  }
  if(i == 2){
    model <- 'm1*'
  }
  plot(x = 1:5,
       y = anavar_data[[i]]$ws_gamma_1,
       main = paste0(popn, ', ', bin_anavar, ' bins, ', model),
       xlab = 'GC content',
       ylab = expression(gamma),
       axes = F,
       pch = 1,
       ylim = c(-3, 6))
  
  points(x = 1:5,
         y = anavar_data[[i]]$sw_gamma_1,
         pch = 2)
  
  points(x = 1:5,
         y = glemin_data[[i]],
         pch = 3)
  
  axis(1, at = 1:5,
       labels = GClabelsDiff,
       las = 2)
  
  legend('topright',
         bty = 'n',
         legend = c('anavar_gamma_WS', 'anavar_gamma_SW', 'glemin_gamma'),
         pch = c(1, 2, 3),
         cex = 0.6)
  axis(2)
  box()
}
dev.off()

