# read the output files produced by anavar for plotting

# # NB Sometimes there's a bug in Kai's code: there is no newline character 
# # separating the parameter names from results, so the highest lnL result
# # is on the same line as the parameter names
# # IVE JUST MANUALLY EDITED ALL OF THE BROKEN RESULTS FILES
# 
# read_BROKEN_anavar_input <- function(filepath, bintype, binnumber, pop, model, algorithm){
#   
#   myCon <- file(filepath)
#   myLines <- readLines(myCon)
#   
#   brokenHeader <- myLines[7]
#   subHeader <- sub(pattern = 'lnL', replacement = 'lnLBREAK', x = brokenHeader)
#   splitHeader <- strsplit(subHeader, split = 'BREAK')[[1]]
#   
#   fixedHeader <- strsplit(splitHeader[1], split = '\\s+')[[1]]
#   fixedBestResult <- as.numeric(strsplit(splitHeader[2], split = '\\s+')[[1]])
#   
#   info <- c(bintype, binnumber, pop, model, algorithm)
#   infonames <- c('bintype', 'binnumber', 'pop', 'model', 'algorithm')
#   
#   myDF <- data.frame(t(c(info, fixedBestResult)), stringsAsFactors = FALSE)
#   names(myDF) <- c(infonames, fixedHeader)
#   
#   close(myCon)
#   
#   # myDF <- data.frame(bestresult)
#   return(myDF)
# }
# 
# read_anavar_input <- function(filepath, bintype, binnumber, pop, model, algorithm){
# 
#   myCon <- file(filepath)
#   myLines <- readLines(myCon)
# 
#   temp <- strsplit(myLines[7], split = '\\s+')[[1]]
# 
#   myHeader <- strsplit(myLines[7], split = '\\s+')[[1]][c(1:12, length(temp))]
#   myBestResult <- strsplit(myLines[8], split = '\\s+')[[1]][c(1:12, length(temp))]
# 
#   info <- c(bintype, binnumber, pop, model, algorithm)
#   infonames <- c('bintype', 'binnumber', 'pop', 'model', 'algorithm')
# 
#   myDF <- data.frame(t(c(info, myBestResult)), stringsAsFactors = FALSE)
#   names(myDF) <- c(infonames, myHeader)
# 
#   close(myCon)
# 
#   return(myDF)
# }
# 
# # test <- read_BROKEN_anavar_input('~/Desktop/results_LBFGS.txt', bintype = 'species', binnumber = '1', pop = 'MD', model = 'M0', algorithm = 'LGBTS')
# # test$ws_gamma_1
# 
# myList <- vector('list', length = 300)
# i = 1
# for(bintype in c('mean', 'diff', 'species')){
#   for(binnumber in 1:5){
#     for(pop in c('MD', 'ZI69')){
#       for(model in c('none', 'M0',  'M1', 'M0star', 'M1star')){
#         for(algorithm in c('LBFGS', 'SLSQP')){
# 
#           fp <- paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/anavar/',
#                        bintype, '/',
#                        pop, '/',
#                        binnumber, '/',
#                        model, '/',
#                        'results_', algorithm, '.txt', collapse = '')
# 
#           myList[[i]] <- read_anavar_input(filepath = fp,
#                                                   bintype = bintype,
#                                                   binnumber = binnumber,
#                                                   pop = pop,
#                                                   model = model,
#                                                   algorithm = algorithm)
# 
# 
#           i <- i+1
# 
#         }}}}}
# 
# anavar_results <- do.call(rbind, myList)
# 
# # test for broken headers:
# for(i in seq_along(myList)){
#   l <- myList[[i]]
#   if(i==1){
#     n_old = names(l)
#     next
#   }
#   n_new = names(l)
#   if (!(identical(n_new, n_old))){
#     print(i)
#   }
#   n_old = n_new
# }


# # let's save this:
# save(anavar_results, file = '/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/data/anavar_results.Rdata')

load('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/data/anavar_results.Rdata')
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_binindices_GCcontent_DIFFERENCES.RData')
load('~/Dropbox/Biology_projects/dros_X_A_mut/data/SI_binindices_GCcontent.RData')

# a function for plotting the anavar results
plot_anavar_results <- function(pop_in, bintype_in, algorithm_in){
  
  if(pop_in == 'MD'){
    if(bintype_in == 'mean'){
      GClabels <- sprintf("%.2f", signif(GCcontent.concat.5bins.sim.meanbins.SI.A, 2))
    }
    if(bintype_in == 'species'){
      GClabels <- sprintf("%.2f", signif(GCcontent.concat.5bins.sim.simbins.SI.A, 2))
    }
    if(bintype_in == 'diff'){
      GClabels <- 1:5
    }
  }
  if(pop_in == 'ZI69'){
    if(bintype_in == 'mean'){
      GClabels <- sprintf("%.2f", signif(GCcontent.concat.5bins.mel.meanbins.SI.A, 2))
    }
    if(bintype_in == 'species'){
      GClabels <- sprintf("%.2f", signif(GCcontent.concat.5bins.mel.melbins.SI.A, 2))
    }
    if(bintype_in == 'diff'){
      GClabels <- 1:5
    }
  }
  
  data <- subset(anavar_results, 
                 subset = anavar_results$bintype == bintype_in &
                   anavar_results$pop == pop_in &
                   anavar_results$algorithm == algorithm_in)
  
  m1 <- subset(data, subset = data$model == 'M1')
  m1star <- subset(data, subset = data$model == 'M1star')
  none <- subset(data, subset = data$model == 'none')
  
  pdf(paste0('/Users/ben/Dropbox/Biology_projects/dros_X_A_mut/plots/gamma_anavar_',
             # pdf(paste0('/Users/ben/Desktop/gamma_anavar_',
             model, '_', pop_in, '_', algorithm_in, '_', bintype_in, '_bins.pdf',
             collapse = ''))
  par(lwd = 2,
      cex = 1.6)
  plot(x = 1:5,
       y = none$sw_gamma_1,
       main = paste0(pop_in, ' ', bintype_in, ' bins'),
       xlab = 'GC content',
       ylab = expression(gamma),
       axes = F,
       pch = 6,
       ylim = c(-3, 6))
  axis(1, at = 1:5,
       labels = GClabels,
       las = 2)
  points(x = 1:5,
         y = none$ws_gamma_1,
         pch = 7)
  points(x = 1:5,
         y = m1star$ws_gamma_1,
         pch = 4)
  points(x = 1:5,
         y = m1$ws_gamma_1,
         pch = 3)
  legend('topleft',
         bty = 'n',
         legend = c('none_gamma_SW', 'none_gamma_WS', 'm1star_gamma_WS', 'm1_gamma_WS'),
         pch = c(6, 7, 4, 3),
         cex = 0.6)
  
  axis(2)
  box()
  dev.off()
}



plot_anavar_results(pop_in = 'MD', bintype_in = 'mean', algorithm_in = 'LBFGS')
plot_anavar_results(pop_in = 'ZI69', bintype_in = 'mean', algorithm_in = 'LBFGS')
plot_anavar_results(pop_in = 'MD', bintype_in = 'species', algorithm_in = 'LBFGS')
plot_anavar_results(pop_in = 'ZI69', bintype_in = 'species', algorithm_in = 'LBFGS')
plot_anavar_results(pop_in = 'MD', bintype_in = 'diff', algorithm_in = 'LBFGS')
plot_anavar_results(pop_in = 'ZI69', bintype_in = 'diff', algorithm_in = 'LBFGS')













# eyeballing stuff in terms of text:
subset(anavar_results, subset = 
              anavar_results$pop == 'MD' &
              anavar_results$binnumber == '5' &
              anavar_results$bintype == 'diff')

# some rough plots:

plot(x = 1:5,
     y = subset(anavar_results, subset = anavar_results$model == 'none' &
                  anavar_results$pop == 'MD' &
                  anavar_results$algorithm == 'LBFGS' &
                  anavar_results$bintype == 'mean')$ws_gamma_1,
     xlab = 'bin',
     ylab = 'gamma')

plot(x = 1:5,
     y = subset(anavar_results, subset = anavar_results$model == 'none' &
                  anavar_results$pop == 'MD' &
                  anavar_results$algorithm == 'LBFGS' &
                  anavar_results$bintype == 'species')$ws_gamma_1,
     xlab = 'bin',
     ylab = 'gamma')

plot(x = 1:5,
     y = subset(anavar_results, subset = anavar_results$model == 'none' &
                  anavar_results$pop == 'MD' &
                  anavar_results$algorithm == 'LBFGS' &
                  anavar_results$bintype == 'diff')$ws_gamma_1,
     xlab = 'bin',
     ylab = 'gamma')

plot(x = 1:5,
     y = subset(anavar_results, subset = anavar_results$model == 'none' &
                  anavar_results$pop == 'ZI69' &
                  anavar_results$algorithm == 'LBFGS' &
                  anavar_results$bintype == 'mean')$ws_gamma_1,
     xlab = 'bin',
     ylab = 'gamma')

plot(x = 1:5,
     y = subset(anavar_results, subset = anavar_results$model == 'none' &
                  anavar_results$pop == 'ZI69' &
                  anavar_results$algorithm == 'LBFGS' &
                  anavar_results$bintype == 'species')$ws_gamma_1,
     xlab = 'bin',
     ylab = 'gamma')

plot(x = 1:5,
     y = subset(anavar_results, subset = anavar_results$model == 'none' &
                  anavar_results$pop == 'ZI69' &
                  anavar_results$algorithm == 'LBFGS' &
                  anavar_results$bintype == 'diff')$ws_gamma_1,
     xlab = 'bin',
     ylab = 'gamma')



plot(x = 1:5,
     y = subset(anavar_results, subset = anavar_results$model == 'none' &
                  anavar_results$pop == 'ZI69' &
                  anavar_results$algorithm == 'LBFGS' &
                  anavar_results$bintype == 'diff')$neu_e_1,
     xlab = 'bin',
     ylab = 'error')







#
