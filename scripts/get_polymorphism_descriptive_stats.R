setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# calculate descriptive population genetics statistics for within-population
# diversity and between-species/population divergence from the SI sites
# that we use for all the analyses

load('../data/SI_alignments.RData')
source('../scripts/R_helper_functions.R')
source('../scripts/sumstats.R')


SFS.MD21 <- lapply(alignments.intersection.MD21.SI.A, get_SFS_from_alignment, 'MD')
SFS.ZI69 <- lapply(alignments.intersection.ZI69.SI.A, get_SFS_from_alignment, 'ZI')
SFS.RG21 <- lapply(alignments.intersection.RG21.SI.A, get_SFS_from_alignment, 'RG')

SFS.MD21.sum <- sum_SFS(SFS.MD21)
SFS.ZI69.sum <- sum_SFS(SFS.ZI69)
SFS.RG21.sum <- sum_SFS(SFS.RG21)

MD21.stats <- c('pi' = pi_from_uSFS(SFS.MD21.sum),
                'theta' = theta_W_from_uSFS(SFS.MD21.sum),
                'D' = tajimas_D_from_uSFS(SFS.MD21.sum),
                'prop_single' = SFS.MD21.sum[2] / sum(SFS.MD21.sum[2:(length(SFS.MD21.sum) - 1)]))

ZI69.stats <- c('pi' = pi_from_uSFS(SFS.ZI69.sum),
                'theta' = theta_W_from_uSFS(SFS.ZI69.sum),
                'D' = tajimas_D_from_uSFS(SFS.ZI69.sum),
                'prop_single' = SFS.ZI69.sum[2] / sum(SFS.ZI69.sum[2:(length(SFS.ZI69.sum) - 1)]))

RG21.stats <- c('pi' = pi_from_uSFS(SFS.RG21.sum),
                'theta' = theta_W_from_uSFS(SFS.RG21.sum),
                'D' = tajimas_D_from_uSFS(SFS.RG21.sum),
                'prop_single' = SFS.RG21.sum[2] / sum(SFS.RG21.sum[2:(length(SFS.RG21.sum) - 1)]))
  

MD21.stats.bs <- bootstrap_sumstats(SFS.MD21)
ZI69.stats.bs <- bootstrap_sumstats(SFS.ZI69)
RG21.stats.bs <- bootstrap_sumstats(SFS.RG21)

lapply(MD21.stats.bs, FUN = quantile, c(0.025, 0.975))
lapply(ZI69.stats.bs, FUN = quantile, c(0.025, 0.975))
lapply(RG21.stats.bs, FUN = quantile, c(0.025, 0.975))

save(list = ls(pattern = '\\.stats|^SFS'),
     file = '../data/population_sumstats.RData')


