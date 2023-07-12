all_aging_cpgs <- read_csv('~/Downloads/aging_cpgs_all_mouse.csv') #from Thompson et al. 2018 (the elastic net clock)
all_aging_cpgs$End <- all_aging_cpgs$Start
all_aging_cpgs <- makeGRangesFromDataFrame(all_aging_cpgs, keep.extra.columns = TRUE)

cpgs_in_cgi <-  all_aging_cpgs[unique(queryHits(findOverlaps(all_aging_cpgs, cpg)))]
cpgs_not_in_cgi <-  all_aging_cpgs[-unique(queryHits(findOverlaps(all_aging_cpgs, cpg)))]

###observed
res_non_proms <- res_granges[-queryHits(findOverlaps(res_granges, res_proms))]
res2_non_proms <- res2_granges[-queryHits(findOverlaps(res2_granges, res2_proms))]
atac_B_KS1_granges_non_proms <- atac_B_KS1_granges[-queryHits(findOverlaps(atac_B_KS1_granges, atac_B_KS1_granges_proms))]
atac_B_KS2_granges_non_proms <- atac_B_KS2_granges[-queryHits(findOverlaps(atac_B_KS2_granges, atac_B_KS2_granges_proms))]
atac_T_KS1_granges_non_proms <- atac_T_KS1_granges[-queryHits(findOverlaps(atac_T_KS1_granges, atac_T_KS1_granges_proms))]
atac_T_KS2_granges_non_proms <- atac_T_KS2_granges[-queryHits(findOverlaps(atac_T_KS2_granges, atac_T_KS2_granges_proms))]


aging_cpgs_in_proms_rank_KS1_neurons <- getWilcoxStat(res_proms, all_aging_cpgs)
aging_cpgs_in_proms_rank_KS2_neurons <- getWilcoxStat(res2_proms, all_aging_cpgs)
aging_cpgs_not_in_proms_rank_KS1_neurons <- getWilcoxStat(res_non_proms, all_aging_cpgs)
aging_cpgs_not_in_proms_rank_KS2_neurons <- getWilcoxStat(res2_non_proms, all_aging_cpgs)

aging_cpgs_in_proms_rank_KS1_B <- getWilcoxStat(atac_B_KS1_granges_proms, all_aging_cpgs)
aging_cpgs_in_proms_rank_KS2_B <- getWilcoxStat(atac_B_KS2_granges_proms, all_aging_cpgs)
aging_cpgs_not_in_proms_rank_KS1_B <- getWilcoxStat(atac_B_KS1_granges_non_proms, all_aging_cpgs)
aging_cpgs_not_in_proms_rank_KS2_B <- getWilcoxStat(atac_B_KS2_granges_non_proms, all_aging_cpgs)

aging_cpgs_in_proms_rank_KS1_T <- getWilcoxStat(atac_T_KS1_granges_proms, all_aging_cpgs)
aging_cpgs_in_proms_rank_KS2_T <- getWilcoxStat(atac_T_KS2_granges_proms, all_aging_cpgs)
aging_cpgs_not_in_proms_rank_KS1_T <- getWilcoxStat(atac_T_KS1_granges_non_proms, all_aging_cpgs)
aging_cpgs_not_in_proms_rank_KS2_T <- getWilcoxStat(atac_T_KS1_granges_non_proms, all_aging_cpgs)

###perm
perm_dist_aging_cpgs_in_proms_KS1_neurons <- getWilcoxNullDistribution(res_proms, all_aging_cpgs)
perm_dist_aging_cpgs_in_proms_KS2_neurons <- getWilcoxNullDistribution(res2_proms, all_aging_cpgs)
perm_dist_aging_cpgs_not_in_proms_KS1_neurons <- getWilcoxNullDistribution(res_non_proms, 
                                                                           all_aging_cpgs)
perm_dist_aging_cpgs_not_in_proms_KS2_neurons <- getWilcoxNullDistribution(res2_non_proms, 
                                                                           all_aging_cpgs)

perm_dist_aging_cpgs_in_proms_KS1_B <- getWilcoxNullDistribution(atac_B_KS1_granges_proms, all_aging_cpgs)
perm_dist_aging_cpgs_in_proms_KS2_B <- getWilcoxNullDistribution(atac_B_KS2_granges_proms, all_aging_cpgs)
perm_dist_aging_cpgs_not_in_proms_KS1_B <- getWilcoxNullDistribution(atac_B_KS1_granges_non_proms, 
                                                                     all_aging_cpgs)
perm_dist_aging_cpgs_not_in_proms_KS2_B <- getWilcoxNullDistribution(atac_B_KS2_granges_non_proms, 
                                                                     all_aging_cpgs)

perm_dist_aging_cpgs_in_proms_KS1_T <- getWilcoxNullDistribution(atac_T_KS1_granges_proms, all_aging_cpgs)
perm_dist_aging_cpgs_in_proms_KS2_T <- getWilcoxNullDistribution(atac_T_KS2_granges_proms, all_aging_cpgs)
perm_dist_aging_cpgs_not_in_proms_KS1_T <- getWilcoxNullDistribution(atac_T_KS1_granges_non_proms, 
                                                                     all_aging_cpgs)
perm_dist_aging_cpgs_not_in_proms_KS2_T <- getWilcoxNullDistribution(atac_T_KS2_granges_non_proms, 
                                                                     all_aging_cpgs)
save(perm_dist_aging_cpgs_in_proms_KS1_neurons, perm_dist_aging_cpgs_in_proms_KS2_neurons, 
     perm_dist_aging_cpgs_not_in_proms_KS1_neurons, perm_dist_aging_cpgs_not_in_proms_KS2_neurons, 
     perm_dist_aging_cpgs_in_proms_KS1_B, perm_dist_aging_cpgs_in_proms_KS2_B,
     perm_dist_aging_cpgs_not_in_proms_KS1_B, perm_dist_aging_cpgs_not_in_proms_KS2_B,
     perm_dist_aging_cpgs_in_proms_KS1_T, perm_dist_aging_cpgs_in_proms_KS2_T,
     perm_dist_aging_cpgs_not_in_proms_KS1_T, perm_dist_aging_cpgs_not_in_proms_KS2_T,
     file = "KS_aging_cpg_permutation_distributions.rda")

#plot wilcox null and observed
makeWilcoxPlot <- function(perm_dist, observed_value, main_lab){
  perm_ranks <- perm_dist
  observed_rank <- observed_value
  
  x_left_lim <- min(min(perm_ranks), observed_rank)
  x_right_lim <- max(max(perm_ranks), observed_rank)
  
  x_left_lim_label <- gsub("e\\+0","e", formatC(x_left_lim, format = "e", digits = 1))
  x_right_lim_label <- gsub("e\\+0","e", formatC(x_right_lim, format = "e", digits = 1))
  
  y_upper_lim <- max(density(perm_ranks)$y)
  y_upper_lim_label <- gsub("e\\-0","e-", formatC(y_upper_lim, format = "e", digits = 1))
  hist(perm_ranks, col = "cornflowerblue", lty = 0, 
       breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", cex.lab = 1, yaxt = 'n',
       main = main_lab, cex.main = 1, font.main = 1, 
       xlim = c(x_left_lim-0.05, x_right_lim +0.05), xaxt = 'n')
  axis(1, at = c(x_left_lim-0.05, x_right_lim +0.05), cex.axis = 1, 
       labels = c(x_left_lim_label, x_right_lim_label))
  axis(2, at = c(0, y_upper_lim), cex.axis = 1, 
       labels = c("0", y_upper_lim_label))
  abline(v = observed_rank, col = alpha("red", 0.6), lwd = 2.5)
}



###make plots
quartz(file = "aging_cpg_in_proms_ranks_KS1.pdf", width = 6, height = 2, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
makeWilcoxPlot(perm_dist_aging_cpgs_in_proms_KS1_neurons, aging_cpgs_in_proms_rank_KS1_neurons, 
               main_lab = "KS1 neurons -\nage CpGs in promoters")
legend <- legend("top", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
                 cex = 1, lty = "solid", lwd = 2.5)

makeWilcoxPlot(perm_dist_aging_cpgs_in_proms_KS1_B, aging_cpgs_in_proms_rank_KS1_B, 
               main_lab = "KS1 B -\nage CpGs in promoters")
makeWilcoxPlot(perm_dist_aging_cpgs_in_proms_KS1_T, aging_cpgs_in_proms_rank_KS1_T, 
               main_lab = "KS1 T -\nage CpGs in promoters")
dev.off()

quartz(file = "aging_cpg_in_proms_ranks_KS2.pdf", width = 6, height = 2, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
makeWilcoxPlot(perm_dist_aging_cpgs_in_proms_KS2_neurons, aging_cpgs_in_proms_rank_KS2_neurons, 
               main_lab = "KS2 neurons -\nage CpGs in promoters")
legend <- legend("top", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
                 cex = 1, lty = "solid", lwd = 2.5)

makeWilcoxPlot(perm_dist_aging_cpgs_in_proms_KS2_B, aging_cpgs_in_proms_rank_KS2_B, 
               main_lab = "KS2 B -\nage CpGs in promoters")
makeWilcoxPlot(perm_dist_aging_cpgs_in_proms_KS2_T, aging_cpgs_in_proms_rank_KS2_T, 
               main_lab = "KS2 T -\nage CpGs in promoters")
dev.off()

quartz(file = "aging_cpg_not_in_proms_ranks_KS1.pdf", width = 6, height = 2, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
makeWilcoxPlot(perm_dist_aging_cpgs_not_in_proms_KS1_neurons, aging_cpgs_not_in_proms_rank_KS1_neurons, 
               main_lab = "KS1 neurons -\nage CpGs outside promoters")
legend <- legend("top", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
                 cex = 1, lty = "solid", lwd = 2.5)

makeWilcoxPlot(perm_dist_aging_cpgs_not_in_proms_KS1_B, aging_cpgs_not_in_proms_rank_KS1_B, 
               main_lab = "KS1 B -\nage CpGs outside promoters")
makeWilcoxPlot(perm_dist_aging_cpgs_not_in_proms_KS1_T, aging_cpgs_not_in_proms_rank_KS1_T, 
               main_lab = "KS1 T -\nage CpGs outside promoters")
dev.off()

quartz(file = "aging_cpg_not_in_proms_ranks_KS2.pdf", width = 6, height = 2, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
makeWilcoxPlot(perm_dist_aging_cpgs_not_in_proms_KS2_neurons, aging_cpgs_not_in_proms_rank_KS2_neurons, 
               main_lab = "KS2 neurons -\nage CpGs outside promoters")
legend <- legend("top", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
                 cex = 1, lty = "solid", lwd = 2.5)

makeWilcoxPlot(perm_dist_aging_cpgs_not_in_proms_KS2_B, aging_cpgs_not_in_proms_rank_KS2_B, 
               main_lab = "KS2 B -\nage CpGs outside promoters")
makeWilcoxPlot(perm_dist_aging_cpgs_not_in_proms_KS2_T, aging_cpgs_not_in_proms_rank_KS2_T, 
               main_lab = "KS2 T -\nage CpGs outside promoters")
dev.off()





###alternative, splitting cpgs based on whether they overlap a CpG island or not
aging_cpgs_in_cgi_rank_KS1_neurons <- getWilcoxStat(res_granges, cpgs_in_cgi)
aging_cpgs_in_cgi_rank_KS2_neurons <- getWilcoxStat(res2_granges, cpgs_in_cgi)
aging_cpgs_not_in_cgi_rank_KS1_neurons <- getWilcoxStat(res_granges, cpgs_not_in_cgi)
aging_cpgs_not_in_cgi_rank_KS2_neurons <- getWilcoxStat(res2_granges, cpgs_not_in_cgi)

aging_cpgs_in_cgi_rank_KS1_B <- getWilcoxStat(atac_B_KS1_granges, cpgs_in_cgi)
aging_cpgs_in_cgi_rank_KS2_B <- getWilcoxStat(atac_B_KS2_granges, cpgs_in_cgi)
aging_cpgs_not_in_cgi_rank_KS1_B <- getWilcoxStat(atac_B_KS1_granges, cpgs_not_in_cgi)
aging_cpgs_not_in_cgi_rank_KS2_B <- getWilcoxStat(atac_B_KS2_granges, cpgs_not_in_cgi)

aging_cpgs_in_cgi_rank_KS1_T <- getWilcoxStat(atac_T_KS1_granges, cpgs_in_cgi)
aging_cpgs_in_cgi_rank_KS2_T <- getWilcoxStat(atac_T_KS2_granges, cpgs_in_cgi)
aging_cpgs_not_in_cgi_rank_KS1_T <- getWilcoxStat(atac_T_KS1_granges, cpgs_not_in_cgi)
aging_cpgs_not_in_cgi_rank_KS2_T <- getWilcoxStat(atac_T_KS2_granges, cpgs_not_in_cgi)

##
perm_dist_aging_cpgs_in_cgi_KS1_neurons <- getWilcoxNullDistribution(res_granges, cpgs_in_cgi)
perm_dist_aging_cpgs_in_cgi_KS2_neurons <- getWilcoxNullDistribution(res2_granges, cpgs_in_cgi)
perm_dist_aging_cpgs_not_in_cgi_KS1_neurons <- getWilcoxNullDistribution(res_granges, cpgs_not_in_cgi)
perm_dist_aging_cpgs_not_in_cgi_KS2_neurons <- getWilcoxNullDistribution(res2_granges, cpgs_not_in_cgi)

perm_dist_aging_cpgs_in_cgi_KS1_B <- getWilcoxNullDistribution(atac_B_KS1_granges, cpgs_in_cgi)
perm_dist_aging_cpgs_in_cgi_KS2_B <- getWilcoxNullDistribution(atac_B_KS2_granges, cpgs_in_cgi)
perm_dist_aging_cpgs_not_in_cgi_KS1_B <- getWilcoxNullDistribution(atac_B_KS1_granges, cpgs_not_in_cgi)
perm_dist_aging_cpgs_not_in_cgi_KS2_B <- getWilcoxNullDistribution(atac_B_KS2_granges, cpgs_not_in_cgi)

perm_dist_aging_cpgs_in_cgi_KS1_T <- getWilcoxNullDistribution(atac_T_KS1_granges, cpgs_in_cgi)
perm_dist_aging_cpgs_in_cgi_KS2_T <- getWilcoxNullDistribution(atac_T_KS2_granges, cpgs_in_cgi)
perm_dist_aging_cpgs_not_in_cgi_KS1_T <- getWilcoxNullDistribution(atac_T_KS1_granges, cpgs_not_in_cgi)
perm_dist_aging_cpgs_not_in_cgi_KS2_T <- getWilcoxNullDistribution(atac_T_KS2_granges, cpgs_not_in_cgi)






