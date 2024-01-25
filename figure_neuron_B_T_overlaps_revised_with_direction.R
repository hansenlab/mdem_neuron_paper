atac_B_KS1_granges$qvalue <- qvalue(atac_B_KS1_granges$pvalue, fdr.level = 0.1)$qvalues
atac_B_KS2_granges$qvalue <- qvalue(atac_B_KS2_granges$pvalue, fdr.level = 0.1)$qvalues


###KS1: For log2FoldChange > 0
##genome-wide
KS1_B_vs_neurons_pos <- getOverlapPi0(res_granges, 
                                      atac_B_KS1_granges[which(atac_B_KS1_granges$log2FoldChange > 0 
                                                               & atac_B_KS1_granges$qvalue <= 0.1)])
KS1_B_vs_T_pos <- getOverlapPi0(atac_T_KS1_granges, 
                                atac_B_KS1_granges[which(atac_B_KS1_granges$log2FoldChange > 0 
                                                         & atac_B_KS1_granges$qvalue <= 0.1)])

#null
B_granges <- atac_B_KS1_granges[which(atac_B_KS1_granges$log2FoldChange > 0)]
permutation_dist_KS1_B_vs_neurons_pos <- getOverlapNullDistribution2(res_granges, B_granges, 
                                                                     B_granges$qvalue, 0.1)
permutation_dist_KS1_B_vs_T_pos <- getOverlapNullDistribution2(atac_T_KS1_granges, B_granges, 
                                                               B_granges$qvalue, 0.1)


##promoters
KS1_B_vs_neurons_proms_pos <- getOverlapPi02(res_proms, atac_B_KS1_granges_proms[
  which(atac_B_KS1_granges_proms$log2FoldChange > 0 & atac_B_KS1_granges_proms$qvalue <= 0.1)], 
  lambda_value = 0.5)
KS1_B_vs_T_proms_pos <- getOverlapPi02(atac_T_KS1_granges_proms, atac_B_KS1_granges_proms[
  which(atac_B_KS1_granges_proms$log2FoldChange > 0 & atac_B_KS1_granges_proms$qvalue <= 0.1)], 
  lambda_value = 0.5)

#null
B_granges_proms <- atac_B_KS1_granges_proms[which(atac_B_KS1_granges_proms$log2FoldChange > 0)]
permutation_dist_KS1_B_vs_neurons_proms_pos <- getOverlapNullDistribution3(res_proms, B_granges_proms, 
                                                                           B_granges_proms$qvalue, 0.1, 
                                                                           lambda_value = 0.5)
permutation_dist_KS1_B_vs_T_proms_pos <- getOverlapNullDistribution3(atac_T_KS1_granges_proms, B_granges_proms, 
                                                                     B_granges_proms$qvalue, 0.1, 
                                                                     lambda_value = 0.5)

###KS1: For log2FoldChange < 0 
##genome-wide
KS1_B_vs_neurons_neg <- getOverlapPi0(res_granges, 
                                      atac_B_KS1_granges[which(atac_B_KS1_granges$log2FoldChange < 0 
                                                               & atac_B_KS1_granges$qvalue <= 0.1)])
KS1_B_vs_T_neg <- getOverlapPi0(atac_T_KS1_granges, 
                                atac_B_KS1_granges[which(atac_B_KS1_granges$log2FoldChange < 0 
                                                         & atac_B_KS1_granges$qvalue <= 0.1)])

#null
B_granges <- atac_B_KS1_granges[which(atac_B_KS1_granges$log2FoldChange < 0)]
permutation_dist_KS1_B_vs_neurons_neg <- getOverlapNullDistribution2(res_granges, B_granges, 
                                                                     B_granges$qvalue, 0.1)
permutation_dist_KS1_B_vs_T_neg <- getOverlapNullDistribution2(atac_T_KS1_granges, B_granges, 
                                                               B_granges$qvalue, 0.1)

##promoters
KS1_B_vs_neurons_proms_neg <- getOverlapPi02(res_proms, atac_B_KS1_granges_proms[
  which(atac_B_KS1_granges_proms$log2FoldChange < 0 & atac_B_KS1_granges_proms$qvalue <= 0.1)], 
  lambda_value = 0.5)
KS1_B_vs_T_proms_neg <- getOverlapPi02(atac_T_KS1_granges_proms, atac_B_KS1_granges_proms[
  which(atac_B_KS1_granges_proms$log2FoldChange < 0 & atac_B_KS1_granges_proms$qvalue <= 0.1)], 
  lambda_value = 0.5)

#null
B_granges_proms <- atac_B_KS1_granges_proms[which(atac_B_KS1_granges_proms$log2FoldChange < 0)]
permutation_dist_KS1_B_vs_neurons_proms_neg <- getOverlapNullDistribution3(res_proms, B_granges_proms, 
                                                                           B_granges_proms$qvalue, 0.1, 
                                                                           lambda_value = 0.5)
permutation_dist_KS1_B_vs_T_proms_neg <- getOverlapNullDistribution3(atac_T_KS1_granges_proms, B_granges_proms, 
                                                                     B_granges_proms$qvalue, 0.1, 
                                                                     lambda_value = 0.5)

###KS2: For log2FoldChange > 0 
##genome-wide
KS2_B_vs_neurons_pos <- getOverlapPi02(res2_granges, 
                                       atac_B_KS2_granges[which(atac_B_KS2_granges$log2FoldChange > 0 
                                                                & atac_B_KS2_granges$qvalue <= 0.1)], 
                                       lambda_value = 0.4)
KS2_B_vs_T_pos <- getOverlapPi02(atac_T_KS2_granges, 
                                 atac_B_KS2_granges[which(atac_B_KS2_granges$log2FoldChange > 0 
                                                          & atac_B_KS2_granges$qvalue <= 0.1)],
                                 lambda_value = 0.5)

#null
B_granges <- atac_B_KS2_granges[which(atac_B_KS2_granges$log2FoldChange > 0)]
permutation_dist_KS2_B_vs_neurons_pos <- getOverlapNullDistribution3(res2_granges, B_granges, 
                                                                     B_granges$qvalue, 0.1, lambda_value = 0.4)
permutation_dist_KS2_B_vs_T_pos <- getOverlapNullDistribution3(atac_T_KS2_granges, B_granges, 
                                                               B_granges$qvalue, 0.1, lambda_value = 0.5)


##promoters (this is not plotted because there's not enough promoters with significant disruption and logFC > 0)
KS2_B_vs_neurons_proms_pos <- getOverlapPi02(res2_proms, atac_B_KS2_granges_proms[
  which(atac_B_KS2_granges_proms$log2FoldChange > 0 & atac_B_KS2_granges_proms$qvalue <= 0.1)], 
  lambda_value = 0.5)
KS2_B_vs_T_proms_pos <- getOverlapPi02(atac_T_KS2_granges_proms, atac_B_KS2_granges_proms[
  which(atac_B_KS2_granges_proms$log2FoldChange > 0 & atac_B_KS2_granges_proms$qvalue <= 0.1)], 
  lambda_value = 0.5)

#null
B_granges_proms <- atac_B_KS2_granges_proms[which(atac_B_KS2_granges_proms$log2FoldChange > 0)]
permutation_dist_KS2_B_vs_neurons_proms_pos <- getOverlapNullDistribution3(res2_proms, B_granges_proms, 
                                                                           B_granges_proms$qvalue, 0.1, 
                                                                           lambda_value = 0.5)
permutation_dist_KS2_B_vs_T_proms_pos <- getOverlapNullDistribution3(atac_T_KS2_granges_proms, B_granges_proms, 
                                                                     B_granges_proms$qvalue, 0.1, 
                                                                     lambda_value = 0.5)


###KS2: For log2FoldChange < 0 
##genome-wide
KS2_B_vs_neurons_neg <- getOverlapPi02(res2_granges, 
                                       atac_B_KS2_granges[which(atac_B_KS2_granges$log2FoldChange < 0 
                                                                & atac_B_KS2_granges$qvalue <= 0.1)], 
                                       lambda_value = 0.5)
KS2_B_vs_T_neg <- getOverlapPi02(atac_T_KS2_granges, 
                                 atac_B_KS2_granges[which(atac_B_KS2_granges$log2FoldChange < 0 
                                                          & atac_B_KS2_granges$qvalue <= 0.1)], 
                                 lambda_value = 0.5)

#null
B_granges <- atac_B_KS2_granges[which(atac_B_KS2_granges$log2FoldChange < 0)]
permutation_dist_KS2_B_vs_neurons_neg <- getOverlapNullDistribution3(res2_granges, B_granges, 
                                                                     B_granges$qvalue, 0.1, 
                                                                     lambda_value = 0.5)
permutation_dist_KS2_B_vs_T_neg <- getOverlapNullDistribution3(atac_T_KS2_granges, B_granges, 
                                                               B_granges$qvalue, 0.1, 
                                                               lambda_value = 0.5)

##promoters
KS2_B_vs_neurons_proms_neg <- getOverlapPi02(res2_proms, atac_B_KS2_granges_proms[
  which(atac_B_KS2_granges_proms$log2FoldChange < 0 & atac_B_KS2_granges_proms$qvalue <= 0.1)], 
  lambda_value = 0.5)
KS2_B_vs_T_proms_neg <- getOverlapPi02(atac_T_KS2_granges_proms, atac_B_KS2_granges_proms[
  which(atac_B_KS2_granges_proms$log2FoldChange < 0 & atac_B_KS2_granges_proms$qvalue <= 0.1)], 
  lambda_value = 0.5)

#null
B_granges_proms <- atac_B_KS2_granges_proms[which(atac_B_KS2_granges_proms$log2FoldChange < 0)]
permutation_dist_KS2_B_vs_neurons_proms_neg <- getOverlapNullDistribution3(res2_proms, B_granges_proms, 
                                                                           B_granges_proms$qvalue, 0.1, 
                                                                           lambda_value = 0.5)
permutation_dist_KS2_B_vs_T_proms_neg <- getOverlapNullDistribution3(atac_T_KS2_granges_proms, B_granges_proms, 
                                                                     B_granges_proms$qvalue, 0.1, 
                                                                     lambda_value = 0.5)


save(permutation_dist_KS1_B_vs_neurons_proms_pos, permutation_dist_KS1_B_vs_neurons_proms_neg,
     permutation_dist_KS1_B_vs_T_proms_pos, permutation_dist_KS1_B_vs_T_proms_neg,
     permutation_dist_KS2_B_vs_neurons_proms_neg,
     permutation_dist_KS2_B_vs_T_proms_neg,
     permutation_dist_KS1_B_vs_neurons_pos, permutation_dist_KS1_B_vs_neurons_neg,
     permutation_dist_KS1_B_vs_T_pos, permutation_dist_KS1_B_vs_T_neg,
     permutation_dist_KS2_B_vs_neurons_pos, permutation_dist_KS2_B_vs_neurons_neg,
     permutation_dist_KS2_B_vs_T_pos, permutation_dist_KS2_B_vs_T_neg, 
     file = "neuron_permutation_distributions_pi0_pos_and_neg_dec2023.rda")


###figure
quartz(file = "neuron_B_T_overlaps.pdf", height = 2.2, width = 3.6, pointsize = 8, type = "pdf")
par(mar = c(4, 5, 1, 1) + 0.1)
plot(0.75, KS1_B_vs_neurons_proms_pos/mean(permutation_dist_KS1_B_vs_neurons_proms_pos), cex = 1.25, col = "orange", pch = 19,
     bty = 'l', main = "", xlab = "", ylab = "accessibility disruption\noverlap w/ B cells", xlim = c(0.8, 8.5), 
     ylim = c(0, 5.9), xaxt = 'n', yaxt = 'n', cex.axis = 0.8)

points(1.25, KS1_B_vs_neurons_proms_neg/mean(permutation_dist_KS1_B_vs_neurons_proms_neg), cex = 1.25, col = "orange", pch = 15)

points(1.75, KS1_B_vs_T_proms_pos/mean(permutation_dist_KS1_B_vs_T_proms_pos), cex = 1.25, col = alpha("black", 0.62), pch = 19)
points(2.25, KS1_B_vs_T_proms_neg/mean(permutation_dist_KS1_B_vs_T_proms_neg), cex = 1.25, col = alpha("black", 0.62), pch = 15)

#points(2.75, KS2_B_vs_neurons_proms_pos/mean(permutation_dist_KS2_B_vs_neurons_proms_pos), cex = 1.25, col = "orange", pch = 19)
points(3, KS2_B_vs_neurons_proms_neg/mean(permutation_dist_KS2_B_vs_neurons_proms_neg), cex = 1.25, col = "orange", pch = 15)

#points(3.75, KS2_B_vs_T_proms_pos/mean(permutation_dist_KS2_B_vs_T_proms_pos), cex = 1.25, col = alpha("black", 0.62), pch = 19)
points(4, KS2_B_vs_T_proms_neg/mean(permutation_dist_KS2_B_vs_T_proms_neg), cex = 1.25, col = alpha("black", 0.62), pch = 15)

points(4.75, KS1_B_vs_neurons_pos/mean(permutation_dist_KS1_B_vs_neurons_pos), cex = 1.25, col = "orange", pch = 19)
points(5.25, KS1_B_vs_neurons_neg/mean(permutation_dist_KS1_B_vs_neurons_neg), cex = 1.25, col = "orange", pch = 15)

points(5.75, KS1_B_vs_T_pos/mean(permutation_dist_KS1_B_vs_T_pos), cex = 1.25, col = alpha("black", 0.62), pch = 19)
points(6.25, KS1_B_vs_T_neg/mean(permutation_dist_KS1_B_vs_T_neg), cex = 1.25, col = alpha("black", 0.62), pch = 15)

points(6.75, KS2_B_vs_neurons_pos/mean(permutation_dist_KS2_B_vs_neurons_pos), cex = 1.25, col = "orange", pch = 19)
points(7.25, KS2_B_vs_neurons_neg/mean(permutation_dist_KS2_B_vs_neurons_neg), cex = 1.25, col = "orange", pch = 15)

points(7.75, KS2_B_vs_T_pos/mean(permutation_dist_KS2_B_vs_T_pos), cex = 1.25, col = alpha("black", 0.62), pch = 19)
points(8.25, KS2_B_vs_T_neg/mean(permutation_dist_KS2_B_vs_T_neg), cex = 1.25, col = alpha("black", 0.62), pch = 15)

abline(v = c(2.5, 4.5, 6.5), lty = "longdash", col = rgb(0,0,0,0.7))
abline(h = 1, lty = "longdash", col = rgb(0,0,0,0.7))

#legend("topright", 
#       legend = c("Neurons log2(FC) > 0", "Neurons log2(FC) < 0", "T cells log2(FC) > 0", "T cells log2(FC) < 0"), 
#       col = c("orange", "orange", alpha("black", 0.62), alpha("black", 0.62)), 
#       pch = c(19, 15, 19, 15), bty = 'n',
#       cex = 0.93)

axis(2, at = c(1, 3, 5))
axis(1, at = seq(1.5, 7.5, by = 2), labels = c("KS1", "KS2", "KS1", "KS2"))
axis(1, at = c(2.5, 6.5), labels = c("promoters", "genome-wide"), col = NA, line = 1.5)
dev.off()


