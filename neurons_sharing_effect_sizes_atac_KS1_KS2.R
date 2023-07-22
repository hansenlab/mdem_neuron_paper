###effect sign concordance between neurons, B and T cells
KS1_B_vs_neurons <- getEffectSizeConcordance(res_granges, atac_B_KS1_granges[which(atac_B_KS1_granges$qvalue <= 0.1)])
KS1_B_vs_T <- getEffectSizeConcordance(atac_T_KS1_granges, atac_B_KS1_granges[which(atac_B_KS1_granges$qvalue <= 0.1)])
KS1_B_vs_neurons_proms <- getEffectSizeConcordance(res_proms, atac_B_KS1_granges_proms[which(atac_B_KS1_granges_proms$qvalue <= 0.1)])
KS1_B_vs_T_proms <- getEffectSizeConcordance(atac_T_KS1_granges_proms, atac_B_KS1_granges_proms[which(atac_B_KS1_granges_proms$qvalue <= 0.1)])


#KS2
KS2_B_vs_neurons <- getEffectSizeConcordance(res2_granges, atac_B_KS2_granges[which(atac_B_KS2_granges$qvalue <= 0.1)])
KS2_B_vs_T <- getEffectSizeConcordance(atac_T_KS2_granges, atac_B_KS2_granges[which(atac_B_KS2_granges$qvalue <= 0.1)])
KS2_B_vs_neurons_proms <- getEffectSizeConcordance(res2_proms, 
                                        atac_B_KS2_granges_proms[which(atac_B_KS2_granges_proms$qvalue <= 0.1)])
KS2_B_vs_T_proms <- getEffectSizeConcordance(atac_T_KS2_granges_proms, 
                                  atac_B_KS2_granges_proms[which(atac_B_KS2_granges_proms$qvalue <= 0.1)])

####
#####null distributions based on sampling random locations with same balance of positive/negative logFC
#KS1
permutation_dist_KS1_B_vs_neurons <- replicate(1000, getRandomEffectSizeConcordance(res_granges, atac_B_KS1_granges, 
                                                                                    atac_B_KS1_granges$qvalue, 0.1))
permutation_dist_KS1_B_vs_T <- replicate(1000, getRandomEffectSizeConcordance(atac_T_KS1_granges, atac_B_KS1_granges, 
                                                                              atac_B_KS1_granges$qvalue, 0.1))
permutation_dist_KS1_B_vs_neurons_proms <- replicate(1000, getRandomEffectSizeConcordance(res_proms, atac_B_KS1_granges_proms, 
                                                                                          atac_B_KS1_granges_proms$qvalue, 0.1))
permutation_dist_KS1_B_vs_T_proms <- replicate(1000, getRandomEffectSizeConcordance(atac_T_KS1_granges_proms, atac_B_KS1_granges_proms, 
                                                                                    atac_B_KS1_granges_proms$qvalue, 0.1))

#KS2
permutation_dist_KS2_B_vs_neurons <- replicate(1000, getRandomEffectSizeConcordance(res2_granges, atac_B_KS2_granges, 
                                                                                    atac_B_KS2_granges$qvalue, 0.1))
permutation_dist_KS2_B_vs_T <- replicate(1000, getRandomEffectSizeConcordance(atac_T_KS2_granges, atac_B_KS2_granges, 
                                                                              atac_B_KS2_granges$qvalue, 0.1))
permutation_dist_KS2_B_vs_neurons_proms <- replicate(1000, getRandomEffectSizeConcordance(res2_proms, atac_B_KS2_granges_proms, 
                                                                       atac_B_KS2_granges_proms$qvalue, 0.1))
permutation_dist_KS2_B_vs_T_proms <- replicate(1000, getRandomEffectSizeConcordance(atac_T_KS2_granges_proms, atac_B_KS2_granges_proms, 
                                                                 atac_B_KS2_granges_proms$qvalue, 0.1))

save(permutation_dist_KS1_B_vs_neurons_proms, permutation_dist_KS1_B_vs_T_proms, 
     permutation_dist_KS2_B_vs_neurons_proms, permutation_dist_KS2_B_vs_T_proms, 
     permutation_dist_KS1_B_vs_neurons, permutation_dist_KS1_B_vs_T, 
     permutation_dist_KS2_B_vs_neurons, permutation_dist_KS2_B_vs_T, file = "neuron_permutation_distributions_effect_size_july2023.rda")


###figures
#KS1
quartz(file = "KS1_genome_wide_effect_sizes.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
permutation_dist <- permutation_dist_KS1_B_vs_neurons
plot(density(permutation_dist), 
     col = alpha("cornflowerblue"), lwd = 2.5, main = "KS1 B vs neurons (genome-wide)", bty = 'l', xlab = "% regions w/ accessibility\ndisruption in same direction", 
     font.main = 1, yaxt = 'n', xaxt = 'n', xlim = c(0.42, 1))
abline(v = KS1_B_vs_neurons, lwd = 2.5, col = alpha("red", 0.57))
axis(1, at = c(0.5, 0.75, 1))
axis(2, at = c(0, round(max(density(permutation_dist)$y))))

permutation_dist <- permutation_dist_KS1_B_vs_T
plot(density(permutation_dist), 
     col = alpha("cornflowerblue"), lwd = 2.5, main = "KS1 B vs T (genome-wide)", bty = 'l', xlab = "% regions w/ accessibility\ndisruption in same direction", 
     font.main = 1, yaxt = 'n', xaxt = 'n', xlim = c(0.42, 1))
abline(v = KS1_B_vs_T, lwd = 2.5, col = alpha("red", 0.57))
axis(1, at = c(0.5, 0.75, 1))
axis(2, at = c(0, round(max(density(permutation_dist)$y))))
legend <- legend("topright", legend = c("B cells top significant hits", 
                                        "random"), col = c(alpha("red", 0.57), 
                                                                  "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.75)
dev.off()

quartz(file = "KS1_proms_effect_sizes.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
permutation_dist <- permutation_dist_KS1_B_vs_neurons_proms
plot(density(permutation_dist), 
     col = alpha("cornflowerblue"), lwd = 2.5, main = "KS1 B vs neurons (promoters only)", bty = 'l', xlab = "% regions w/ accessibility\ndisruption in same direction", 
     font.main = 1, yaxt = 'n', xaxt = 'n', xlim = c(0.42, 1))
abline(v = KS1_B_vs_neurons_proms, lwd = 2.5, col = alpha("red", 0.57))
axis(1, at = c(0.5, 0.75, 1))
axis(2, at = c(0, round(max(density(permutation_dist)$y))))

permutation_dist <- permutation_dist_KS1_B_vs_T_proms
plot(density(permutation_dist), 
     col = alpha("cornflowerblue"), lwd = 2.5, main = "KS1 B vs T (promoters only)", bty = 'l', xlab = "% regions w/ accessibility\ndisruption in same direction", 
     font.main = 1, yaxt = 'n', xaxt = 'n', xlim = c(0.42, 1))
abline(v = KS1_B_vs_T_proms, lwd = 2.5, col = alpha("red", 0.57))
axis(1, at = c(0.5, 0.75, 1))
axis(2, at = c(0, round(max(density(permutation_dist)$y))))
legend <- legend("topright", legend = c("B cells top significant hits", 
                                        "random"), col = c(alpha("red", 0.57), 
                                                           "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.75)
dev.off()

#KS2
quartz(file = "KS2_genome_wide_effect_sizes.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
permutation_dist <- permutation_dist_KS2_B_vs_neurons
plot(density(permutation_dist), 
     col = alpha("cornflowerblue"), lwd = 2.5, main = "KS2 B vs neurons (genome-wide)", bty = 'l', xlab = "% regions w/ accessibility\ndisruption in same direction", 
     font.main = 1, yaxt = 'n', xaxt = 'n', xlim = c(0.42, 1))
abline(v = KS2_B_vs_neurons, lwd = 2.5, col = alpha("red", 0.57))
axis(1, at = c(0.5, 0.75, 1))
axis(2, at = c(0, round(max(density(permutation_dist)$y))))

permutation_dist <- permutation_dist_KS2_B_vs_T
plot(density(permutation_dist), 
     col = alpha("cornflowerblue"), lwd = 2.5, main = "KS2 B vs T (genome-wide)", bty = 'l', xlab = "% regions w/ accessibility\ndisruption in same direction", 
     font.main = 1, yaxt = 'n', xaxt = 'n', xlim = c(0.42, 1))
abline(v = KS2_B_vs_T, lwd = 2.5, col = alpha("red", 0.57))
axis(1, at = c(0.5, 0.75, 1))
axis(2, at = c(0, round(max(density(permutation_dist)$y))))
legend <- legend("topright", legend = c("B cells top significant hits", 
                                        "random"), col = c(alpha("red", 0.57), 
                                                           "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.75)
dev.off()

quartz(file = "KS2_proms_effect_sizes.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
permutation_dist <- permutation_dist_KS2_B_vs_neurons_proms
plot(density(permutation_dist), 
     col = alpha("cornflowerblue"), lwd = 2.5, main = "KS2 B vs neurons (promoters only)", bty = 'l', xlab = "% regions w/ accessibility\ndisruption in same direction", 
     font.main = 1, yaxt = 'n', xaxt = 'n', xlim = c(0.42, 1))
abline(v = KS2_B_vs_neurons_proms, lwd = 2.5, col = alpha("red", 0.57))
axis(1, at = c(0.5, 0.75, 1))
axis(2, at = c(0, round(max(density(permutation_dist)$y))))

permutation_dist <- permutation_dist_KS2_B_vs_T_proms
plot(density(permutation_dist), 
     col = alpha("cornflowerblue"), lwd = 2.5, main = "KS2 B vs T (promoters only)", bty = 'l', xlab = "% regions w/ accessibility\ndisruption in same direction", 
     font.main = 1, yaxt = 'n', xaxt = 'n', xlim = c(0.42, 1))
abline(v = KS2_B_vs_T_proms, lwd = 2.5, col = alpha("red", 0.57))
axis(1, at = c(0.5, 0.75, 1))
axis(2, at = c(0, round(max(density(permutation_dist)$y))))
legend <- legend("topright", legend = c("B cells top significant hits", 
                                        "random"), col = c(alpha("red", 0.57), 
                                                           "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.75)
dev.off()


quartz(file = "neuron_B_T_overlaps_effect_sizes.pdf", height = 2.2, width = 3.6, pointsize = 8, type = "pdf")
par(mar = c(4, 5, 1, 1) + 0.1)
plot(1, KS1_B_vs_neurons_proms/mean(permutation_dist_KS1_B_vs_neurons_proms), cex = 1.25, col = "orange", pch = 19,
     bty = 'l', main = "", xlab = "", ylab = "concordance in direction of\naccessibility disruption w/ B cells", xlim = c(0.8, 8.2), 
     ylim = c(0.9, 2), xaxt = 'n', yaxt = 'n', cex.lab = 0.75)
points(2, KS1_B_vs_T_proms/mean(permutation_dist_KS1_B_vs_T_proms), cex = 1.25, col = alpha("black", 0.62), pch = 19)
points(3, KS2_B_vs_neurons_proms/mean(permutation_dist_KS2_B_vs_neurons_proms), cex = 1.25, col = "orange", pch = 19)
points(4, KS2_B_vs_T_proms/mean(permutation_dist_KS2_B_vs_T_proms), cex = 1.25, col = alpha("black", 0.62), pch = 19)
points(5, KS1_B_vs_neurons/mean(permutation_dist_KS1_B_vs_neurons), cex = 1.25, col = "orange", pch = 19)
points(6, KS1_B_vs_T/mean(permutation_dist_KS1_B_vs_T), cex = 1.25, col = alpha("black", 0.62), pch = 19)
points(7, KS2_B_vs_neurons/mean(permutation_dist_KS2_B_vs_neurons), cex = 1.25, col = "orange", pch = 19)
points(8, KS2_B_vs_T/mean(permutation_dist_KS2_B_vs_T), cex = 1.25, col = alpha("black", 0.62), pch = 19)
abline(v = 4.5, lty = "longdash", col = rgb(0,0,0,0.7))
abline(h = 1, lty = "longdash", col = rgb(0,0,0,0.7))
legend <- legend("topleft", legend = c("neurons", "T cells"), col = c("orange", alpha("black", 0.62)), pch = 19, 
                 cex = 0.93)
axis(2, at = c(1, 2))
axis(1, at = seq(1.5, 7.5, by = 2), labels = c("KS1", "KS2", "KS1", "KS2"))
axis(1, at = c(2.5, 6.5), labels = c("promoters", "genome-wide"), col = NA, line = 1.5)
dev.off()



###effect size comparisons:
quartz(file = "KS1_KS2_magnitude_effect_sizes.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
indices_top <- which(res_granges$qvalue < 0.1)
indices_overlap <- unique(queryHits(findOverlaps(res_granges, 
                                                 atac_B_KS1_granges[which(atac_B_KS1_granges$qvalue < 0.1)])))
plot(density(-res_granges$log2FoldChange[indices_top]), 
     col = "orange", lwd = 2.5, main = "genome-wide", bty = 'l', 
     xlab = "log2(fold change) in KS1 neurons", 
     font.main = 1, yaxt = 'n', xaxt = 'n', xlim = c(-0.93, 0.93))
lines(density(-res_granges$log2FoldChange[indices_overlap]), 
      col = "deep pink", lwd = 2.5)
lines(density(-res_granges$log2FoldChange[-unique(c(indices_top, indices_overlap))]), col = "cornflowerblue", 
      lwd = 2.5)
axis(1, at = c(-0.4, 0, 0.4))
axis(2, at = c(0, round(max(density(-res_granges$log2FoldChange[indices_top])$y, 
                            density(-res_granges$log2FoldChange[indices_overlap])$y, 
                            density(-res_granges$log2FoldChange[-unique(c(indices_top, indices_overlap))])$y))))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
legend <- legend("topleft", legend = c("top significant\nneuronal peaks", 
                                       "neuronal peaks\noverlapping top B\nsignificant\npeaks", 
                                       "other neuronal peaks"), 
                 col = c("orange", "deep pink", "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.52)



indices_top <- which(res2_granges$qvalue < 0.1)
indices_overlap <- unique(queryHits(findOverlaps(res2_granges, 
                                    atac_B_KS2_granges[which(atac_B_KS2_granges$qvalue < 0.1)])))
plot(density(-res2_granges$log2FoldChange[indices_overlap]), 
     col = "deep pink", lwd = 2.5, main = "genome-wide", bty = 'l', 
     xlab = "log2(fold change) in KS2 neurons", 
     font.main = 1, yaxt = 'n', xaxt = 'n', xlim = c(-1.4, 1.4), ylim = c(0, 2))
lines(density(-res2_granges$log2FoldChange[indices_top]), 
      col = "orange", lwd = 2.5)
lines(density(-res2_granges$log2FoldChange[-unique(c(indices_top, indices_overlap))]), col = "cornflowerblue", 
      lwd = 2.5)
axis(1, at = c(-0.5, 0, 0.5))
axis(2, at = c(0, round(max(density(-res2_granges$log2FoldChange[indices_top])$y, 
                            density(-res2_granges$log2FoldChange[indices_overlap])$y, 
                            density(-res2_granges$log2FoldChange[-unique(c(indices_top, indices_overlap))])$y))))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()



#indices_top <- which(atac_T_KS1_granges$pvalue < quantile(atac_T_KS1_granges$pvalue, 0.05))
#indices_overlap <- unique(queryHits(findOverlaps(atac_T_KS1_granges, 
#                                                 atac_B_KS1_granges[which(atac_B_KS1_granges$padj < 0.1)])))
#plot(density(-atac_T_KS1_granges$log2FoldChange[indices_overlap]), 
#     col = "deep pink", lwd = 2.5, main = "KS1 T cells (genome-wide)", bty = 'l', 
#     xlab = "log2(fold change)", 
#     font.main = 1, yaxt = 'n', xaxt = 'n')
#lines(density(-atac_T_KS1_granges$log2FoldChange[indices_top]), 
#      col = "orange", lwd = 2.5)
#lines(density(-atac_T_KS1_granges$log2FoldChange[-unique(c(indices_top, indices_overlap))]), col = "cornflowerblue", 
#      lwd = 2.5)
#axis(1, at = c(-0.4, 0, 0.4))
#axis(2, at = c(0, round(max(density(-atac_T_KS1_granges$log2FoldChange[indices_top])$y, 
#                            density(-atac_T_KS1_granges$log2FoldChange[indices_overlap])$y, 
#                            density(-atac_T_KS1_granges$log2FoldChange[-unique(c(indices_top, indices_overlap))])$y))))
#abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
#legend <- legend("topleft", legend = c("top significant\nT cell peaks", 
#                                       "T cell peaks\noverlapping top B\nsignificant\npeaks", 
#                                       "other T cell peaks"), 
#                 col = c("orange", "deep pink", "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.65)




