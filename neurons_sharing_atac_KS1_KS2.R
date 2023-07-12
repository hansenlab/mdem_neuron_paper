#####observed pi0 quantifying sharing
#KS1
KS1_B_vs_neurons <- getOverlapPi0(res_granges, atac_B_KS1_granges[which(atac_B_KS1_granges$padj < 0.1)])
KS1_B_vs_T <- getOverlapPi0(atac_T_KS1_granges, atac_B_KS1_granges[which(atac_B_KS1_granges$padj < 0.1)])
KS1_B_vs_neurons_proms <- getOverlapPi0(res_proms, atac_B_KS1_granges_proms[which(atac_B_KS1_granges_proms$padj < 0.1)])
KS1_B_vs_T_proms <- getOverlapPi0(atac_T_KS1_granges_proms, atac_B_KS1_granges_proms[which(atac_B_KS1_granges_proms$padj < 0.1)])


#KS2
KS2_B_vs_neurons <- getOverlapPi0(res2_granges, atac_B_KS2_granges[which(atac_B_KS2_granges$pvalue <= 
                                                                   quantile(atac_B_KS2_granges$pvalue, 0.01))])
KS2_B_vs_T <- getOverlapPi0(atac_T_KS2_granges, atac_B_KS2_granges[which(atac_B_KS2_granges$pvalue <= 
                                                             quantile(atac_B_KS2_granges$pvalue, 0.01))])
KS2_B_vs_neurons_proms <- getOverlapPi02(res2_proms, 
                                        atac_B_KS2_granges_proms[which(atac_B_KS2_granges_proms$pvalue <= 
                                                             quantile(atac_B_KS2_granges_proms$pvalue, 0.05))])
KS2_B_vs_T_proms <- getOverlapPi02(atac_T_KS2_granges_proms, 
                                  atac_B_KS2_granges_proms[which(atac_B_KS2_granges_proms$pvalue <= 
                                                                   quantile(atac_B_KS2_granges_proms$pvalue, 0.05))])

#####null distributions based on sampling location sets of same size
#KS1
permutation_dist_KS1_B_vs_neurons <- getOverlapNullDistribution2(res_granges, atac_B_KS1_granges, 
                                                                 atac_B_KS1_granges$padj, 0.1)
permutation_dist_KS1_B_vs_T <- getOverlapNullDistribution2(atac_T_KS1_granges, atac_B_KS1_granges, 
                                                           atac_B_KS1_granges$padj, 0.1)
permutation_dist_KS1_B_vs_neurons_proms <- getOverlapNullDistribution2(res_proms, atac_B_KS1_granges_proms, 
                                                                       atac_B_KS1_granges_proms$padj, 0.1)
permutation_dist_KS1_B_vs_T_proms <- getOverlapNullDistribution2(atac_T_KS1_granges_proms, atac_B_KS1_granges_proms, 
                                                                 atac_B_KS1_granges_proms$padj, 0.1)

#KS2
permutation_dist_KS2_B_vs_neurons <- getOverlapNullDistribution2(res2_granges, atac_B_KS2_granges, 
                                                           atac_B_KS2_granges$pvalue, 
                                                           quantile(atac_B_KS2_granges$pvalue, 0.01))
permutation_dist_KS2_B_vs_T <- getOverlapNullDistribution2(atac_T_KS2_granges, atac_B_KS2_granges, 
                                                           atac_B_KS2_granges$pvalue, 
                                                           quantile(atac_B_KS2_granges$pvalue, 0.01))
permutation_dist_KS2_B_vs_neurons_proms <- getOverlapNullDistribution3(res2_proms, atac_B_KS2_granges_proms, 
                                                           atac_B_KS2_granges_proms$pvalue, 
                                                           quantile(atac_B_KS2_granges_proms$pvalue, 0.05))
permutation_dist_KS2_B_vs_T_proms <- getOverlapNullDistribution3(atac_T_KS2_granges_proms, atac_B_KS2_granges_proms, 
                                                           atac_B_KS2_granges_proms$pvalue, 
                                                           quantile(atac_B_KS2_granges_proms$pvalue, 0.05))

save(permutation_dist_KS1_B_vs_neurons_proms, permutation_dist_KS1_B_vs_T_proms, 
     permutation_dist_KS2_B_vs_neurons_proms, permutation_dist_KS2_B_vs_T_proms, 
     permutation_dist_KS1_B_vs_neurons, permutation_dist_KS1_B_vs_T, 
     permutation_dist_KS2_B_vs_neurons, permutation_dist_KS2_B_vs_T, file = "neuron_permutation_distributions_pi0.rda")

###figures
#KS1
quartz(file = "KS1_genome_wide.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(res_granges$pvalue[unique(queryHits(findOverlaps(res_granges, 
                              atac_B_KS1_granges[which(atac_B_KS1_granges$padj < 0.1)])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (genome-wide)", bty = 'l', xlab = "p-value (neurons)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 3.5))
lines(density(res_granges$pvalue[unique(queryHits(findOverlaps(res_granges, 
                                            atac_B_KS1_granges[-which(atac_B_KS1_granges$padj < 0.1)])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 3.5))

plot(density(atac_T_KS1_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS1_granges, 
                                                              atac_B_KS1_granges[which(atac_B_KS1_granges$padj < 0.1)])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (genome-wide)", bty = 'l', xlab = "p-value (T cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 3.25))
lines(density(atac_T_KS1_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS1_granges, 
                                                               atac_B_KS1_granges[-which(atac_B_KS1_granges$padj < 0.1)])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 3))
legend <- legend("topright", legend = c("B cells top significant hits", 
                                        "B cells other"), col = c(alpha("red", 0.57), 
                                                                  "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.75)

dev.off()

quartz(file = "KS1_proms.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(res_proms$pvalue[unique(queryHits(findOverlaps(res_proms, 
                                         atac_B_KS1_granges_proms[which(atac_B_KS1_granges_proms$padj < 0.1)])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (promoters)", bty = 'l', xlab = "p-value (neurons)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 5.5))
lines(density(res_proms$pvalue[unique(queryHits(findOverlaps(res_proms, 
                                         atac_B_KS1_granges_proms[-which(atac_B_KS1_granges_proms$padj < 0.1)])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 5.5))

plot(density(atac_T_KS1_granges_proms$pvalue[unique(queryHits(findOverlaps(atac_T_KS1_granges_proms, 
                                         atac_B_KS1_granges_proms[which(atac_B_KS1_granges_proms$padj < 0.1)])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (promoters)", bty = 'l', xlab = "p-value (T cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 4.25))
lines(density(atac_T_KS1_granges_proms$pvalue[unique(queryHits(findOverlaps(atac_T_KS1_granges_proms, 
                                        atac_B_KS1_granges_proms[-which(atac_B_KS1_granges_proms$padj < 0.1)])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 4))
dev.off()

#KS2
quartz(file = "KS2_genome_wide.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(res2_granges$pvalue[unique(queryHits(findOverlaps(res2_granges, 
                                                              atac_B_KS2_granges[which(atac_B_KS2_granges$pvalue < quantile(atac_B_KS2_granges$pvalue, 0.01))])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (genome-wide)", bty = 'l', xlab = "p-value (neurons)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 4))
lines(density(res2_granges$pvalue[unique(queryHits(findOverlaps(res2_granges, 
                                                               atac_B_KS2_granges[-which(atac_B_KS2_granges$pvalue < quantile(atac_B_KS2_granges$pvalue, 0.01))])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 4))

plot(density(atac_T_KS2_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS2_granges, 
                                                                     atac_B_KS2_granges[which(atac_B_KS2_granges$pvalue < quantile(atac_B_KS2_granges$pvalue, 0.01))])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (genome-wide)", bty = 'l', xlab = "p-value (T cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 4.4))
lines(density(atac_T_KS2_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS2_granges, 
                                                                      atac_B_KS2_granges[-which(atac_B_KS2_granges$pvalue < quantile(atac_B_KS2_granges$pvalue, 0.01))])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 4))
dev.off()

quartz(file = "KS2_proms.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(res2_proms$pvalue[unique(queryHits(findOverlaps(res2_proms, 
                                                    atac_B_KS2_granges_proms[which(atac_B_KS2_granges_proms$pvalue < quantile(atac_B_KS2_granges_proms$pvalue, 0.05))])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (promoters)", bty = 'l', xlab = "p-value (neurons)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 6.5))
lines(density(res2_proms$pvalue[unique(queryHits(findOverlaps(res2_proms, 
                                                    atac_B_KS2_granges_proms[-which(atac_B_KS2_granges_proms$pvalue < quantile(atac_B_KS2_granges_proms$pvalue, 0.05))])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 6.5))

plot(density(atac_T_KS2_granges_proms$pvalue[unique(queryHits(findOverlaps(atac_T_KS2_granges_proms, 
                                                    atac_B_KS2_granges_proms[which(atac_B_KS2_granges_proms$pvalue < quantile(atac_B_KS2_granges_proms$pvalue, 0.05))])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (promoters)", bty = 'l', xlab = "p-value (T cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 4))
lines(density(atac_T_KS2_granges_proms$pvalue[unique(queryHits(findOverlaps(atac_T_KS2_granges_proms, 
                                                    atac_B_KS2_granges_proms[-which(atac_B_KS2_granges_proms$pvalue < quantile(atac_B_KS2_granges_proms$pvalue, 0.05))])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 4))
dev.off()


quartz(file = "neuron_B_T_overlaps.pdf", height = 2.6, width = 3.2, pointsize = 8, type = "pdf")
par(mar = c(4, 5, 1, 1) + 0.1)
plot(1, KS1_B_vs_neurons_proms/mean(permutation_dist_KS1_B_vs_neurons_proms), cex = 1.25, col = "orange", pch = 19,
     bty = 'l', main = "", xlab = "", ylab = "accessibility disruption overlap\nw/ B cells over null expectation", xlim = c(0.8, 8.2), 
     ylim = c(0.7, 12.5), xaxt = 'n', yaxt = 'n', cex.axis = 0.8)
points(2, KS1_B_vs_T_proms/mean(permutation_dist_KS1_B_vs_T_proms), cex = 1.25, col = alpha("black", 0.62), pch = 19)
points(3, KS2_B_vs_neurons_proms/mean(permutation_dist_KS2_B_vs_neurons_proms), cex = 1.25, col = "orange", pch = 19)
points(4, KS2_B_vs_T_proms/mean(permutation_dist_KS2_B_vs_T_proms), cex = 1.25, col = alpha("black", 0.62), pch = 19)
points(5, KS1_B_vs_neurons/mean(permutation_dist_KS1_B_vs_neurons), cex = 1.25, col = "orange", pch = 19)
points(6, KS1_B_vs_T/mean(permutation_dist_KS1_B_vs_T), cex = 1.25, col = alpha("black", 0.62), pch = 19)
points(7, KS2_B_vs_neurons/mean(permutation_dist_KS2_B_vs_neurons), cex = 1.25, col = "orange", pch = 19)
points(8, KS2_B_vs_T/mean(permutation_dist_KS2_B_vs_T), cex = 1.25, col = alpha("black", 0.62), pch = 19)
abline(v = 4.5, lty = "longdash", col = rgb(0,0,0,0.7))
abline(h = 1, lty = "longdash", col = rgb(0,0,0,0.7))
legend <- legend("topright", legend = c("neurons", "T cells"), col = c("orange", alpha("black", 0.62)), pch = 19, 
                 cex = 0.93)
axis(2, at = c(1, 6, 11))
axis(1, at = seq(1.5, 7.5, by = 2), labels = c("KS1", "KS2", "KS1", "KS2"))
axis(1, at = c(2.5, 6.5), labels = c("promoters", "genome-wide"), col = NA, line = 1.5)
dev.off()















#############supplemental figure, conditioning on neurons to show that the result is symmetric
###figures
#KS1
quartz(file = "KS1_genome_wide_alt.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(atac_B_KS1_granges$pvalue[unique(queryHits(
  findOverlaps(atac_B_KS1_granges, res_granges[which(res_granges$pvalue <= quantile(res_granges$pvalue, 0.05))])))], 
  from = 0, to = 1, bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (genome-wide)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 2.25))
lines(density(atac_B_KS1_granges$pvalue[unique(queryHits(
  findOverlaps(atac_B_KS1_granges, res_granges[-which(res_granges$pvalue <= quantile(res_granges$pvalue, 0.05))])))], 
  from = 0, to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 2))
legend <- legend("topright", legend = c("neurons top significant hits", 
                                        "neurons other"), col = c(alpha("red", 0.57), 
                                                                  "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.75)


plot(density(atac_B_KS1_granges$pvalue[unique(queryHits(findOverlaps(atac_B_KS1_granges, 
                                                                     atac_T_KS1_granges[which(atac_T_KS1_granges$padj < 0.1)])))], from = 0, to = 1, 
             bw = 0.035), 
     col = "forest green", lwd = 2.5, main = "KS1 (genome-wide)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 5))
lines(density(atac_B_KS1_granges$pvalue[unique(queryHits(findOverlaps(atac_B_KS1_granges, 
                                                                      atac_T_KS1_granges[-which(atac_T_KS1_granges$padj < 0.1)])))], from = 0, 
              to = 1, bw = 0.035), col = alpha("black", 0.5), lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 5))
legend <- legend("topright", legend = c("T cells top significant hits", 
                                        "T cells other"), col = c("forest green", 
                                                                  alpha("black", 0.5)), lwd = 2.5, bty = 'n', cex = 0.75)

dev.off()

quartz(file = "KS1_proms_alt.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(atac_B_KS1_granges_proms$pvalue[
  unique(queryHits(findOverlaps(atac_B_KS1_granges_proms, 
                                res_proms[which(res_proms$pvalue <= quantile(res_proms$pvalue, 0.05))])))], 
  from = 0, to = 1, bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (promoters)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 2.2))
lines(density(atac_B_KS1_granges_proms$pvalue[
  unique(queryHits(findOverlaps(atac_B_KS1_granges_proms, 
                                res_proms[-which(res_proms$pvalue <= quantile(res_proms$pvalue, 0.05))])))], 
  from = 0, to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 2))

plot(density(atac_B_KS1_granges_proms$pvalue[unique(queryHits(findOverlaps(atac_B_KS1_granges_proms, 
                                                                           atac_T_KS1_granges_proms[which(atac_T_KS1_granges_proms$padj < 0.1)])))], from = 0, to = 1, 
             bw = 0.035), 
     col = "forest green", lwd = 2.5, main = "KS1 (promoters)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 7))
lines(density(atac_B_KS1_granges_proms$pvalue[unique(queryHits(findOverlaps(atac_B_KS1_granges_proms, 
                                                                            atac_T_KS1_granges_proms[-which(atac_T_KS1_granges_proms$padj < 0.1)])))], from = 0, 
              to = 1, bw = 0.035), col = alpha("black", 0.5), lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 7))
dev.off()

#KS2
quartz(file = "KS2_genome_wide_alt.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(atac_B_KS2_granges$pvalue[unique(queryHits(findOverlaps(atac_B_KS2_granges, 
                                                               res2_granges[which(res2_granges$pvalue <= quantile(res2_granges$pvalue , 0.05))])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (genome-wide)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 2))
lines(density(atac_B_KS2_granges$pvalue[unique(queryHits(findOverlaps(atac_B_KS2_granges, 
                                                                      res2_granges[-which(res2_granges$pvalue <= quantile(res2_granges$pvalue , 0.05))])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 2))

plot(density(atac_B_KS2_granges$pvalue[unique(queryHits(findOverlaps(atac_B_KS2_granges, 
                                                                     atac_T_KS2_granges[which(atac_T_KS2_granges$pvalue < quantile(atac_T_KS2_granges$pvalue, 0.01))])))], from = 0, to = 1, 
             bw = 0.035), 
     col = "forest green", lwd = 2.5, main = "KS2 (genome-wide)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 4))
lines(density(atac_B_KS2_granges$pvalue[unique(queryHits(findOverlaps(atac_B_KS2_granges, 
                                                                      atac_T_KS2_granges[-which(atac_T_KS2_granges$pvalue < quantile(atac_T_KS2_granges$pvalue, 0.01))])))], from = 0, 
              to = 1, bw = 0.035), col = alpha("black", 0.5), lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 4))
dev.off()

quartz(file = "KS2_proms_alt.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(atac_B_KS2_granges_proms$pvalue[unique(queryHits(findOverlaps(atac_B_KS2_granges_proms, 
                                                             res2_proms[which(res2_proms$pvalue <= quantile(res2_proms$pvalue , 0.05))])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (promoters)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 1.5))
lines(density(atac_B_KS2_granges_proms$pvalue[unique(queryHits(findOverlaps(atac_B_KS2_granges_proms, 
                                                                            res2_proms[-which(res2_proms$pvalue <= quantile(res2_proms$pvalue , 0.05))])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 1.5))

plot(density(atac_B_KS2_granges_proms$pvalue[unique(queryHits(findOverlaps(atac_B_KS2_granges_proms, 
                                                                           atac_T_KS2_granges_proms[which(atac_T_KS2_granges_proms$pvalue < quantile(atac_T_KS2_granges_proms$pvalue, 0.05))])))], from = 0, to = 1, 
             bw = 0.035), 
     col = "forest green", lwd = 2.5, main = "KS2 (promoters)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 4))
lines(density(atac_B_KS2_granges_proms$pvalue[unique(queryHits(findOverlaps(atac_B_KS2_granges_proms, 
                                                                            atac_T_KS2_granges_proms[-which(atac_T_KS2_granges_proms$pvalue < quantile(atac_T_KS2_granges_proms$pvalue, 0.05))])))], from = 0, 
              to = 1, bw = 0.035), col = alpha("black", 0.5), lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 4))
dev.off()





























