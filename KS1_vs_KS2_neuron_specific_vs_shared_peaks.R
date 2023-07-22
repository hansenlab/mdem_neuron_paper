overlaps <- findOverlaps(res_granges, atac_B_KS1_granges)
res_unique <- res_granges[-unique(queryHits(overlaps))]
overlaps2 <- findOverlaps(res_unique, atac_T_KS1_granges)
res_unique <- res_unique[-unique(queryHits(overlaps2))]

overlaps <- findOverlaps(res2_granges, atac_B_KS2_granges)
res2_unique <- res2_granges[-unique(queryHits(overlaps))]
overlaps2 <- findOverlaps(res2_unique, atac_T_KS2_granges)
res2_unique <- res2_unique[-unique(queryHits(overlaps2))]


quartz(file = "KS1_KS2_unique_vs_shared.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(res_granges$pvalue[-unique(queryHits(findOverlaps(res_granges, res_unique)))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (genome-wide)", bty = 'l', xlab = "p-value", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 3))
lines(density(res_unique$pvalue, from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 3))
legend <- legend("topright", legend = c("ATAC peaks\nshared with B or T cells", 
                                        "ATAC peaks unique to neurons"), col = c(alpha("red", 0.57), 
                                                                     "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.75)

plot(density(res2_granges$pvalue[-unique(queryHits(findOverlaps(res2_granges, res2_unique)))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (genome-wide)", bty = 'l', xlab = "p-value", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 3.5))
lines(density(res2_unique$pvalue, from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 3.5))

dev.off()
