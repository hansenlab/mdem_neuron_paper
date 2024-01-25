dist_1 <- replicate(54000, rnorm(75, 0, 1))
dist_2 <- replicate(6000, rnorm(75, 0.5, 1))

matrix_1 <- dist_1[, 1:10000]
matrix_2 <- cbind(dist_2[, 1:2000], dist_1[, 10001:18000])

matrix_3 <- dist_1[, 18001:28000]
matrix_4 <- cbind(dist_2[, 2001:3400], dist_1[, 28001:36000], dist_2[, 3401:4000])

matrix_5 <- dist_1[, 36001:46000]
matrix_6 <- cbind(dist_2[, 4001:4460], dist_1[, 46001:53000], 
                  dist_2[, 4461:6000], dist_1[, 53001:54000])

library(matrixTests)
pvals_1 <- col_t_welch(matrix_1, matrix_2)$pvalue
adj_pvals_2 <- qvalue(pvals_1, fdr.level = 0.1)$qvalue
pvals_2 <- col_t_welch(matrix_3, matrix_4)$pvalue
adj_pvals_2 <- qvalue(pvals_2, fdr.level = 0.1)$qvalue
pvals_3 <- col_t_welch(matrix_5, matrix_6)$pvalue
adj_pvals_3 <- qvalue(pvals_3, fdr.level = 0.1)$qvalue

quartz(file = "conceptual_figure_neurons_mdem_density_1.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(density(pvals_1[which(adj_pvals_2 <= 0.1)], from = 0, to = 1), 
     bty = 'l', lwd = 2.5, col = "forest green", main = "", xlab = "pvalues (cell type A)", ylim = c(0, 4.2),
     xaxt = 'n', yaxt = 'n')
lines(density(pvals_1[which(adj_pvals_1 > 0.1)], from = 0, to = 1), col = 1, 
      lwd = 2)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 4))
legend <- legend("topright", legend = c("differential in cell type B", "null in cell type B"), 
                 lwd = 2, col = c("forest green", 1), bty = 'n', cex = 0.84)
dev.off()

quartz(file = "conceptual_figure_neurons_mdem_density_2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(density(pvals_1[which(adj_pvals_3 > 0.1)], from = 0, to = 1), 
     bty = 'l', lwd = 2.5, col = 1, main = "", xlab = "pvalues (cell type A)", ylim = c(0, 2), 
     xaxt = 'n', yaxt = 'n')
lines(density(pvals_1[which(adj_pvals_3 <= 0.1)], from = 0, to = 1), col = "forest green", 
      lwd = 2)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 2))
dev.off()
