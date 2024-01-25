load(file = "rna_differential_results_Nov2023.rda")
res_B_KS1$qvalue <- qvalue(res_B_KS1$pvalue, fdr.level = 0.1)$qvalues
res_T_KS1$qvalue <- qvalue(res_T_KS1$pvalue, fdr.level = 0.1)$qvalues
res_B_KS2$qvalue <- qvalue(res_B_KS2$pvalue, fdr.level = 0.1)$qvalues
res_T_KS2$qvalue <- qvalue(res_T_KS2$pvalue, fdr.level = 0.1)$qvalues

quartz(file = "KS1_KS2_RNA_B_vs_T.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(res_T_KS1$pvalue[which(rownames(res_T_KS1) %in% 
                                      rownames(res_B_KS1)[which(res_B_KS1$qvalue <= 0.1)])], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (T cells)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 5))
lines(density(res_T_KS1$pvalue[-which(rownames(res_T_KS1) %in% 
                                       rownames(res_B_KS1)[which(res_B_KS1$qvalue <= 0.1)])], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 5))

plot(density(res_T_KS2$pvalue[which(rownames(res_T_KS2) %in% 
                                      rownames(res_B_KS2)[which(res_B_KS2$qvalue <= 0.1)])], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (T cells)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 5))
lines(density(res_T_KS2$pvalue[-which(rownames(res_T_KS2) %in% 
                                        rownames(res_B_KS2)[which(res_B_KS2$qvalue <= 0.1)])], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 5))
legend <- legend("topright", legend = c("B cells top significant hits", 
                                        "B cells other"), col = c(alpha("red", 0.57), 
                                                                  "cornflowerblue"), lwd = 2.5, 
                 bty = 'n', cex = 0.75)

dev.off()

