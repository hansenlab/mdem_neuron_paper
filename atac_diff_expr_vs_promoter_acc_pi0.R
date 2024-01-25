quartz(file = "KS1_KS2_RNA_vs_ATAC.pdf", height = 2.2, width = 8.8, pointsize = 8, type = "pdf")
par(mfrow = c(1, 4))
plot(sapply(seq(500, 10000, by = 250), function(xx) 
  1 - pi0est(res_B_KS1$P.Value[which(rownames(res_B_KS1) %in% gene_level_results_B_KS1$gene_ids[
    which(rank(gene_level_results_B_KS1$pval) < xx)])], pi0.method = "bootstrap")$pi0), 
  pch = 19, col = alpha("red", 0.57), xlab = "promoter differential accessibility threshold", 
  ylab = "% differentially expressed genes", xaxt = 'n', yaxt = 'n', bty = 'l', main = "KS1 (B cells)", 
  font.main = 1, ylim = c(0.1, 0.52))
axis(1, at = c(1, length(seq(500, 10000, by = 250))), labels = c("500", "10000"))
axis(2, at = c(0, 0.1, 0.5))

plot(sapply(seq(500, 10000, by = 250), function(xx) 
  1 - pi0est(res_T_KS1$P.Value[which(rownames(res_T_KS1) %in% gene_level_results_T_KS1$gene_ids[
    which(rank(gene_level_results_T_KS1$pval) < xx)])], pi0.method = "bootstrap")$pi0), 
  pch = 19, col = alpha("red", 0.57), xaxt = 'n', yaxt = 'n', bty = 'l', main = "KS1 (T cells)", 
  font.main = 1, ylim = c(0.1, 0.52), xlab = "", ylab = "")
axis(1, at = c(1, length(seq(500, 10000, by = 250))), labels = c("500", "10000"))
axis(2, at = c(0, 0.1, 0.5))

plot(sapply(seq(500, 10000, by = 250), function(xx) 
  1 - pi0est(res_B_KS2$P.Value[which(rownames(res_B_KS2) %in% gene_level_results_B_KS2$gene_ids[
    which(rank(gene_level_results_B_KS2$pval) < xx)])], pi0.method = "bootstrap")$pi0), 
  pch = 19, col = alpha("red", 0.57), xaxt = 'n', yaxt = 'n', bty = 'l', main = "KS2 (B cells)", 
  font.main = 1, ylim = c(0.08, 0.52), xlab = "", ylab = "")
axis(1, at = c(1, length(seq(500, 10000, by = 250))), labels = c("500", "10000"))
axis(2, at = c(0, 0.1, 0.5))

plot(sapply(seq(500, 10000, by = 250), function(xx) 
  1 - pi0est(res_T_KS2$P.Value[which(rownames(res_T_KS2) %in% gene_level_results_T_KS2$gene_ids[
    which(rank(gene_level_results_T_KS2$pval) < xx)])], pi0.method = "bootstrap")$pi0), 
  pch = 19, col = alpha("red", 0.57), xaxt = 'n', yaxt = 'n', bty = 'l', main = "KS2 (T cells)", 
  font.main = 1, ylim = c(0.08, 0.52), xlab = "", ylab = "")
axis(1, at = c(1, length(seq(500, 10000, by = 250))), labels = c("500", "10000"))
axis(2, at = c(0, 0.1, 0.5))
dev.off()







