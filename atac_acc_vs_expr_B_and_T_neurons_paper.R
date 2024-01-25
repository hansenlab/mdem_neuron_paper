##########
KS1_top_promoters_B <- gene_level_results_B_KS1$gene_ids[
  which(gene_level_results_B_KS1$pval < quantile(gene_level_results_B_KS1$pval, 0.1))]
KS2_top_promoters_B <- gene_level_results_B_KS2$gene_ids[
  which(gene_level_results_B_KS2$pval < quantile(gene_level_results_B_KS2$pval, 0.1))]
KS1_top_promoters_T <- gene_level_results_T_KS1$gene_ids[
  which(gene_level_results_T_KS1$pval < quantile(gene_level_results_T_KS1$pval, 0.1))]
KS2_top_promoters_T <- gene_level_results_T_KS2$gene_ids[
  which(gene_level_results_T_KS2$pval < quantile(gene_level_results_T_KS2$pval, 0.1))]

#
top_proms_KS1_de_genes_B <- qvalue(res_B_KS1$P.Value[which(rownames(res_B_KS1) %in% KS1_top_promoters_B)], pi0.method = "bootstrap", fdr.level = 0.1)
top_proms_KS1_de_genes_B <- res_B_KS1[which(rownames(res_B_KS1) %in% KS1_top_promoters_B)[which(
  top_proms_KS1_de_genes_B$significant == TRUE)], ]

top_proms_KS2_de_genes_B <- qvalue(res_B_KS2$P.Value[which(rownames(res_B_KS2) %in% KS2_top_promoters_B)], pi0.method = "bootstrap", fdr.level = 0.1)
top_proms_KS2_de_genes_B <- res_B_KS2[which(rownames(res_B_KS2) %in% KS2_top_promoters_B)[which(
  top_proms_KS2_de_genes_B$significant == TRUE)], ]

top_proms_KS1_de_genes_T <- qvalue(res_T_KS1$P.Value[which(rownames(res_T_KS1) %in% KS1_top_promoters_T)], pi0.method = "bootstrap", fdr.level = 0.1)
top_proms_KS1_de_genes_T <- res_T_KS1[which(rownames(res_T_KS1) %in% KS1_top_promoters_T)[which(
  top_proms_KS1_de_genes_T$significant == TRUE)], ]

top_proms_KS2_de_genes_T <- qvalue(res_T_KS2$P.Value[which(rownames(res_T_KS2) %in% KS2_top_promoters_T)], pi0.method = "bootstrap", fdr.level = 0.1)
top_proms_KS2_de_genes_T <- res_T_KS2[which(rownames(res_T_KS2) %in% KS2_top_promoters_T)[which(
  top_proms_KS2_de_genes_T$significant == TRUE)], ]

#in the following lines, the `gene_level_results` objects are from the "atac_pathways_neurons_updated.R" script
proms_logFC_KS1_B <- sapply(rownames(top_proms_KS1_de_genes_B), function(xx) 
  -gene_level_results_B_KS1$logFC[which(gene_level_results_B_KS1$gene_ids == xx)]) 

proms_logFC_KS2_B <- sapply(rownames(top_proms_KS2_de_genes_B), function(xx) 
  -gene_level_results_B_KS2$logFC[which(gene_level_results_B_KS2$gene_ids == xx)])

proms_logFC_KS1_T <- sapply(rownames(top_proms_KS1_de_genes_T), function(xx) 
  -gene_level_results_T_KS1$logFC[which(gene_level_results_T_KS1$gene_ids == xx)])

proms_logFC_KS2_T <- sapply(rownames(top_proms_KS2_de_genes_T), function(xx) 
  -gene_level_results_T_KS2$logFC[which(gene_level_results_T_KS2$gene_ids == xx)])

#
pos_KS1_top <- -gene_level_results_neurons_KS1$logFC[which(gene_level_results_neurons_KS1$logFC < 0 & 
          gene_level_results_neurons_KS1$pval < quantile(gene_level_results_neurons_KS1$pval, 0.1))]
neuron_threshold_pos_KS1 <- quantile(pos_KS1_top, 0.1)

pos_KS2_top <- -gene_level_results_neurons_KS2$logFC[which(gene_level_results_neurons_KS2$logFC < 0 & 
          gene_level_results_neurons_KS2$pval < quantile(gene_level_results_neurons_KS2$pval, 0.1))]
neuron_threshold_pos_KS2 <- quantile(pos_KS2_top, 0.1)

neg_KS2_top <- -gene_level_results_neurons_KS2$logFC[
  which(gene_level_results_neurons_KS2$logFC > 0 & 
          gene_level_results_neurons_KS2$pval < quantile(gene_level_results_neurons_KS2$pval, 0.1))]
neuron_threshold_neg_KS2 <- quantile(neg_KS2_top, 0.9)
#

quartz(file = "KS1_KS2_RNA_vs_ATAC_B_scatterplot.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
plot(proms_logFC_KS1_B, 
     -top_proms_KS1_de_genes_B$log2FoldChange, 
     pch = 19, col = alpha("red", 0.57), 
     xlab = "promoter accessibility logFC", ylab = "gene expression logFC", main = "KS1", bty = 'l', yaxt = 'n', xaxt = 'n', 
     font.main = 1)
abline(h = 0, lty = "longdash", col = rgb(0,0,0,0.7))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
abline(v = c(neuron_threshold_pos_KS1), col = rgb(0,0,0,0.7))
axis(2, at = c(-2, 0, 6))
axis(1, at = c(-0.75, 0, 1.2))

plot(proms_logFC_KS2_B, 
     -top_proms_KS2_de_genes_B$log2FoldChange, 
     pch = 19, col = alpha("red", 0.57), 
     xlab = "promoter accessibility logFC", ylab = "gene expression logFC", main = "KS2", bty = 'l', yaxt = 'n', xaxt = 'n', 
     font.main = 1, xlim = c(-0.9, 1.1)) #set xlim to ensure we depict the vertical line corresponding to neurons
abline(h = 0, lty = "longdash", col = rgb(0,0,0,0.7))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(-1, 0, 1.5))
axis(1, at = c(-0.45, 0, 0.9))
abline(v = c(neuron_threshold_pos_KS2, neuron_threshold_neg_KS2), col = rgb(0,0,0,0.7))
dev.off()

quartz(file = "KS1_KS2_RNA_vs_ATAC_T_scatterplot.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
plot(proms_logFC_KS1_T, 
     -top_proms_KS1_de_genes_T$log2FoldChange, 
     pch = 19, col = alpha("red", 0.57), 
     xlab = "promoter accessibility logFC", ylab = "gene expression logFC", main = "KS1 (T cells)", bty = 'l', yaxt = 'n', xaxt = 'n', 
     font.main = 1)
abline(h = 0, lty = "longdash", col = rgb(0,0,0,0.7))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
abline(v = c(neuron_threshold_pos_KS1), col = rgb(0,0,0,0.7))
axis(2, at = c(-2, 0, 2))
axis(1, at = c(-0.75, 0, 1.5))

plot(proms_logFC_KS2_T, 
     -top_proms_KS2_de_genes_T$log2FoldChange, 
     pch = 19, col = alpha("red", 0.57), 
     xlab = "promoter accessibility logFC", ylab = "gene expression logFC", main = "KS2 (T cells)", bty = 'l', yaxt = 'n', xaxt = 'n', 
     font.main = 1)
abline(h = 0, lty = "longdash", col = rgb(0,0,0,0.7))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(-2, 0, 2))
axis(1, at = c(-1, 0, 1))
abline(v = c(neuron_threshold_pos_KS2, neuron_threshold_neg_KS2), col = rgb(0,0,0,0.7))
dev.off()


###
concordance_percent <- c(prop.table(table(as.factor(proms_logFC_KS1_B * (-top_proms_KS1_de_genes_B$log2FoldChange) > 0)))[2], 
                         prop.table(table(as.factor(proms_logFC_KS2_B * (-top_proms_KS2_de_genes_B$log2FoldChange) > 0)))[2], 
                         prop.table(table(as.factor(proms_logFC_KS1_T * (-top_proms_KS1_de_genes_T$log2FoldChange) > 0)))[2],
                         prop.table(table(as.factor(proms_logFC_KS2_T * (-top_proms_KS2_de_genes_T$log2FoldChange) > 0)))[2])


# Create a matrix from the data
data_mat <- matrix(concordance_percent, nrow = 2, byrow = TRUE)

# Assign row and column names to the matrix
dimnames(data_mat) <- list(c("B cells", "T cells"), 
                           c("KS1", "KS2"))

# Create a bar plot
quartz(file = "KS1_KS2_RNA_vs_ATAC_concordance.pdf", height = 2.2, width = 1.8, pointsize = 8, type = "pdf")
barplot(data_mat, beside=TRUE, legend.text=TRUE, args.legend = list(x = "topright", bty = 'n', border = 'NA'), 
        xlab = "", ylab = "% promoter-gene pairs with concordant disruption", col = c("palevioletred1", "gray48"), 
        border = "NA", yaxt = 'n')
axis(2, at = c(0, 0.9), labels = c("0", "0.9"))
dev.off()
