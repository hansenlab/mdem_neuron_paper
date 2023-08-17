load(file = "em_adjmatrix.rda")
###this loads all EM genes from Boukas et al. 2019, Genome Research 
###(genes are in order of coexpression status; 1:74 are highly co-expressed, 158:270 non-coexpressed and the rest are in between)

###EM genes
quartz(file = "atac_em_genes_neurons_KS1.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(density(gene_level_results_neurons_KS1$pval[which(toupper(gene_level_results_neurons_KS1$gene_names) %in% 
                                                         rownames(adjmatrix)[1:74])], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (neurons)", bty = 'l', xlab = "p-value", 
     font.main = 1, yaxt = 'n', xaxt = 'n')
lines(density(gene_level_results_neurons_KS1$pval[which(toupper(gene_level_results_neurons_KS1$gene_names) %in% 
                                                          rownames(adjmatrix)[158:270])], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
lines(density(gene_level_results_neurons_KS1$pval[-which(toupper(gene_level_results_neurons_KS1$gene_names) %in% 
                                                           rownames(adjmatrix))], from = 0, 
              to = 1, bw = 0.035), col = "black", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 9))
legend <- legend("topright", legend = c("highly co-expressed EM genes", 
                                        "non-co-expressed EM genes", 
                                        "other genes"), col = c(alpha("red", 0.57), "cornflowerblue", "black"),
                 lwd = 2.5, bty = 'n', cex = 0.57)
dev.off()

quartz(file = "atac_em_genes_neurons_KS2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(density(gene_level_results_neurons_KS2$pval[which(toupper(gene_level_results_neurons_KS2$gene_names) %in% 
                                                         rownames(adjmatrix)[1:74])], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (neurons)", bty = 'l', xlab = "p-value", 
     font.main = 1, yaxt = 'n', xaxt = 'n')
lines(density(gene_level_results_neurons_KS2$pval[which(toupper(gene_level_results_neurons_KS2$gene_names) %in% 
                                                          rownames(adjmatrix)[158:270])], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
lines(density(gene_level_results_neurons_KS2$pval[-which(toupper(gene_level_results_neurons_KS2$gene_names) %in% 
                                                           rownames(adjmatrix))], from = 0, 
              to = 1, bw = 0.035), col = "black", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 9.5))
dev.off()

quartz(file = "atac_em_genes_B.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(gene_level_results_B_KS1$pval[-which(toupper(gene_level_results_B_KS1$gene_names) %in% 
                                                    rownames(adjmatrix))], from = 0, to = 1), 
     col = "black", lwd = 2.5, main = "KS1 (B cells)", bty = 'l', xlab = "p-value", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 1.9))
lines(density(gene_level_results_B_KS1$pval[which(toupper(gene_level_results_B_KS1$gene_names) %in% 
                                                    rownames(adjmatrix)[158:270])], from = 0, 
              to = 1), col = "cornflowerblue", lwd = 2.5)
lines(density(gene_level_results_B_KS1$pval[which(toupper(gene_level_results_B_KS1$gene_names) %in% 
                                                    rownames(adjmatrix)[1:74])], from = 0, 
              to = 1), col = alpha("red", 0.57), lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 1.75))
legend <- legend("topright", legend = c("highly co-expressed EM genes", 
                                        "non-co-expressed EM genes", 
                                        "other genes"), col = c(alpha("red", 0.57), "cornflowerblue", "black"),
                 lwd = 2.5, bty = 'n', cex = 0.57)

plot(density(gene_level_results_B_KS2$pval[which(toupper(gene_level_results_B_KS2$gene_names) %in% 
                                                   rownames(adjmatrix)[1:74])], from = 0, to = 1), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (B cells)", bty = 'l', xlab = "p-value", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 1.8))
lines(density(gene_level_results_B_KS2$pval[which(toupper(gene_level_results_B_KS2$gene_names) %in% 
                                                    rownames(adjmatrix)[158:270])], from = 0, 
              to = 1), col = "cornflowerblue", lwd = 2.5)
lines(density(gene_level_results_B_KS2$pval[-which(toupper(gene_level_results_B_KS2$gene_names) %in% 
                                                     rownames(adjmatrix))], from = 0, 
              to = 1), col = "black", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 1.6))
dev.off()

quartz(file = "atac_em_genes_T.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(gene_level_results_T_KS1$pval[-which(toupper(gene_level_results_T_KS1$gene_names) %in% 
                                                    rownames(adjmatrix))], from = 0, to = 1), 
     col = "black", lwd = 2.5, main = "KS1 (T cells)", bty = 'l', xlab = "p-value", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 1.8))
lines(density(gene_level_results_T_KS1$pval[which(toupper(gene_level_results_T_KS1$gene_names) %in% 
                                                    rownames(adjmatrix)[158:270])], from = 0, 
              to = 1), col = "cornflowerblue", lwd = 2.5)
lines(density(gene_level_results_T_KS1$pval[which(toupper(gene_level_results_T_KS1$gene_names) %in% 
                                                    rownames(adjmatrix)[1:74])], from = 0, 
              to = 1), col = alpha("red", 0.57), lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 1.6))
legend <- legend("topright", legend = c("highly co-expressed EM genes", 
                                        "non-co-expressed EM genes", 
                                        "other genes"), col = c(alpha("red", 0.57), "cornflowerblue", "black"),
                 lwd = 2.5, bty = 'n', cex = 0.57)

plot(density(gene_level_results_T_KS2$pval[-which(toupper(gene_level_results_T_KS2$gene_names) %in% 
                                                    rownames(adjmatrix))], from = 0, to = 1), 
     col = "black", lwd = 2.5, main = "KS2 (T cells)", bty = 'l', xlab = "p-value", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 1.4))
lines(density(gene_level_results_T_KS2$pval[which(toupper(gene_level_results_T_KS2$gene_names) %in% 
                                                    rownames(adjmatrix)[158:270])], from = 0, 
              to = 1), col = "cornflowerblue", lwd = 2.5)
lines(density(gene_level_results_T_KS2$pval[which(toupper(gene_level_results_T_KS2$gene_names) %in% 
                                                    rownames(adjmatrix)[1:74])], from = 0, 
              to = 1), col = alpha("red", 0.57), lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 1.3))
dev.off()




###tables
overlaps <- findOverlaps(res_granges, proms_mouse)
results_proms <- res_granges[queryHits(overlaps)]
results_proms$gene_id <- proms_mouse$gene_id[subjectHits(overlaps)]
results_proms$gene_name <- proms_mouse$gene_name[subjectHits(overlaps)]

gene_level <- split(results_proms, results_proms$gene_name)
gene_level_em_KS1 <- gene_level[which(toupper(names(gene_level)) %in% rownames(adjmatrix))]
gene_level_em_KS1 <- endoapply(gene_level_em_KS1, function(xx) unique(xx))
gene_level_em_KS1_df <- annoGR2DF(unlist(gene_level_em_KS1))
gene_level_em_KS1_df <- -gene_level_em_KS1_df$log2FoldChange
#
overlaps <- findOverlaps(res2_granges, proms_mouse)
results_proms <- res2_granges[queryHits(overlaps)]
results_proms$gene_id <- proms_mouse$gene_id[subjectHits(overlaps)]
results_proms$gene_name <- proms_mouse$gene_name[subjectHits(overlaps)]

gene_level <- split(results_proms, results_proms$gene_name)
gene_level_em_KS2 <- gene_level[which(toupper(names(gene_level)) %in% rownames(adjmatrix))]
gene_level_em_KS2 <- endoapply(gene_level_em_KS2, function(xx) unique(xx))
gene_level_em_KS2_df <- annoGR2DF(unlist(gene_level_em_KS2))
gene_level_em_KS2_df <- -gene_level_em_KS2_df$log2FoldChange

###write the csv files
write_csv(res_df[order(res_df$qvalue), 
                 c("chr", "start", "end", "width", "log2FoldChange", "pvalue", "qvalue")], 
          "KS1_neurons_EM_genes.csv")



