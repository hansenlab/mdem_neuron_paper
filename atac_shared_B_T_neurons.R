###find peaks disrupted in all three cell types
#KS1
qobj_KS1_neurons_B <- qvalue(p = res_granges$pvalue[unique(queryHits(findOverlaps(res_granges, 
                                                                                  atac_B_KS1_granges[which(atac_B_KS1_granges$qvalue < 0.1)])))], 
                             fdr.level = 0.1, pi0.method = "bootstrap")
KS1_neurons_B <- res_granges[unique(queryHits(findOverlaps(res_granges, 
                                                           atac_B_KS1_granges[which(atac_B_KS1_granges$qvalue < 0.1)])))][
                                                             which(qobj_KS1_neurons_B$significant == TRUE)]


unique_B <- atac_B_KS1_granges[-unique(queryHits(findOverlaps(atac_B_KS1_granges, KS1_neurons_B)))]


#
qobj_KS1_neurons_B_T <- qvalue(p = atac_T_KS1_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS1_granges, 
                                                                                           KS1_neurons_B)))], 
                               fdr.level = 0.1, pi0.method = "bootstrap")
KS1_neurons_B_T <- atac_T_KS1_granges[unique(queryHits(findOverlaps(atac_T_KS1_granges, 
                                                                    KS1_neurons_B)))][which(qobj_KS1_neurons_B_T$significant == TRUE)]


shared_genes_df <- getGeneLevelResults(KS1_neurons_B_T)
shared_genes_df <- shared_genes_df[, 1:2] #keeps only gene ids and gene name
all_genes <- unique(Reduce(intersect, list(gene_level_results_neurons_KS1$gene_ids, 
                                                  gene_level_results_B_KS1$gene_ids, 
                                                  gene_level_results_T_KS1$gene_ids)))



#KS2
qobj_KS2_neurons_B <- qvalue(p = res2_granges$pvalue[unique(queryHits(findOverlaps(res2_granges, 
                                                                                  atac_B_KS2_granges[which(atac_B_KS2_granges$pvalue <= 
                                                                                                             quantile(atac_B_KS2_granges$pvalue, 0.01))])))], 
                             fdr.level = 0.1, pi0.method = "bootstrap")
KS2_neurons_B <- res2_granges[unique(queryHits(findOverlaps(res2_granges, 
                                                           atac_B_KS2_granges[which(atac_B_KS2_granges$pvalue <= 
                                                                                      quantile(atac_B_KS2_granges$pvalue, 0.01))])))][
                                                             which(qobj_KS2_neurons_B$significant == TRUE)]


unique_B_2 <- atac_B_KS2_granges[-unique(queryHits(findOverlaps(atac_B_KS2_granges, KS2_neurons_B)))]


#
qobj_KS2_neurons_B_T <- qvalue(p = atac_T_KS2_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS2_granges, 
                                                                                           KS2_neurons_B)))], 
                               fdr.level = 0.1, pi0.method = "bootstrap")
KS2_neurons_B_T <- atac_T_KS2_granges[unique(queryHits(findOverlaps(atac_T_KS2_granges, 
                               KS2_neurons_B)))][which(qobj_KS2_neurons_B_T$significant == TRUE)]

shared_genes2_df <- getGeneLevelResults(KS2_neurons_B_T)
shared_genes2_df <- shared_genes2_df[, 1:2] #keeps only gene ids and gene name
all_genes2 <- unique(Reduce(intersect, list(gene_level_results_neurons_KS2$gene_ids, 
                                                  gene_level_results_B_KS2$gene_ids, 
                                                  gene_level_results_T_KS2$gene_ids)))

shared_genes <- shared_genes_df$gene_ids
shared_genes2 <- shared_genes2_df$gene_ids
n11 <- length(intersect(shared_genes, shared_genes2))
n12 <- length(shared_genes[-which(shared_genes %in% shared_genes2)])
n21 <- length(shared_genes2[-which(shared_genes2 %in% shared_genes)])
n22 <- length(c(all_genes[-which(all_genes %in% shared_genes | all_genes %in% shared_genes2)], 
                all_genes2[-which(all_genes2 %in% shared_genes | all_genes2 %in% shared_genes2)]))
fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))



###pathway enrichment analysis for shared genes
library(goseq)
library(reactome.db)
library(org.Mm.eg.db)
library(clusterProfiler)
reactome_shared_KS1 <- getReactomeEnrichedPathways(shared_genes, all_genes)
reactome_shared_KS2 <- getReactomeEnrichedPathways(shared_genes2, all_genes2)
write_csv(reactome_shared_KS1[1:20, ], "neuron_mdem_biorxiv_1/supp_tables/SuppTable15_KS1_shared_3_cell_types_top_pathways.csv")
write_csv(reactome_shared_KS2[1:20, ], "neuron_mdem_biorxiv_1/supp_tables/SuppTable16_KS2_shared_3_cell_types_top_pathways.csv")

###figures
#KS1
quartz(file = "KS1_T_vs_B_and_neurons.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mfrow = c(1,3))
plot(density(atac_T_KS1_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS1_granges, KS1_neurons_B)))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (genome-wide)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 5))
axis(2, at = c(0, 5))
axis(1, at = c(0, 0.5, 1))
lines(density(atac_T_KS1_granges$pvalue[
  unique(queryHits(findOverlaps(atac_T_KS1_granges, unique_B[which(unique_B$padj < 0.1)])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
lines(density(atac_T_KS1_granges$pvalue[
  unique(queryHits(findOverlaps(atac_T_KS1_granges, unique_B[-which(unique_B$padj < 0.1)])))], from = 0, 
  to = 1, bw = 0.035), col = "black", lwd = 2.5)
dev.off()

quartz(file = "KS2_T_vs_B_and_neurons.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mfrow = c(1,3))
plot(density(atac_T_KS2_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS2_granges, KS2_neurons_B)))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (genome-wide)", bty = 'l', xlab = "p-value (B cells)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 5))
axis(2, at = c(0, 5))
axis(1, at = c(0, 0.5, 1))
lines(density(atac_T_KS2_granges$pvalue[
  unique(queryHits(findOverlaps(atac_T_KS2_granges, unique_B_2[which(unique_B_2$pvalue <= 
                                                                       quantile(atac_B_KS2_granges$pvalue, 0.01))])))], from = 0, 
  to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
lines(density(atac_T_KS2_granges$pvalue[
  unique(queryHits(findOverlaps(atac_T_KS2_granges, unique_B_2[-which(unique_B_2$pvalue <= 
                                                                       quantile(atac_B_KS2_granges$pvalue, 0.01))])))], from = 0, 
  to = 1, bw = 0.035), col = "black", lwd = 2.5)
dev.off()






quartz(file = "KS1_T_vs_B_and_neurons.pdf", height = 2.2, width = 6.5, pointsize = 8, type = "pdf")
par(mfrow = c(1,3))
hist(atac_T_KS1_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS1_granges, KS1_neurons_B)))], 
     freq = FALSE, breaks = 40, ylim = c(0, 18),
     col = alpha("red", 0.5), xlab = "KS1 T cell p-values", main = "regions disrupted in B and neurons", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 9, 18))
axis(1, at = c(0, 0.5, 1))
abline(h = 18)
hist(atac_T_KS1_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS1_granges, unique_B[which(unique_B$qvalue < 0.1)])))], 
     freq = FALSE, breaks = 40, ylim = c(0, 18), 
     col = alpha("cornflowerblue", 0.75), main = "regions disrupted in B only", xlab = "KS1 T cell p-values", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 9, 18))
axis(1, at = c(0, 0.5, 1))
abline(h = 18)
hist(atac_T_KS1_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS1_granges, unique_B[-which(unique_B$qvalue < 0.1)])))], 
     freq = FALSE, breaks = 40, ylim = c(0, 18), 
     col = alpha("forest green", 0.75), main = "regions disrupted in neither B nor neurons", xlab = "KS1 T cell p-values", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 9, 18))
axis(1, at = c(0, 0.5, 1))
abline(h = 18)
dev.off()

quartz(file = "KS2_T_vs_B_and_neurons.pdf", height = 2.2, width = 6.5, pointsize = 8, type = "pdf")
par(mfrow = c(1,3))
hist(atac_T_KS2_granges$pvalue[unique(queryHits(findOverlaps(atac_T_KS2_granges, KS2_neurons_B)))], 
     freq = FALSE, breaks = 40, ylim = c(0, 18),
     col = alpha("red", 0.5), xlab = "KS2 T cell p-values", main = "regions disrupted in B and neurons", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 9, 18))
axis(1, at = c(0, 0.5, 1))
abline(h = 18)
hist(atac_T_KS2_granges$pvalue[
  unique(queryHits(findOverlaps(atac_T_KS2_granges, unique_B_2[which(unique_B_2$pvalue <= 
                                                                       quantile(atac_B_KS2_granges$pvalue, 0.01))])))], 
     freq = FALSE, breaks = 40, ylim = c(0, 18), 
     col = alpha("cornflowerblue", 0.75), main = "regions disrupted in B only", xlab = "KS1 T cell p-values", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 9, 18))
axis(1, at = c(0, 0.5, 1))
abline(h = 18)
hist(atac_T_KS2_granges$pvalue[
  unique(queryHits(findOverlaps(atac_T_KS2_granges, unique_B_2[-which(unique_B_2$pvalue <= 
                                                                        quantile(atac_B_KS2_granges$pvalue, 0.01))])))], 
     freq = FALSE, breaks = 40, ylim = c(0, 18), 
     col = alpha("forest green", 0.75), main = "regions disrupted in neither B nor neurons", xlab = "KS1 T cell p-values", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 9, 18))
axis(1, at = c(0, 0.5, 1))
abline(h = 18)
dev.off()



















