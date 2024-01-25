
getGeneLevelResults <- function(genome_wide_ranges){
  overlaps <- findOverlaps(genome_wide_ranges, proms_mouse)
  results_proms <- genome_wide_ranges[queryHits(overlaps)]
  results_proms$gene_id <- proms_mouse$gene_id[subjectHits(overlaps)]
  results_proms$gene_name <- proms_mouse$gene_name[subjectHits(overlaps)]
  
  gene_level <- split(results_proms, results_proms$gene_id)
  gene_level_df <- data.frame(gene_ids = names(gene_level), 
                              gene_names = sapply(gene_level, function(xx) unique(xx$gene_name)),
                              pval = sapply(gene_level, function(xx) min(xx$pvalue)), 
                              padj = sapply(gene_level, function(xx) min(xx$padj)),
                              logFC = sapply(gene_level, function(xx) xx$log2FoldChange[which.min(xx$pvalue)]))
  gene_level_df
}

gene_level_results_neurons_KS1 <- getGeneLevelResults(res_granges)
gene_level_results_neurons_KS2 <- getGeneLevelResults(res2_granges)
gene_level_results_B_KS1 <- getGeneLevelResults(atac_B_KS1_granges)
gene_level_results_B_KS2 <- getGeneLevelResults(atac_B_KS2_granges)
gene_level_results_T_KS1 <- getGeneLevelResults(atac_T_KS1_granges)
gene_level_results_T_KS2  <- getGeneLevelResults(atac_T_KS2_granges)


ks1_pathways_neurons <- getReactomeEnrichedPathways(gene_level_results_neurons_KS1$gene_ids[which(gene_level_results_neurons_KS1$pval < quantile(gene_level_results_neurons_KS1$pval, 0.05))], 
                                                    gene_level_results_neurons_KS1$gene_ids)
ks2_pathways_neurons <- getReactomeEnrichedPathways(gene_level_results_neurons_KS2$gene_ids[which(gene_level_results_neurons_KS2$pval < quantile(gene_level_results_neurons_KS2$pval, 0.05))], 
                                                    gene_level_results_neurons_KS2$gene_ids)
ks1_pathways_B <- getReactomeEnrichedPathways(gene_level_results_B_KS1$gene_ids[which(gene_level_results_B_KS1$pval < quantile(gene_level_results_B_KS1$pval, 0.01))], 
                                              gene_level_results_B_KS1$gene_ids)
ks2_pathways_B <- getReactomeEnrichedPathways(gene_level_results_B_KS2$gene_ids[which(gene_level_results_B_KS2$pval < quantile(gene_level_results_B_KS2$pval, 0.01))], 
                                              gene_level_results_B_KS2$gene_ids)

ks1_pathways_T <- getReactomeEnrichedPathways(gene_level_results_T_KS1$gene_ids[which(gene_level_results_T_KS1$pval < quantile(gene_level_results_T_KS1$pval, 0.01))], 
                                              gene_level_results_T_KS1$gene_ids)


ks2_pathways_T <- getReactomeEnrichedPathways(gene_level_results_T_KS2$gene_ids[which(gene_level_results_T_KS2$pval < quantile(gene_level_results_T_KS2$pval, 0.01))], 
                                              gene_level_results_T_KS2$gene_ids)

write_csv(ks1_pathways_neurons[1:20, ], "KS1_atac_neurons_top_pathways.csv")
write_csv(ks2_pathways_neurons[1:20, ], "KS2_atac_neurons_top_pathways.csv")


################
################KEGG longevity pathway
library(KEGGREST)
kegg_pathway_id <- "mmu04213"
pathway_info <- keggGet(kegg_pathway_id)
gene_list <- pathway_info[[1]]$GENE
gene_list <- gene_list[seq(2, length(gene_list), by = 2)]
longevity_gene_names <- gsub(";.*", "", gene_list)

###neurons
indices <- which(gene_level_results_neurons_KS1$gene_names %in% longevity_gene_names)
wilcox_stat_neurons_KS1 <- wilcox.test(gene_level_results_neurons_KS1$pval[indices], 
                                       gene_level_results_neurons_KS1$pval[-indices])$statistic
null_distribution_neurons_KS1 <- replicate(1000, {
  random_indices <- sample(1:length(gene_level_results_neurons_KS1$pval), length(indices))
  random_wilcox_stat <- wilcox.test(gene_level_results_neurons_KS1$pval[random_indices], 
                                    gene_level_results_neurons_KS1$pval[-random_indices])$statistic
  random_wilcox_stat
})
indices <- which(gene_level_results_neurons_KS2$gene_names %in% longevity_gene_names)
wilcox_stat_neurons_KS2 <- wilcox.test(gene_level_results_neurons_KS2$pval[indices], 
                                       gene_level_results_neurons_KS2$pval[-indices])$statistic
null_distribution_neurons_KS2 <- replicate(1000, {
  random_indices <- sample(1:length(gene_level_results_neurons_KS2$pval), length(indices))
  random_wilcox_stat <- wilcox.test(gene_level_results_neurons_KS2$pval[random_indices], 
                                    gene_level_results_neurons_KS2$pval[-random_indices])$statistic
  random_wilcox_stat
})


###B cells
indices <- which(gene_level_results_B_KS1$gene_names %in% longevity_gene_names)
wilcox_stat_B_KS1 <- wilcox.test(gene_level_results_B_KS1$pval[indices], 
                                 gene_level_results_B_KS1$pval[-indices])$statistic
null_distribution_B_KS1 <- replicate(1000, {
  random_indices <- sample(1:length(gene_level_results_B_KS1$pval), length(indices))
  random_wilcox_stat <- wilcox.test(gene_level_results_B_KS1$pval[random_indices], 
                                    gene_level_results_B_KS1$pval[-random_indices])$statistic
  random_wilcox_stat
})
indices <- which(gene_level_results_B_KS2$gene_names %in% longevity_gene_names)
wilcox_stat_B_KS2 <- wilcox.test(gene_level_results_B_KS2$pval[indices], 
                                 gene_level_results_B_KS2$pval[-indices])$statistic
null_distribution_B_KS2 <- replicate(1000, {
  random_indices <- sample(1:length(gene_level_results_B_KS2$pval), length(indices))
  random_wilcox_stat <- wilcox.test(gene_level_results_B_KS2$pval[random_indices], 
                                    gene_level_results_B_KS2$pval[-random_indices])$statistic
  random_wilcox_stat
})

###T cells
indices <- which(gene_level_results_T_KS1$gene_names %in% longevity_gene_names)
wilcox_stat_T_KS1 <- wilcox.test(gene_level_results_T_KS1$pval[indices], 
                                 gene_level_results_T_KS1$pval[-indices])$statistic
null_distribution_T_KS1 <- replicate(1000, {
  random_indices <- sample(1:length(gene_level_results_T_KS1$pval), length(indices))
  random_wilcox_stat <- wilcox.test(gene_level_results_T_KS1$pval[random_indices], 
                                    gene_level_results_T_KS1$pval[-random_indices])$statistic
  random_wilcox_stat
})
indices <- which(gene_level_results_T_KS2$gene_names %in% longevity_gene_names)
wilcox_stat_T_KS2 <- wilcox.test(gene_level_results_T_KS2$pval[indices], 
                                 gene_level_results_T_KS2$pval[-indices])$statistic
null_distribution_T_KS2 <- replicate(1000, {
  random_indices <- sample(1:length(gene_level_results_T_KS2$pval), length(indices))
  random_wilcox_stat <- wilcox.test(gene_level_results_T_KS2$pval[random_indices], 
                                    gene_level_results_T_KS2$pval[-random_indices])$statistic
  random_wilcox_stat
})

quartz(file = "aging_pathway_ranks_KS1.pdf", width = 6, height = 2, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
makeWilcoxPlot(null_distribution_neurons_KS1, wilcox_stat_neurons_KS1, 
               main_lab = "KS1 neurons")
legend <- legend("topright", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
                 cex = 1, lty = "solid", lwd = 2.5)

makeWilcoxPlot(null_distribution_B_KS1, wilcox_stat_B_KS1, 
               main_lab = "KS1 B cells")
makeWilcoxPlot(null_distribution_B_KS1, wilcox_stat_T_KS1, 
               main_lab = "KS1 T cells")
dev.off()

quartz(file = "aging_pathway_ranks_KS2.pdf", width = 6, height = 2, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
makeWilcoxPlot(null_distribution_neurons_KS2, wilcox_stat_neurons_KS2, 
               main_lab = "KS2 neurons")
legend <- legend("topright", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
                 cex = 1, lty = "solid", lwd = 2.5)

makeWilcoxPlot(null_distribution_B_KS2, wilcox_stat_B_KS2, 
               main_lab = "KS2 B cells")
makeWilcoxPlot(null_distribution_B_KS2, wilcox_stat_T_KS2, 
               main_lab = "KS2 T cells")
dev.off()