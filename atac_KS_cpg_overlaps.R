library(rtracklayer)
session <- browserSession()
genome(session) <- "mm10"
#cpg <- session[["CpG Islands"]]
query <- ucscTableQuery(session, "CpG Islands", GRangesForUCSCGenome("mm10")) #update way of getting the cgis
cpg <- getTable(query)
cpg <- makeGRangesFromDataFrame(cpg, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)


getProportionOfDEgenesInCgi <- function(prom_granges, quantile_values, what = c("less", "greater", "in between")){
  if (what == "less"){
    granges_to_use <- prom_granges[which(prom_granges$pvalue <= 
                                           quantile(prom_granges$pvalue, quantile_values[1], na.rm = TRUE))]
    overlaps <- findOverlaps(granges_to_use, cpg)
    length(unique(queryHits(overlaps)))/length(granges_to_use)
    
  } else if (what == "in between"){
    granges_to_use <- prom_granges[which(prom_granges$pvalue <= 
                                           quantile(prom_granges$pvalue, quantile_values[2], na.rm = TRUE) & 
                                           prom_granges$pvalue >
                                           quantile(prom_granges$pvalue, quantile_values[1], na.rm = TRUE))]
    overlaps <- findOverlaps(granges_to_use, cpg)
    length(unique(queryHits(overlaps)))/length(granges_to_use)
  } else {
    granges_to_use <- prom_granges[which(prom_granges$pvalue > 
                                           quantile(prom_granges$pvalue, quantile_values[1], na.rm = TRUE))]
    overlaps <- findOverlaps(granges_to_use, cpg)
    length(unique(queryHits(overlaps)))/length(granges_to_use)
  }
}
getNumberOfRangesInDecile <- function(prom_granges, quantile_values, what = c("less", "greater", "in between")){
  if (what == "less"){
    granges_to_use <- prom_granges[which(prom_granges$pvalue <= 
                                           quantile(prom_granges$pvalue, quantile_values[1], na.rm = TRUE))]
    length(granges_to_use)
    
  } else if (what == "in between"){
    granges_to_use <- prom_granges[which(prom_granges$pvalue <= 
                                           quantile(prom_granges$pvalue, quantile_values[2], na.rm = TRUE) & 
                                           prom_granges$pvalue >
                                           quantile(prom_granges$pvalue, quantile_values[1], na.rm = TRUE))]
    length(granges_to_use)
  } else {
    granges_to_use <- prom_granges[which(prom_granges$pvalue > 
                                           quantile(prom_granges$pvalue, quantile_values[1], na.rm = TRUE))]
    length(granges_to_use)
  }
}

makeCGIvsAccessibilityPlot <- function(ranges, main_lab, color_to_use, y_limit){
  cpg_overlap_vec <- c(getProportionOfDEgenesInCgi(ranges, 0.1, "less"), 
                                   getProportionOfDEgenesInCgi(ranges, c(0.1, 0.2), "in between"),
                                   getProportionOfDEgenesInCgi(ranges, c(0.2, 0.3), "in between"),
                                   getProportionOfDEgenesInCgi(ranges, c(0.3, 0.4), "in between"),
                                   getProportionOfDEgenesInCgi(ranges, c(0.4, 0.5), "in between"),
                                   getProportionOfDEgenesInCgi(ranges, c(0.5, 0.6), "in between"),
                                   getProportionOfDEgenesInCgi(ranges, c(0.6, 0.7), "in between"),
                                   getProportionOfDEgenesInCgi(ranges, c(0.7, 0.8), "in between"),
                                   getProportionOfDEgenesInCgi(ranges, c(0.8, 0.9), "in between"),
                                   getProportionOfDEgenesInCgi(ranges, 0.9, "greater"))
  
  number_of_proms <- c(getNumberOfRangesInDecile(ranges, 0.1, "less"), 
                       getNumberOfRangesInDecile(ranges, c(0.1, 0.2), "in between"),
                       getNumberOfRangesInDecile(ranges, c(0.2, 0.3), "in between"),
                       getNumberOfRangesInDecile(ranges, c(0.3, 0.4), "in between"),
                       getNumberOfRangesInDecile(ranges, c(0.4, 0.5), "in between"),
                       getNumberOfRangesInDecile(ranges, c(0.5, 0.6), "in between"),
                       getNumberOfRangesInDecile(ranges, c(0.6, 0.7), "in between"),
                       getNumberOfRangesInDecile(ranges, c(0.7, 0.8), "in between"),
                       getNumberOfRangesInDecile(ranges, c(0.8, 0.9), "in between"),
                       getNumberOfRangesInDecile(ranges, 0.9, "greater"))
  
  plot(cpg_overlap_vec, pch = 19, cex = 0, col = color_to_use, bty = 'l', 
       main = main_lab, font.main = 1, ylab = "% of promoters overlapping CGIs", 
       xlab = "differential accessibility p-val decile", 
       xaxt = 'n', ylim = y_limit, yaxt = 'n')
  for(i in 1:10){
    right <- i + 0.2
    left <- i - 0.2
    segments(left, cpg_overlap_vec[i], right, cpg_overlap_vec[i], col = color_to_use, lwd = 1.2)
  }
  for(i in 1:10){
    upper <- cpg_overlap_vec[i] + 
      1.96*sqrt(cpg_overlap_vec[i]*(1-cpg_overlap_vec[i])/number_of_proms[i])
    lower <- cpg_overlap_vec[i] - 
      1.96*sqrt(cpg_overlap_vec[i]*(1-cpg_overlap_vec[i])/number_of_proms[i])
    segments(i, lower, i, upper, col = color_to_use, lwd = 1.2)
  }
  axis(1, at = 1:10, cex.axis = 0.62)
  axis(2, at = c(0.4, 0.8))
}


quartz(file = "atac_cpg_overlaps_KS1_neurons.pdf", height = 2.4, width = 2, pointsize = 8, type = "pdf")
makeCGIvsAccessibilityPlot(res_proms, "KS1 - neurons", "orange", c(0.15, 0.95))
dev.off()
quartz(file = "atac_cpg_overlaps_KS1_B.pdf", height = 2.4, width = 2, pointsize = 8, type = "pdf")
makeCGIvsAccessibilityPlot(atac_B_KS1_granges_proms, "KS1 - B", alpha("black", 0.62), c(0.15, 0.95))
dev.off()
quartz(file = "atac_cpg_overlaps_KS1_T.pdf", height = 2.4, width = 2, pointsize = 8, type = "pdf")
makeCGIvsAccessibilityPlot(atac_T_KS1_granges_proms, "KS1 - T", alpha("black", 0.62), c(0.15, 0.95))
dev.off()
quartz(file = "atac_cpg_overlaps_KS2_neurons.pdf", height = 2.4, width = 2, pointsize = 8, type = "pdf")
makeCGIvsAccessibilityPlot(res2_proms, "KS2 - neurons", "orange", c(0.25, 0.95))
dev.off()
quartz(file = "atac_cpg_overlaps_KS2_B.pdf", height = 2.4, width = 2, pointsize = 8, type = "pdf")
makeCGIvsAccessibilityPlot(atac_B_KS2_granges_proms, "KS2 - B", alpha("black", 0.62), c(0.25, 0.95))
dev.off()
quartz(file = "atac_cpg_overlaps_KS2_T.pdf", height = 2.4, width = 2, pointsize = 8, type = "pdf")
makeCGIvsAccessibilityPlot(atac_T_KS2_granges_proms, "KS2 - T", alpha("black", 0.62), c(0.25, 0.95))
dev.off()


###polycomb regions only


#####logFC for top disrupted CpG
top_cpg_ks1 <- findOverlaps(cpg, res_proms[which(res_proms$pvalue < quantile(res_proms$pvalue, 0.1))])
bottom_cpg_ks1 <- findOverlaps(cpg, res_proms[which(res_proms$pvalue > quantile(res_proms$pvalue, 0.9))])


ks2_top_ks1_cpg_overlaps <- unique(queryHits(findOverlaps(res2_proms, cpg[unique(queryHits(top_cpg_ks1))])))
ks2_bottom_ks1_cpg_overlaps <- unique(queryHits(findOverlaps(res2_proms, cpg[unique(queryHits(bottom_cpg_ks1))])))

ks2_ks1_top_cpg <- getOverlapPi0(res2_proms, cpg[unique(queryHits(top_cpg_ks1))])
ks2_ks1_bottom_cpg <- getOverlapPi0(res2_proms, cpg[unique(queryHits(bottom_cpg_ks1))])

permutation_dist_KS2_topKS1_cpg <- getOverlapNullDistribution1(res2_proms, cpg[unique(queryHits(top_cpg_ks1))])
permutation_dist_KS2_bottomKS1_cpg <- getOverlapNullDistribution1(res2_proms, cpg[unique(queryHits(bottom_cpg_ks1))])


top_ks1 <- res_proms[which(res_proms$pvalue < quantile(res_proms$pvalue, 0.1))]
top_ks1_cpg <- top_ks1[unique(queryHits(findOverlaps(top_ks1, cpg)))]
top_ks1_cpg$log2FoldChange <- -top_ks1_cpg$log2FoldChange #so that positive logFC means increased in mutants compared to wt

top_ks2 <- res2_proms[which(res2_proms$pvalue < quantile(res2_proms$pvalue, 0.1))]
top_ks2_cpg <- top_ks2[unique(queryHits(findOverlaps(top_ks2, cpg)))]
top_ks2_cpg$log2FoldChange <- -top_ks2_cpg$log2FoldChange #so that positive logFC means increased in mutants compared to wt


quartz(file = "atac_neuron_logFC.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
hist(top_ks1_cpg$log2FoldChange, breaks = 45, freq = FALSE, main = "KS1 - neurons", xlab = "log2FoldChange", 
     xaxt = 'n', yaxt = 'n', bty = 'l', col = "orange", font.main = 1, xlim = c(-0.75, 0.75), lty = 0)
axis(1, at = c(-0.5, 0, 0.5))
axis(2, at = c(0, 3.5, 7))
hist(top_ks2_cpg$log2FoldChange, breaks = 40, freq = FALSE, main = "KS2 - neurons", xlab = "log2FoldChange", 
     xaxt = 'n', yaxt = 'n', bty = 'l', col = "orange", font.main = 1, lty = 0)
axis(1, at = c(-0.75, 0, 0.75))
axis(2, at = c(0, 4.5))
dev.off()


quartz(file = "atac_neurons_prc2.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))

genome_wide_ranges <- res_granges
overlaps <- findOverlaps(genome_wide_ranges, proms_mouse)
results_proms <- genome_wide_ranges[queryHits(overlaps)]
results_proms$gene_id <- proms_mouse$gene_id[subjectHits(overlaps)]

res_proms_ezh2 <- results_proms[which(results_proms$gene_id %in% mouse_ezh2_targets_strong)]
cpg_proms_ezh2 <- cpg[unique(queryHits(findOverlaps(cpg, res_proms_ezh2)))]
cpg_proms_non_ezh2 <- cpg[-unique(queryHits(findOverlaps(cpg, res_proms_ezh2)))]
plot(density(res_proms$pvalue[unique(queryHits(findOverlaps(res_proms, cpg_proms_non_ezh2)))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (neurons)", bty = 'l', xlab = "p-value", 
     font.main = 1, yaxt = 'n', xaxt = 'n')
lines(density(res_proms$pvalue[unique(queryHits(findOverlaps(res_proms, cpg_proms_ezh2)))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
lines(density(res_proms$pvalue[-unique(queryHits(findOverlaps(res_proms, cpg)))], from = 0, 
              to = 1, bw = 0.035), col = "black", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 5.5))
legend <- legend("topright", legend = c("CpG island promoters\n(non-PRC2 targets)", 
                                        "CpG island promoters\n(PRC2 targets)", 
                                        "non-CpG island promoters"), col = c(alpha("red", 0.57), "cornflowerblue", "black"),
                 lwd = 2.5, bty = 'n', cex = 0.57)


genome_wide_ranges <- res2_granges
overlaps <- findOverlaps(genome_wide_ranges, proms_mouse)
results_proms <- genome_wide_ranges[queryHits(overlaps)]
results_proms$gene_id <- proms_mouse$gene_id[subjectHits(overlaps)]

res2_proms_ezh2 <- results_proms[which(results_proms$gene_id %in% mouse_ezh2_targets_strong)]
cpg_proms2_ezh2 <- cpg[unique(queryHits(findOverlaps(cpg, res2_proms_ezh2)))]
cpg_proms2_non_ezh2 <- cpg[-unique(queryHits(findOverlaps(cpg, res2_proms_ezh2)))]
plot(density(res2_proms$pvalue[unique(queryHits(findOverlaps(res2_proms, cpg_proms2_non_ezh2)))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS2 (neurons)", bty = 'l', xlab = "p-value", 
     font.main = 1, yaxt = 'n', xaxt = 'n')
lines(density(res2_proms$pvalue[unique(queryHits(findOverlaps(res2_proms, cpg_proms2_ezh2)))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
lines(density(res2_proms$pvalue[-unique(queryHits(findOverlaps(res2_proms, cpg)))], from = 0, 
              to = 1, bw = 0.035), col = "black", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 5.5))
dev.off()


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


