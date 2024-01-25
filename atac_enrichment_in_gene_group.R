getEnrichmentInGeneGroup <- function(genes_target, genes_other, gene_group){
  n11 <- length(which(genes_target %in% gene_group))
  n12 <- length(which(!(genes_target %in% gene_group)))
  n21 <- length(which(genes_other %in% gene_group))
  n22 <- length(which(!(genes_other %in% gene_group)))
  fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))
}

###goal is to find gene names whose promoters contain CGIs, in order to then test for expression changes in these
cgi_ids_neurons_KS1 <- res_proms$gene_id[unique(queryHits(findOverlaps(res_proms, cpg)))]
cgi_ids_B_KS1 <- atac_B_KS1_granges_proms$gene_id[unique(queryHits(findOverlaps(atac_B_KS1_granges_proms, cpg)))]
cgi_ids_T_KS1 <- atac_T_KS1_granges_proms$gene_id[unique(queryHits(findOverlaps(atac_T_KS1_granges_proms, cpg)))]

cgi_ids_neurons_KS2 <- res2_proms$gene_id[unique(queryHits(findOverlaps(res2_proms, cpg)))]
cgi_ids_B_KS2 <- atac_B_KS2_granges_proms$gene_id[unique(queryHits(findOverlaps(atac_B_KS2_granges_proms, cpg)))]
cgi_ids_T_KS2 <- atac_T_KS2_granges_proms$gene_id[unique(queryHits(findOverlaps(atac_T_KS2_granges_proms, cpg)))]

################################
quartz(file = "KS1_KS2_ATAC_top_vs_bottom_cpg.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
enrichment_test <- getEnrichmentInGeneGroup(gene_level_results_neurons_KS1$gene_ids[
  which(rank(gene_level_results_neurons_KS1$pval) <= 500)], 
  gene_level_results_neurons_KS1$gene_ids[
    which(rank(gene_level_results_neurons_KS1$pval) >= max(rank(gene_level_results_neurons_KS1$pval)) - 500)], 
  cgi_ids_neurons_KS1)
x.pos <- 1
plot(x.pos, enrichment_test$estimate, 
     pch = 19, cex = 1.35, col = alpha("orange", 1), main = "", bty = 'l', ylab = "enrichment in gene group", 
     xlab = "", ylim = c(0, 16), xlim = c(0.8, 6.2), xaxt = 'n', yaxt = 'n')
#segments(x.pos - 0.2, enrichment_test$estimate, 
#         x.pos + 0.2, enrichment_test$estimate, col = alpha("orange", 1), lwd = 1.4)
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("orange", 1), lwd = 1)

abline(h = 1, lty = "dashed", lwd = 1.25)

enrichment_test <- getEnrichmentInGeneGroup(gene_level_results_B_KS1$gene_ids[
  which(rank(gene_level_results_B_KS1$pval) <= 500)], 
  gene_level_results_B_KS1$gene_ids[
    which(rank(gene_level_results_B_KS1$pval) >= max(rank(gene_level_results_B_KS1$pval)) - 500)], 
  cgi_ids_B_KS1)
x.pos <- 2
points(x.pos, enrichment_test$estimate, pch = 19, cex = 1.35, col = alpha("orange", 1))
#segments(x.pos - 0.2, enrichment_test$estimate, 
#         x.pos + 0.2, enrichment_test$estimate, col = alpha("orange", 1), lwd = 1.4)
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("orange", 1), lwd = 1)

enrichment_test <- getEnrichmentInGeneGroup(gene_level_results_T_KS1$gene_ids[
  which(rank(gene_level_results_T_KS1$pval) <= 500)], 
  gene_level_results_T_KS1$gene_ids[
    which(rank(gene_level_results_T_KS1$pval) >= max(rank(gene_level_results_T_KS1$pval)) - 500)], 
  cgi_ids_T_KS1)
x.pos <- 3
points(x.pos, enrichment_test$estimate, pch = 19, cex = 1.35, col = alpha("orange", 1))
#segments(x.pos - 0.2, enrichment_test$estimate, 
#         x.pos + 0.2, enrichment_test$estimate, col = alpha("orange", 1), lwd = 1.4)
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("orange", 1), lwd = 1)

enrichment_test <- getEnrichmentInGeneGroup(gene_level_results_neurons_KS2$gene_ids[
  which(rank(gene_level_results_neurons_KS2$pval) <= 500)], 
  gene_level_results_neurons_KS2$gene_ids[
    which(rank(gene_level_results_neurons_KS2$pval) >= max(rank(gene_level_results_neurons_KS2$pval)) - 500)], 
  cgi_ids_neurons_KS2)
x.pos <- 4
points(x.pos, enrichment_test$estimate, pch = 19, cex = 1.35, col = alpha("orange", 1))
#segments(x.pos - 0.2, enrichment_test$estimate, 
#         x.pos + 0.2, enrichment_test$estimate, col = alpha("orange", 1), lwd = 1.4)
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("orange", 1), lwd = 1)

enrichment_test <- getEnrichmentInGeneGroup(gene_level_results_B_KS2$gene_ids[
  which(rank(gene_level_results_B_KS2$pval) <= 500)], 
  gene_level_results_B_KS2$gene_ids[
    which(rank(gene_level_results_B_KS2$pval) >= max(rank(gene_level_results_B_KS2$pval)) - 500)], 
  cgi_ids_B_KS2)
x.pos <- 5
points(x.pos, enrichment_test$estimate, pch = 19, cex = 1.35, col = alpha("orange", 1))
#segments(x.pos - 0.2, enrichment_test$estimate, 
#         x.pos + 0.2, enrichment_test$estimate, col = alpha("orange", 1), lwd = 1.4)
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("orange", 1), lwd = 1)

enrichment_test <- getEnrichmentInGeneGroup(gene_level_results_T_KS2$gene_ids[
  which(rank(gene_level_results_T_KS2$pval) <= 500)], 
  gene_level_results_T_KS2$gene_ids[
    which(rank(gene_level_results_T_KS2$pval) >= max(rank(gene_level_results_T_KS2$pval)) - 500)], 
  cgi_ids_T_KS2)
x.pos <- 6
points(x.pos, enrichment_test$estimate, pch = 19, cex = 1.35, col = alpha("orange", 1))
#segments(x.pos - 0.2, enrichment_test$estimate, 
#         x.pos + 0.2, enrichment_test$estimate, col = alpha("orange", 1), lwd = 1.4)
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("orange", 1), lwd = 1)

axis(2, at = c(1, 8, 15))
dev.off()


##
quartz(file = "KS1_KS2_ATAC_top_vs_bottom_em.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
enrichment_test <- getEnrichmentInGeneGroup(toupper(gene_level_results_neurons_KS1$gene_names[
  which(rank(gene_level_results_neurons_KS1$pval) <= 3000)]), 
  toupper(gene_level_results_neurons_KS1$gene_names[
    which(rank(gene_level_results_neurons_KS1$pval) >= max(rank(gene_level_results_neurons_KS1$pval)) - 3000)]),
  rownames(adjmatrix)[1:74])
x.pos <- 1
plot(1, enrichment_test$estimate, 
     pch = 19, cex = 1.35, col = alpha("deep pink", 1), main = "", bty = 'l', ylab = "enrichment in gene group", 
     xlab = "", ylim = c(0, 40), xlim = c(0.8, 6.2), xaxt = 'n', yaxt = 'n')
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("deep pink", 1), lwd = 1)

abline(h = 1, lty = "dashed", lwd = 1.25)

enrichment_test <- getEnrichmentInGeneGroup(toupper(gene_level_results_B_KS1$gene_names[
  which(rank(gene_level_results_B_KS1$pval) <= 3000)]), 
  toupper(gene_level_results_B_KS1$gene_names[
    which(rank(gene_level_results_B_KS1$pval) >= max(rank(gene_level_results_B_KS1$pval)) - 3000)]),
  rownames(adjmatrix)[1:74])
x.pos <- 2
points(x.pos, enrichment_test$estimate, pch = 19, cex = 1.35, col = alpha("deep pink", 1))
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("deep pink", 1), lwd = 1)

enrichment_test <- getEnrichmentInGeneGroup(toupper(gene_level_results_T_KS1$gene_names[
  which(rank(gene_level_results_T_KS1$pval) <= 3000)]), 
  toupper(gene_level_results_T_KS1$gene_names[
    which(rank(gene_level_results_T_KS1$pval) >= max(rank(gene_level_results_T_KS1$pval)) - 3000)]),
  rownames(adjmatrix)[1:74])
x.pos <- 3
points(x.pos, enrichment_test$estimate, pch = 19, cex = 1.35, col = alpha("deep pink", 1))
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("deep pink", 1), lwd = 1)

enrichment_test <- getEnrichmentInGeneGroup(toupper(gene_level_results_neurons_KS2$gene_names[
  which(rank(gene_level_results_neurons_KS2$pval) <= 3000)]), 
  toupper(gene_level_results_neurons_KS2$gene_names[
    which(rank(gene_level_results_neurons_KS2$pval) >= max(rank(gene_level_results_neurons_KS2$pval)) - 3000)]),
  rownames(adjmatrix)[1:74])
x.pos <- 4
points(x.pos, enrichment_test$estimate, pch = 19, cex = 1.35, col = alpha("deep pink", 1))
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("deep pink", 1), lwd = 1)

enrichment_test <- getEnrichmentInGeneGroup(toupper(gene_level_results_B_KS2$gene_names[
  which(rank(gene_level_results_B_KS2$pval) <= 3000)]), 
  toupper(gene_level_results_B_KS2$gene_names[
    which(rank(gene_level_results_B_KS2$pval) >= max(rank(gene_level_results_B_KS2$pval)) - 3000)]),
  rownames(adjmatrix)[1:74])
x.pos <- 5
points(x.pos, enrichment_test$estimate, pch = 19, cex = 1.35, col = alpha("deep pink", 1))
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("deep pink", 1), lwd = 1)

enrichment_test <- getEnrichmentInGeneGroup(toupper(gene_level_results_T_KS2$gene_names[
  which(rank(gene_level_results_T_KS2$pval) <= 3000)]), 
  toupper(gene_level_results_T_KS2$gene_names[
    which(rank(gene_level_results_T_KS2$pval) >= max(rank(gene_level_results_T_KS2$pval)) - 3000)]),
  rownames(adjmatrix)[1:74])
x.pos <- 6
points(x.pos, enrichment_test$estimate, pch = 19, cex = 1.35, col = alpha("deep pink", 1))
segments(x0 = x.pos, y0 = enrichment_test$conf.int[1], x1 = x.pos, y1 = enrichment_test$conf.int[2], 
         col = alpha("deep pink", 1), lwd = 1)
#axis(1, at = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5), labels = rep(c("neurons", "B cells", "T cells"), 2))
axis(2, at = c(1, 15, 30))
#legend("topright", legend = c("CpG island genes", 
#                              "highly co-expressed\nepigenetic machinery genes"), pch = 19, 
#       col = c("orange", alpha("deep pink", 0.7)), bty = 'n')
dev.off()


######from RNA seq
res_B_KS1$gene_name <- sapply(rownames(res_B_KS1), 
                              function(xx) unique(proms_mouse$gene_name[which(proms_mouse$gene_id == xx)]))
res_T_KS1$gene_name <- sapply(rownames(res_T_KS1), 
                              function(xx) unique(proms_mouse$gene_name[which(proms_mouse$gene_id == xx)]))
res_B_KS2$gene_name <- sapply(rownames(res_B_KS2), 
                              function(xx) unique(proms_mouse$gene_name[which(proms_mouse$gene_id == xx)]))
res_T_KS2$gene_name <- sapply(rownames(res_T_KS2), 
                              function(xx) unique(proms_mouse$gene_name[which(proms_mouse$gene_id == xx)]))

###cgi
getEnrichmentInGeneGroup(rownames(res_B_KS1)[which(rank(res_B_KS1$pvalue) <= 500)], 
                         rownames(res_B_KS1)[which(res_B_KS1$pvalue >= max(res_B_KS1$pvalue) - 500)], 
                         cgi_ids_B_KS1)
getEnrichmentInGeneGroup(rownames(res_T_KS1)[which(rank(res_T_KS1$pvalue) <= 500)], 
                         rownames(res_T_KS1)[which(res_T_KS1$pvalue >= max(res_T_KS1$pvalue) - 500)], 
                         cgi_ids_T_KS1)
getEnrichmentInGeneGroup(rownames(res_B_KS2)[which(rank(res_B_KS2$pvalue) <= 500)], 
                         rownames(res_B_KS2)[which(res_B_KS2$pvalue >= max(res_B_KS2$pvalue) - 500)], 
                         cgi_ids_B_KS2)
getEnrichmentInGeneGroup(rownames(res_T_KS2)[which(rank(res_T_KS2$pvalue) <= 500)], 
                         rownames(res_T_KS2)[which(res_T_KS2$pvalue >= max(res_T_KS2$pvalue) - 500)], 
                         cgi_ids_T_KS2)

###highly coexpressed EM
getEnrichmentInGeneGroup(toupper(res_B_KS1$gene_name[
  which(rank(res_B_KS1$pvalue) <= 500)]), 
  toupper(res_B_KS1$gene_name[
    which(rank(res_B_KS1$pvalue) >= max(res_B_KS1$pvalue) - 500)]),
  rownames(adjmatrix)[1:74])

getEnrichmentInGeneGroup(toupper(res_T_KS1$gene_name[
  which(rank(res_T_KS1$pvalue) <= 500)]), 
  toupper(res_T_KS1$gene_name[
    which(rank(res_T_KS1$pvalue) >= max(res_T_KS1$pvalue) - 500)]),
  rownames(adjmatrix)[1:74])

getEnrichmentInGeneGroup(toupper(res_B_KS2$gene_name[
  which(rank(res_B_KS2$pvalue) <= 500)]), 
  toupper(res_B_KS2$gene_name[
    which(rank(res_B_KS2$pvalue) >= max(res_B_KS2$pvalue) - 500)]),
  rownames(adjmatrix)[1:74])

getEnrichmentInGeneGroup(toupper(res_T_KS2$gene_name[
  which(rank(res_T_KS2$pvalue) <= 500)]), 
  toupper(res_T_KS2$gene_name[
    which(rank(res_T_KS2$pvalue) >= max(res_T_KS2$pvalue) - 500)]),
  rownames(adjmatrix)[1:74])




