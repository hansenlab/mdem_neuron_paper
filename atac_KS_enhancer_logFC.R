regulatory_regions_mouse <- promoters(edb_mouse, filter = TxBiotypeFilter("protein_coding"), 
                                      upstream = 100000, downstream = 100000, columns = c("gene_name", "tx_id", "gene_id"))
#genome(seqinfo(proms_mouse)) <- "mm10"
genome(seqinfo(regulatory_regions_mouse)) <- "mm10"
#seqlevelsStyle(proms_mouse) <- "ucsc"
seqlevelsStyle(regulatory_regions_mouse) <- "ucsc"
#proms_mouse <- proms_mouse[-which(duplicated(proms_mouse$gene_id))]
regulatory_regions_mouse <- regulatory_regions_mouse[which(seqnames(regulatory_regions_mouse) %in% 
                                                             seqnames(Mmusculus)[1:21])]


getGeneLevelResultsList <- function(genome_wide_ranges){
  genome_wide_ranges <- genome_wide_ranges[-unique(queryHits(findOverlaps(genome_wide_ranges, proms_mouse)))]
  overlaps <- findOverlaps(genome_wide_ranges, regulatory_regions_mouse)
  results_proms <- genome_wide_ranges[queryHits(overlaps)]
  results_proms$gene_id <- regulatory_regions_mouse$gene_id[subjectHits(overlaps)]
  results_proms$gene_name <- regulatory_regions_mouse$gene_name[subjectHits(overlaps)]
  results_proms$corresponding_prom_strand <- as.character(strand(regulatory_regions_mouse)[subjectHits(overlaps)])
  results_proms$corresponding_prom_start <- start(regulatory_regions_mouse)[subjectHits(overlaps)]
  results_proms$corresponding_prom_end <- end(regulatory_regions_mouse)[subjectHits(overlaps)]
  
  gene_level <- split(results_proms, results_proms$gene_id)
  gene_level <- endoapply(gene_level, function(xx) unique(xx))
  gene_level
}

gene_level_results_list_neurons_KS1_regulatory <- getGeneLevelResultsList(res_granges)
gene_level_results_list_neurons_KS2_regulatory <- getGeneLevelResultsList(res2_granges)
gene_level_results_list_B_KS1_regulatory <- getGeneLevelResultsList(atac_B_KS1_granges)
gene_level_results_list_B_KS2_regulatory <- getGeneLevelResultsList(atac_B_KS2_granges)
gene_level_results_list_T_KS1_regulatory <- getGeneLevelResultsList(atac_T_KS1_granges)
gene_level_results_list_T_KS2_regulatory  <- getGeneLevelResultsList(atac_T_KS2_granges)


all_regulatory <- list(gene_level_results_list_neurons_KS1_regulatory, 
              gene_level_results_list_neurons_KS2_regulatory, 
              gene_level_results_list_B_KS1_regulatory, 
              gene_level_results_list_B_KS2_regulatory,
              gene_level_results_list_T_KS1_regulatory,
              gene_level_results_list_T_KS2_regulatory)

all_prom_results <- list(gene_level_results_neurons_KS1, 
                  gene_level_results_neurons_KS2, 
                  gene_level_results_B_KS1, 
                  gene_level_results_B_KS2,
                  gene_level_results_T_KS1,
                  gene_level_results_T_KS2)

top_positive <- lapply(1:6, function(xx) {
  reg <- all_regulatory[[xx]]
  prom_level <- all_prom_results[[xx]]
  reg[which(names(reg) %in% prom_level$gene_ids[
    which(prom_level$pval <= quantile(prom_level$pval, 0.1) & 
            prom_level$logFC > 0)])]
})

bottom_positive <- lapply(1:6, function(xx) {
  reg <- all_regulatory[[xx]]
  prom_level <- all_prom_results[[xx]]
  reg[which(names(reg) %in% prom_level$gene_ids[
    which(prom_level$pval > quantile(prom_level$pval, 0.9) & 
            prom_level$logFC > 0)])]
})

top_negative <- lapply(1:6, function(xx) {
  reg <- all_regulatory[[xx]]
  prom_level <- all_prom_results[[xx]]
  reg[which(names(reg) %in% prom_level$gene_ids[
    which(prom_level$pval <= quantile(prom_level$pval, 0.1) & 
            prom_level$logFC < 0)])]
})

bottom_negative <- lapply(1:6, function(xx) {
  reg <- all_regulatory[[xx]]
  prom_level <- all_prom_results[[xx]]
  reg[which(names(reg) %in% prom_level$gene_ids[
    which(prom_level$pval > quantile(prom_level$pval, 0.9) & 
            prom_level$logFC < 0)])]
})
####
names(top_positive) <- names(bottom_positive) <- names(top_negative) <- names(bottom_negative) <- 
  c("KS1_neurons", "KS2_neurons", "KS1_B", "KS2_B", "KS1_T", "KS2_T")


####negative
logFC_min_pval_top_KS1_neurons_positive <- sapply(top_positive$KS1_neurons, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS1_neurons_positive <- sapply(bottom_positive$KS1_neurons, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])

logFC_min_pval_top_KS2_neurons_positive <- sapply(top_positive$KS2_neurons, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS2_neurons_positive <- sapply(bottom_positive$KS2_neurons, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])

#
logFC_min_pval_top_KS1_B_positive <- sapply(top_positive$KS1_B, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS1_B_positive <- sapply(bottom_positive$KS1_B, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])


logFC_min_pval_top_KS2_B_positive <- sapply(top_positive$KS2_B, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS2_B_positive <- sapply(bottom_positive$KS2_B, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])

#
logFC_min_pval_top_KS1_T_positive <- sapply(top_positive$KS1_T, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS1_T_positive <- sapply(bottom_positive$KS1_T, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])

#
logFC_min_pval_top_KS2_T_positive <- sapply(top_positive$KS2_T, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS2_T_positive <- sapply(bottom_positive$KS2_T, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
####


####positive
logFC_min_pval_top_KS1_neurons_negative <- sapply(top_negative$KS1_neurons, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS1_neurons_negative <- sapply(bottom_negative$KS1_neurons, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])

logFC_min_pval_top_KS2_neurons_negative <- sapply(top_negative$KS2_neurons, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS2_neurons_negative <- sapply(bottom_negative$KS2_neurons, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])

#
logFC_min_pval_top_KS1_B_negative <- sapply(top_negative$KS1_B, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS1_B_negative <- sapply(bottom_negative$KS1_B, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])


logFC_min_pval_top_KS2_B_negative <- sapply(top_negative$KS2_B, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS2_B_negative <- sapply(bottom_negative$KS2_B, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])

#
logFC_min_pval_top_KS1_T_negative <- sapply(top_negative$KS1_T, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS1_T_negative <- sapply(bottom_negative$KS1_T, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])

#
logFC_min_pval_top_KS2_T_negative <- sapply(top_negative$KS2_T, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
logFC_min_pval_bottom_KS2_T_negative <- sapply(bottom_negative$KS2_T, function(xx) 
  xx$log2FoldChange[which.min(xx$pvalue)])
####




###########negative
quartz(file = "KS1_enhancers_neg_logFC.pdf", height = 2.2, width = 6.6, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
##KS1
plot(density(-logFC_min_pval_bottom_KS1_neurons_positive, from = -1.5, to = 1.5), lwd = 2, bty = 'l', xlab = "log2(FoldChange)", 
     main = "KS1 neurons", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.5), col = "deep pink", font.main = 1)
#lines(density(-logFC_min_pval_top_KS1_neurons_positive, from = -1.5, to = 1.5), lwd = 2, col = "orange")
abline(v = -logFC_min_pval_top_KS1_neurons_positive, lwd = 1, col = "orange")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 1))

plot(density(-logFC_min_pval_top_KS1_B_positive, from = -1.5, to = 1.5), lwd = 2, bty = 'l', xlab = "", 
     main = "KS1 B", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.5), col = "orange", font.main = 1)
lines(density(-logFC_min_pval_bottom_KS1_B_positive, from = -1.5, to = 1.5), lwd = 2, col = "deep pink")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 1.2))

plot(density(-logFC_min_pval_top_KS1_T_positive, from = -1.5, to = 1.5), lwd = 2, bty = 'l', xlab = "", 
     main = "KS1 T", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.5), col = "orange", font.main = 1)
lines(density(-logFC_min_pval_bottom_KS1_T_positive, from = -1.5, to = 1.5), lwd = 2, col = "deep pink")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 1.2))
dev.off()

#KS2
quartz(file = "KS2_enhancers_neg_logFC.pdf", height = 2.2, width = 6.6, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
plot(density(-logFC_min_pval_top_KS2_neurons_positive, from = -2, to = 1.5), lwd = 2, bty = 'l', xlab = "log2(FoldChange)", 
     main = "KS2 neurons", xaxt = 'n', yaxt = 'n', xlim = c(-2, 1.5), col = "orange", font.main = 1)
lines(density(-logFC_min_pval_bottom_KS2_neurons_positive, from = -2, to = 1.5), lwd = 2, col = "deep pink")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 1.1))

plot(density(-logFC_min_pval_top_KS2_B_positive, from = -1.5, to = 1.5), lwd = 2, bty = 'l', xlab = "", 
     main = "KS2 B", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.5), col = "orange", font.main = 1)
lines(density(-logFC_min_pval_bottom_KS2_B_positive, from = -1.5, to = 1.5), lwd = 2, col = "deep pink")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 3))

plot(density(-logFC_min_pval_top_KS2_T_positive, from = -1.5, to = 1.5), lwd = 2, bty = 'l', xlab = "", 
     main = "KS2 T", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.5), col = "orange", font.main = 1)
lines(density(-logFC_min_pval_bottom_KS2_T_positive, from = -1.5, to = 1.5), lwd = 2, col = "deep pink")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 1.2))
#legend <- legend("topright", legend = c("nearby promoters w/ top 10% p-values", 
#                                       "nearby promoters w/ bottom 10% p-values"), col = c("orange", 
#                                                                                           "deep pink"), lwd = 2.5, bty = 'n', cex = 0.75)

dev.off()

###########positive
##KS1
quartz(file = "KS1_enhancers_pos_logFC.pdf", height = 2.2, width = 6.6, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
plot(density(-logFC_min_pval_top_KS1_neurons_negative, from = -1.5, to = 1.5), lwd = 2, bty = 'l', xlab = "log2(FoldChange)", 
     main = "KS1 neurons", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.5), col = "orange", font.main = 1)
lines(density(-logFC_min_pval_bottom_KS1_neurons_negative, from = -1.5, to = 1.5), lwd = 2, col = "deep pink")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 1.2))

plot(density(-logFC_min_pval_bottom_KS1_B_negative, from = -1.5, to = 1.75), lwd = 2, bty = 'l', xlab = "", 
     main = "KS1 B", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.75), col = "deep pink", font.main = 1)
lines(density(-logFC_min_pval_top_KS1_B_negative, from = -1.5, to = 1.75), lwd = 2, col = "orange")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 1.2))

plot(density(-logFC_min_pval_top_KS1_T_negative, from = -1.5, to = 1.5), lwd = 2, bty = 'l', xlab = "", 
     main = "KS1 T", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.5), col = "orange", font.main = 1)
lines(density(-logFC_min_pval_bottom_KS1_T_negative, from = -1.5, to = 1.5), lwd = 2, col = "deep pink")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 1))
dev.off()

#KS2
quartz(file = "KS2_enhancers_pos_logFC.pdf", height = 2.2, width = 6.6, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
plot(density(-logFC_min_pval_top_KS2_neurons_negative, from = -1.5, to = 1.5), lwd = 2, bty = 'l', xlab = "log2(FoldChange)", 
     main = "KS2 neurons", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.5), col = "orange", font.main = 1)
lines(density(-logFC_min_pval_bottom_KS2_neurons_negative, from = -1.5, to = 1.5), lwd = 2, col = "deep pink")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 2))

plot(density(-logFC_min_pval_top_KS2_B_negative, from = -1.5, to = 1.5), lwd = 2, bty = 'l', xlab = "", 
     main = "KS2 B", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.5), col = "orange", font.main = 1)
lines(density(-logFC_min_pval_bottom_KS2_B_negative, from = -1.5, to = 1.5), lwd = 2, col = "deep pink")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 1.35))

plot(density(-logFC_min_pval_top_KS2_T_negative, from = -1.5, to = 1.5), lwd = 2, bty = 'l', xlab = "", 
     main = "KS2 T", xaxt = 'n', yaxt = 'n', xlim = c(-1.5, 1.5), col = "orange", font.main = 1)
lines(density(-logFC_min_pval_bottom_KS2_T_negative, from = -1.5, to = 1.5), lwd = 2, col = "deep pink")
axis(1, at = c(-1.5, 0, 1.5))
axis(2, at = c(0, 1))
dev.off()


