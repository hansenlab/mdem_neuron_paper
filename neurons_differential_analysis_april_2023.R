load(file = "count_reads_in_combined_peaks_from_condition_specific_bam_files_neurons_kabuki_cohorts.rda")

peaks_by_samples <- getPeaksBySamplesMatrix("count_reads_in_combined_peaks_from_condition_specific_bam_files_neurons_kabuki_cohorts.rda", 
                                            "neurons")
sample_info <- data.frame(genotype = rep("WT", 28), stringsAsFactors = FALSE)
sample_info$genotype[grep('KS1', colnames(peaks_by_samples))] <- "KS1"
sample_info$genotype[grep('KS2', colnames(peaks_by_samples))] <- "KS2"
sample_info$cohort[grep('Cohort8', colnames(peaks_by_samples))] <- "8"
sample_info$cohort[grep('10-', colnames(peaks_by_samples))] <- "1"
sample_info$cohort[grep('Cohort9', colnames(peaks_by_samples))] <- "9"

getDDSObject <- function(peaks_by_samples_matrix, sample_info_df, mutant_genotype, cohort_names){
  features_by_samples_mat <- peaks_by_samples_matrix[, which(sample_info_df$genotype %in% c("WT", mutant_genotype) & 
                                                        sample_info_df$cohort %in% cohort_names)]
  sample_info2 <- sample_info_df[which(sample_info_df$genotype %in% c("WT", mutant_genotype) & 
                                         sample_info_df$cohort %in% cohort_names), , drop = FALSE]
  
  sampleTable <- data.frame(condition = factor(sample_info2$genotype))
  dds <- DESeqDataSetFromMatrix(round(features_by_samples_mat), sampleTable, design = ~condition)
  
  idx <- rowMedians(counts(dds)) > 10
  dds <- dds[idx, ]
  dds
}

###KS1 cohort 9 excluding problematic label
mat <- peaks_by_samples[, -which(colnames(peaks_by_samples) == "Cohort9_2-22_wt")]
sampinf <- sample_info[-which(colnames(peaks_by_samples) == "Cohort9_2-22_wt"), ]
dds <- getDDSObject(mat, sampinf, "KS1", "9")
dat <- counts(dds)
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(dat, mod, mod0)

ddssva <- dds
ddssva$SV1 <- svseq$sv[, 1]
#ddssva$SV2 <- svseq$sv[, 2]
design(ddssva) <- formula(~ SV1 + condition)
ddssva <- DESeq(ddssva)

dds <- ddssva
res <- results(dds)


###KS2 cohort 1 and 8
dds <- getDDSObject(peaks_by_samples, sample_info, "KS2", c("1", "8"))
dat <- counts(dds)
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(dat, mod, mod0)

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
design(ddssva) <- formula(~ SV1 + condition)
ddssva <- DESeq(ddssva)

dds <- ddssva
res2 <- results(dds)





###############
peak_df <- count_reads_in_peaks$annotation
rownames(peak_df) <- peak_df$GeneID

res$chr <- peak_df[rownames(res), "Chr"]
res$start <- peak_df[rownames(res), "Start"]
res$end <- peak_df[rownames(res), "End"]
res_granges <- makeGRangesFromDataFrame(res, keep.extra.columns = TRUE)
prom_overlaps <- findOverlaps(res_granges, proms_mouse)
res_proms <- res_granges[queryHits(prom_overlaps)]
res_proms$gene_id <- proms_mouse$gene_id[subjectHits(prom_overlaps)]
res_proms$gene_name <- proms_mouse$gene_name[subjectHits(prom_overlaps)]
res_proms <- res_proms[-which(duplicated(res_proms))]

res2$chr <- peak_df[rownames(res2), "Chr"]
res2$start <- peak_df[rownames(res2), "Start"]
res2$end <- peak_df[rownames(res2), "End"]
res2_granges <- makeGRangesFromDataFrame(res2, keep.extra.columns = TRUE)
prom_overlaps2 <- findOverlaps(res2_granges, proms_mouse)
res2_proms <- res2_granges[queryHits(prom_overlaps2)]
res2_proms$gene_id <- proms_mouse$gene_id[subjectHits(prom_overlaps2)]
res2_proms$gene_name <- proms_mouse$gene_name[subjectHits(prom_overlaps2)]
res2_proms <- res2_proms[-which(duplicated(res2_proms))]

res_non_proms <- res_granges[-queryHits(prom_overlaps)]
#res_non_proms <- res_non_proms[-which(duplicated(res_non_proms))]

res2_non_proms <- res2_granges[-queryHits(prom_overlaps2)]
#res2_non_proms <- res2_non_proms[-which(duplicated(res2_non_proms))]


###supp tables
res_granges$qvalue <- qvalue(res_granges$pvalue, fdr.level = 0.1)$qvalues
res2_granges$qvalue <- qvalue(res2_granges$pvalue, fdr.level = 0.1)$qvalues

res_df <- annoGR2DF(res_granges)
res2_df <- annoGR2DF(res2_granges)
res_df <- res_df[which(res_df$qvalue <= 0.1), ]
res2_df <- res_df[which(res2_df$qvalue <= 0.1), ]
res_df$log2FoldChange <- -res_df$log2FoldChange
res2_df$log2FoldChange <- -res2_df$log2FoldChange
write_csv(res_df[order(res_df$qvalue), 
                 c("chr", "start", "end", "width", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "qvalue")], 
          "KS1_neurons_atac_most_significant.csv")
write_csv(res2_df[order(res2_df$qvalue), 
                  c("chr", "start", "end", "width", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "qvalue")],
          "KS2_neurons_atac_most_significant.csv")



##############
####some QC figures
####p-val distributions stratified based on overlap of kmt2d chip seq peaks vs not
quartz(file = "KS1_neurons_vs_kmt2d_chip_ranges.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(res_proms$pvalue[unique(queryHits(findOverlaps(res_proms, kmt2d_chip_granges)))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (promoters)", bty = 'l', xlab = "p-value (neurons)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 8))
lines(density(res_proms$pvalue[-unique(queryHits(findOverlaps(res_proms, kmt2d_chip_granges)))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 8))
legend <- legend("topright", legend = c("ATAC peaks overlapping\nKMT2D ChIP-seq peaks", 
                                        "other ATAC peaks"), col = c(alpha("red", 0.57), 
                                                                     "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.75)


plot(density(res_non_proms$pvalue[unique(queryHits(findOverlaps(res_non_proms, kmt2d_chip_granges)))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "KS1 (outside promoters)", bty = 'l', xlab = "p-value (neurons)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 2.5))
lines(density(res_non_proms$pvalue[-unique(queryHits(findOverlaps(res_non_proms, kmt2d_chip_granges)))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 2.5))
dev.off()


quartz(file = "KS1_vs_KS2_neurons.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
plot(density(res2_proms$pvalue[unique(queryHits(findOverlaps(res2_proms, 
                                                               res_granges[which(res_granges$qvalue < 0.1)])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "neurons", bty = 'l', xlab = "p-value (KS2)", 
     font.main = 1, yaxt = 'n', xaxt = 'n')
lines(density(res2_proms$pvalue[unique(queryHits(findOverlaps(res2_proms, 
                                                                res_granges[-which(res_granges$qvalue < 0.1)])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 7.5))
legend <- legend("topright", legend = c("top disrupted KS1 peaks", 
                                        "other KS1 peaks"), col = c(alpha("red", 0.57), 
                                                                    "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.75)

plot(density(res2_non_proms$pvalue[unique(queryHits(findOverlaps(res2_non_proms, 
                                                             res_granges[which(res_granges$qvalue < 0.1)])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "neurons", bty = 'l', xlab = "p-value (KS2)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 3.62))
lines(density(res2_non_proms$pvalue[unique(queryHits(findOverlaps(res2_non_proms, 
                                                              res_granges[-which(res_granges$qvalue < 0.1)])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 3.5))
legend <- legend("topright", legend = c("top disrupted KS1 peaks", 
                                        "other KS1 peaks"), col = c(alpha("red", 0.57), 
                                                                    "cornflowerblue"), lwd = 2.5, bty = 'n', cex = 0.75)

dev.off()


###KS1 overlap with chip calculation
diff <- which(res_proms$pvalue <= quantile(res_proms$pvalue, 0.1))
non_diff <- which(res_proms$pvalue >= quantile(res_proms$pvalue, 0.9))
kmt2d_bound <- unique(queryHits(findOverlaps(res_proms, kmt2d_chip_granges)))
all <- 1:length(res_proms)
kmt2d_non_bound <- all[-unique(queryHits(findOverlaps(res_proms, kmt2d_chip_granges)))]
n11 <- length(intersect(diff, kmt2d_bound))
n12 <- length(intersect(diff, kmt2d_non_bound))
n21 <- length(intersect(non_diff, kmt2d_bound))
n22 <- length(intersect(non_diff, kmt2d_non_bound))
fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))



