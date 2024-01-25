mat <- peaks_by_samples[, -which(colnames(peaks_by_samples) == "Cohort8-12_KS1")]
sampinf <- sample_info[-which(colnames(peaks_by_samples) == "Cohort8-12_KS1"), ]
dds <- getDDSObject(mat, sampinf, "KS1", "8") #getDDSObject() function defined in the "neurons_differential_analysis_april_2023.R" script
dat <- counts(dds)
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(dat, mod, mod0)

ddssva <- dds
ddssva$SV1 <- svseq$sv[, 1]
design(ddssva) <- formula(~ SV1 + condition)
ddssva <- DESeq(ddssva)

dds <- ddssva
res_8 <- results(dds)

res_8$chr <- peak_df[rownames(res_8), "Chr"] #peak_df as in the "neurons_differential_analysis_april_2023.R" script
res_8$start <- peak_df[rownames(res_8), "Start"]
res_8$end <- peak_df[rownames(res_8), "End"]
res_8_granges <- makeGRangesFromDataFrame(res_8, keep.extra.columns = TRUE)

getEffectSizeConcordance(res_8_granges, 
                         res_proms[which(res_proms$pvalue <= quantile(res_proms$pvalue, 0.1))])
getEffectSizeConcordance(res_8_granges, 
                         res_non_proms[which(res_non_proms$pvalue <= quantile(res_non_proms$pvalue, 0.05))])

1 - pi0est(res_8_granges$pvalue[unique(queryHits(findOverlaps(
  res_8_granges, res_proms[which(res_proms$pvalue <= quantile(res_proms$pvalue, 0.1))])))], 
  pi0.method = "bootstrap")$pi0

1 - pi0est(res_8_granges$pvalue[unique(queryHits(findOverlaps(
  res_8_granges, res_proms[which(res_proms$pvalue > quantile(res_proms$pvalue, 0.9))])))], 
  pi0.method = "bootstrap")$pi0

1 - pi0est(res_8_granges$pvalue[unique(queryHits(findOverlaps(
  res_8_granges, res_non_proms[which(res_non_proms$pvalue <= quantile(res_non_proms$pvalue, 0.01))])))], 
  pi0.method = "bootstrap")$pi0

1 - pi0est(res_8_granges$pvalue[unique(queryHits(findOverlaps(
  res_8_granges, res_non_proms[which(res_non_proms$pvalue > quantile(res_non_proms$pvalue, 0.99))])))], 
  pi0.method = "bootstrap")$pi0


#####################

quartz(file = "KS1_ATAC_replication_exp.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
plot(density(res_8_granges$pvalue[unique(queryHits(findOverlaps(
  res_8_granges, res_proms[which(res_proms$pvalue <= quantile(res_proms$pvalue, 0.1))])))], from = 0, to = 1, 
             bw = 0.035), 
     col = alpha("red", 0.57), lwd = 2.5, main = "promoter ATAC peaks\nred pi_0 = 73%", bty = 'l', xlab = "p-value (replication experiment)", 
     font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 4.5))
lines(density(res_8_granges$pvalue[unique(queryHits(findOverlaps(
  res_8_granges, res_proms[which(res_proms$pvalue > quantile(res_proms$pvalue, 0.9))])))], from = 0, 
              to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 4.5))
legend <- legend("topright", legend = c("original experiment: 1st decile pvalues", 
                                        "original experiment: 10th decile pvalues"), col = c(alpha("red", 0.57), 
                                                                  "cornflowerblue"), lwd = 2.5, 
                 bty = 'n', cex = 0.50)


plot(density(res_8_granges$pvalue[unique(queryHits(findOverlaps(
  res_8_granges, res_non_proms[which(res_non_proms$pvalue <= quantile(res_non_proms$pvalue, 0.01))])))], from = 0, to = 1, 
  bw = 0.035), 
  col = alpha("red", 0.57), lwd = 2.5, main = "ATAC peaks outside promoters\nred pi_0 = 26%", bty = 'l', xlab = "p-value (replication experiment)", 
  font.main = 1, yaxt = 'n', xaxt = 'n', ylim = c(0, 2.5))
lines(density(res_8_granges$pvalue[unique(queryHits(findOverlaps(
  res_8_granges, res_non_proms[which(res_non_proms$pvalue > quantile(res_non_proms$pvalue, 0.99))])))], from = 0, 
  to = 1, bw = 0.035), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 2.5))
legend <- legend("topright", legend = c("original experiment: 1st %ile pvalues", 
                                        "original experiment: 100th %ile pvalues"), col = c(alpha("red", 0.57), 
                                                                                            "cornflowerblue"), lwd = 2.5, 
                 bty = 'n', cex = 0.50)
dev.off()

