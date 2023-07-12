####QC for cell types
###the following is based on the peaks_by_samples matrix from the B cell peak locations
###A very similar result is obtained with the T cell peak locations (note that both B and T cell samples 
###have been mapped twice; once onto the peak locations called from the T cell meta-samples, and once onto the 
###peak locations called from the B cell meta-samples)
sample_info <- data.frame(genotype = getGenotypeVector(colnames(peaks_by_samples)), 
                          cell_type = getCellTypeVector(colnames(peaks_by_samples)))
sampleTable <- data.frame(condition = factor(sample_info$cell_type))
dds <- DESeqDataSetFromMatrix(round(peaks_by_samples), sampleTable, design = ~condition)
vstab <- vst(dds)

pcaData <- plotPCA(vstab, returnData = TRUE)

colnames(pcaData)[4] <- "cell_type"
percentVar <- round(100 * attr(pcaData, "percentVar"))


quartz(file = "pca_B_T_cells_atac_ks.pdf", width = 2.4, height = 2, pointsize = 8, type = "pdf")
#par(mfrow = c(1,2))
par(mar = c(4, 4, 1, 1) + 0.1)
#kabuki cohort
plot(pcaData$PC1[which(pcaData$cell_type == "B")], pcaData$PC2[which(pcaData$cell_type == "B")], 
     cex = 1.15, col = "deep pink", xlab = paste0("PC1 (", percentVar[1], "%)"), ylab = paste0("PC2 (", percentVar[2], "%)"), 
     main = "", pch = 19, xlim = c(-50, 50), ylim = c(-3, 15), bty = 'l', xaxt = 'n', yaxt = 'n', 
     font.main = 1)
points(pcaData$PC1[which(pcaData$cell_type == "T")], pcaData$PC2[which(pcaData$cell_type == "T")], 
       cex = 1.15, col = "orange", pch = 19)
legend <- legend("bottom", legend = c("B cells", "T cells"), bty = 'n', cex = 0.82, col = c("deep pink", "orange"), pch = 19)
axis(1, at = c(-40,40))
axis(2, at = c(-1.5, 14))
dev.off()

