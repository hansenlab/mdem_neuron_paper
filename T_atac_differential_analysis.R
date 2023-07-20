load(file = "count_reads_in_combined_peaks_from_condition_specific_bam_files_T_cells_only_kabuki_cohort.rda")
peaks_by_samples <- getPeaksBySamplesMatrix("count_reads_in_combined_peaks_from_condition_specific_bam_files_T_cells_only_kabuki_cohort.rda", 
                                            "blood_all_samples")
###note: this script is for T cells, but the exact same analysis is done for B cells if instead of the 
###count_reads_in_combined_peaks_from_condition_specific_bam_files_T_cells_only_kabuki_cohort.rda file
###I load the count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_kabuki_cohort.rda

genotype <- "KS1" #or KS2
celltype <- "T"  #or B 

sample_info <- data.frame(genotype = getGenotypeVector(colnames(peaks_by_samples)), 
                          cell_type = getCellTypeVector(colnames(peaks_by_samples)))

features_by_samples_mat <- peaks_by_samples[, which(sample_info$genotype %in% c("WT", genotype) & 
                                                      sample_info$cell_type == celltype)] #& sample_info$cohort %in% c("KMT2D_1", "KMT2D_2", "KMT2D_3")

sample_info2 <- sample_info[which(sample_info$genotype %in% c("WT", genotype) & 
                                    sample_info$cell_type == celltype), , drop = FALSE]
#sample_info2$genotype <- droplevels(sample_info2$genotype) #this is not needed for RT

sampleTable <- data.frame(condition = factor(sample_info2$genotype))
dds <- DESeqDataSetFromMatrix(round(features_by_samples_mat), sampleTable, design = ~condition)

idx <- rowMedians(counts(dds)) > 10
dat <- counts(dds)[idx,]
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(dat, mod, mod0)

dds <- dds[idx, ]
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- formula(~ SV1 + SV2 + condition)
ddssva <- DESeq(ddssva)

dds <- ddssva
res <- results(dds)
#if (length(which(is.na(res$padj))) > 0){res <- res[-which(is.na(res$padj)), ]} #optional step that was done in the eLife paper (Luperchio et al., 2021). Mainly excludes non-promoter peaks with low counts and doesn't really affect the results

peak_df <- count_reads_in_peaks$annotation
rownames(peak_df) <- peak_df$GeneID

res$chr <- peak_df[rownames(res), "Chr"]
res$start <- peak_df[rownames(res), "Start"]
res$end <- peak_df[rownames(res), "End"]
atac_T_KS1_granges <- makeGRangesFromDataFrame(res, keep.extra.columns = TRUE)



