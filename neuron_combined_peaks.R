setwd('/dcl01/hansen/data/kabuki_mice/HiSeq349/ATAC-seq/brain_samples_merged_bam')
###run the following on the cluster to get genotype specific "meta-bam" files and then call peaks from these files
###this excludes sample  "Cohort9_2-22_wt" which QC suggests is a KS1 and "Cohort8-12_KS1" which QC suggests is a WT
#samtools merge -f -@ 30 all_wt_samples_updated_Fe2023.bam Cohort8-16*noMito_merged.bam Cohort8-19*noMito_merged.bam Cohort8-21*noMito_merged.bam 10-5-2*noMito_merged.bam 10-5-10*noMito_merged.bam 10-5-4*noMito_merged.bam Cohort9-10*noMito_merged.bam Cohort9-6*noMito_merged.bam Cohort9-2_*noMito_merged.bam Cohort9-2-11*noMito_merged.bam Cohort9-2-12*noMito_merged.bam
#samtools merge -f -@ 30 all_KS1_samples_updated_Fe2023.bam Cohort8-22*noMito_merged.bam Cohort8-11*noMito_merged.bam Cohort8-9*noMito_merged.bam Cohort9-8*noMito_merged.bam Cohort9-14*noMito_merged.bam Cohort9-3*noMito_merged.bam Cohort9-1*noMito_merged.bam Cohort9-2-25*noMito_merged.bam Cohort9-2-20*noMito_merged.bam Cohort9-2-4*noMito_merged.bam Cohort9-2-17*noMito_merged.bam

#macs2 callpeak -t all_KS2_samples.bam -n all_KS2_samples -g mm --keep-dup all -f BAMPE --bdg
#macs2 callpeak -t all_KS1_samples_updated_Fe2023.bam -n all_KS1_samples -g mm --keep-dup all -f BAMPE --bdg
#macs2 callpeak -t all_wt_samples_updated_Fe2023.bam -n all_wt_samples -g mm --keep-dup all -f BAMPE --bdg

library(Rsubread)
getGrangesFromNarrowPeaks <- function(narrowPeak_filepath){
  peaks <- read.delim(paste0(narrowPeak_filepath), header = FALSE)
  colnames(peaks) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "fold.enrichment",
                       "log10.pval", "log10.qval", "peak")
  peaks_granges <- GRanges(seqnames = peaks$chrom, IRanges(peaks$chromStart+1, peaks$chromEnd), 
                           fold.enrichment = peaks$fold.enrichment, log10.pval = peaks$log10.pval, 
                           log10.qval = peaks$log10.qval)
  peaks_granges
}

wt_granges <- getGrangesFromNarrowPeaks('all_wt_samples_peaks.narrowPeak')
ks1_granges <- getGrangesFromNarrowPeaks('all_KS1_samples_peaks.narrowPeak')
ks2_granges <- getGrangesFromNarrowPeaks('all_KS2_samples_peaks.narrowPeak')

combined_granges_from_peaks_called_from_all_fragments_neurons <- reduce(unlist(GRangesList(wt_granges, ks1_granges, ks2_granges)))
mm10_regions_to_exclude_granges <- import("/dcl01/hansen/data/kabuki_mice/HiSeq349/ATAC-seq/mm10-blacklist.v2.bed", format = "BED")
combined_granges_from_peaks_called_from_all_fragments_neurons <- combined_granges_from_peaks_called_from_all_fragments_neurons[-unique(
  queryHits(findOverlaps(combined_granges_from_peaks_called_from_all_fragments_neurons, mm10_regions_to_exclude_granges)))]

save(combined_granges_from_peaks_called_from_all_fragments_neurons, 
     file = "combined_granges_from_peaks_called_from_condition_specific_bams_with_all_fragments_neurons.rda")

combined_granges <- combined_granges_from_peaks_called_from_all_fragments_neurons

combined_granges$id <- paste0(seqnames(combined_granges), "_", 
                              as.character(start(combined_granges)), "_", as.character(end(combined_granges)))
#combined_peaks_ann <- createAnnotationFile(combined_granges) ###this function is longer supported so we run the following

combined_peaks_ann <- data.frame(GeneID = combined_granges$id, 
                                 Chr = seqnames(combined_granges), 
                                 Start = start(combined_granges), 
                                 End = end(combined_granges), 
                                 strand = strand(combined_granges), 
                                 stringsAsFactors = FALSE)

bamFiles_filepaths <- list.files(pattern = "noMito_merged.bam")

count_reads_in_peaks <- featureCounts(bamFiles_filepaths, annot.ext = combined_peaks_ann, 
                                      isPairedEnd = TRUE, minOverlap = 3,
                                      requireBothEndsMapped = TRUE, countChimericFragments = FALSE, 
                                      nthreads = 30, countMultiMappingReads = FALSE)

save(count_reads_in_peaks, file = "count_reads_in_combined_peaks_from_condition_specific_bam_files_neurons_kabuki_cohorts.rda")


####public atac
setwd('/dcl01/hansen/data/kabuki_mice/HiSeq349/ATAC-seq/atac_fastq_from_GSE82010')

public_atac_dentate_gyrus <- getGrangesFromNarrowPeaks('all_samples_peaks.narrowPeak')

percentage_shared_wt <- length(unique(queryHits(findOverlaps(public_atac_dentate_gyrus, wt_granges)))) / 
  length(public_atac_dentate_gyrus)
percentage_shared_KS1 <- length(unique(queryHits(findOverlaps(public_atac_dentate_gyrus, ks1_granges)))) / 
  length(public_atac_dentate_gyrus)
percentage_shared_KS2 <- length(unique(queryHits(findOverlaps(public_atac_dentate_gyrus, ks2_granges)))) / 
  length(public_atac_dentate_gyrus)

quartz(file = "neuronal_atac_comparison_Su_et_al.pdf", height = 2.4, width = 2, pointsize = 8, type = "pdf")
plot(1, 0.9643499, pch = 19, col = alpha("red", 0.75), bty = 'l', ylab = "% Su et al. peaks shared",
     main = "", xaxt = 'n', yaxt = 'n', xlim = c(0.8, 3.2), ylim = c(0, 1), xlab = "", cex = 1.25)
points(2, 0.9736306, pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(3, 0.9035804, pch = 19, col = alpha("red", 0.75), cex = 1.25)
abline(v = c(1.5, 2.5), lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(0.1, 0.50, 0.90), labels = c(10, 50, 90))
axis(1, at = c(1, 2, 3), labels = c("WT", "KS1", "KS2"), cex.axis = 0.8, las = 1)
dev.off()

###encode 3 forebrain 
#run in cluster
setwd('/dcl01/hansen/data/kabuki_mice/HiSeq349/ATAC-seq/encode3_forebrain_P0')
encode_forebrain <- import('encode3RenAtacSignalForebrainP0.bw', format = "BigWig") #downloaded from https://hgdownload.soe.ucsc.edu/gbdb/mm10/encode3/atac/
encode_forebrain <- encode_forebrain[-which(encode_forebrain$score == 0)]

length(unique(queryHits(findOverlaps(wt_granges, encode_forebrain)))) / length(wt_granges)
length(unique(queryHits(findOverlaps(ks1_granges, encode_forebrain)))) / length(ks1_granges)
length(unique(queryHits(findOverlaps(ks2_granges, encode_forebrain)))) / length(ks2_granges)


getOverlapsWithEncode <- function(genomic_ranges, name){
  overlaps_proms <- findOverlaps(encode_forebrain, 
                                 genomic_ranges[unique(queryHits(findOverlaps(genomic_ranges, proms_mouse)))])
  overlaps_non_proms <- findOverlaps(encode_forebrain, 
                                 genomic_ranges[-unique(queryHits(findOverlaps(genomic_ranges, proms_mouse)))])
  encode_overlapping_proms <- encode_forebrain[unique(queryHits(overlaps_proms))]
  encode_overlapping_non_proms <- encode_forebrain[unique(queryHits(overlaps_non_proms))]
  encode_non_overlapping <- encode_forebrain[-unique(c(unique(queryHits(overlaps_proms)), 
                                                       unique(queryHits(overlaps_non_proms))))]
  vec1 <- encode_overlapping_proms$score
  vec2 <- encode_overlapping_non_proms$score
  vec3 <- encode_non_overlapping$score
  
  save(vec1, vec2, vec3, file = paste0(name, ".rda"))
} 


getOverlapsWithEncode(wt_granges, "wt_granges_encode_overlap")
getOverlapsWithEncode(ks1_granges, "ks1_granges_encode_overlap")
getOverlapsWithEncode(ks2_granges, "ks2_granges_encode_overlap")

#now run locally to make plots
quartz(file = "wt_encode_overlap.pdf", height = 2.4, width = 2, pointsize = 8, type = "pdf")
load(file = "wt_granges_encode_overlap.rda")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
boxplot(vec1, frame = FALSE, main = "ENCODE vs wild-type", font.main = 1, 
        lty = "solid", col = "orange", yaxt = 'n', ylim = c(0, 19),
        xlim = c(0.8, 3.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "ENCODE ATAC-seq signal", xlab = "")
boxplot(vec2, frame = FALSE, yaxt = 'n',
        lty = "solid", col = "deep pink", medlty = 1, medlwd = 0.8, at = 2, xaxt = 'n', add = TRUE,
        boxlty = 0, staplelwd = 0, outline = FALSE)
boxplot(vec3, frame = FALSE, yaxt = 'n',
        lty = "solid", col = rgb(0,0,0,0.75), medlty = 1, medlwd = 0.8, at = 3, xaxt = 'n', add = TRUE,
        boxlty = 0, staplelwd = 0, outline = FALSE)

legend("topright", legend = c("overlapping\npromoter peaks", "overlapping\nnon-promoter peaks", 
                                "not overlapping\npeaks"), title = "ENCODE forebrain peaks", 
       fill = c("orange", "deep pink", rgb(0,0,0,0.75)), bty = 'n', border = "white",
     cex = 0.75)

axis(2, at = c(0, 7.5, 15))
dev.off()

quartz(file = "ks1_encode_overlap.pdf", height = 2.4, width = 2, pointsize = 8, type = "pdf")
load(file = "ks1_granges_encode_overlap.rda")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
boxplot(vec1, frame = FALSE, main = "ENCODE vs KS1", font.main = 1, 
        lty = "solid", col = "orange", yaxt = 'n', ylim = c(0, 19),
        xlim = c(0.8, 3.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "ENCODE ATAC-seq signal", xlab = "")
boxplot(vec2, frame = FALSE, yaxt = 'n',
        lty = "solid", col = "deep pink", medlty = 1, medlwd = 0.8, at = 2, xaxt = 'n', add = TRUE,
        boxlty = 0, staplelwd = 0, outline = FALSE)
boxplot(vec3, frame = FALSE, yaxt = 'n',
        lty = "solid", col = rgb(0,0,0,0.75), medlty = 1, medlwd = 0.8, at = 3, xaxt = 'n', add = TRUE,
        boxlty = 0, staplelwd = 0, outline = FALSE)

axis(2, at = c(0, 7.5, 15))
dev.off()

quartz(file = "ks2_encode_overlap.pdf", height = 2.4, width = 2, pointsize = 8, type = "pdf")
load(file = "ks2_granges_encode_overlap.rda")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
boxplot(vec1, frame = FALSE, main = "ENCODE vs KS2", font.main = 1, 
        lty = "solid", col = "orange", yaxt = 'n', ylim = c(0, 19),
        xlim = c(0.8, 3.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "ENCODE ATAC-seq signal", xlab = "")
boxplot(vec2, frame = FALSE, yaxt = 'n',
        lty = "solid", col = "deep pink", medlty = 1, medlwd = 0.8, at = 2, xaxt = 'n', add = TRUE,
        boxlty = 0, staplelwd = 0, outline = FALSE)
boxplot(vec3, frame = FALSE, yaxt = 'n',
        lty = "solid", col = rgb(0,0,0,0.75), medlty = 1, medlwd = 0.8, at = 3, xaxt = 'n', add = TRUE,
        boxlty = 0, staplelwd = 0, outline = FALSE)

axis(2, at = c(0, 7.5, 15))
dev.off()

###for T cells
setwd('/dcl01/hansen/data/kabuki_mice/HiSeq325_ATAC_new/ATAC/dedup_files/merged_files/encode/T_cells')
encode_T_1 <- import('ENCFF297OQF.bed', format = "narrowPeak")
encode_T_2 <- import('ENCFF160JGQ.bed', format = "narrowPeak")
encode_T_3 <- import('ENCFF119WHO.bed', format = "narrowPeak")
encode_T_4 <- import('ENCFF605SFN.bed', format = "narrowPeak")
combined_encode_T <- reduce(unlist(GRangesList(encode_T_1, encode_T_2, encode_T_3, encode_T_4)))

setwd('/dcl01/hansen/data/kabuki_mice/HiSeq325_ATAC_new/ATAC/dedup_files/merged_files')
wt_granges <- getGrangesFromNarrowPeaks('all_wt_samples_T_peaks.narrowPeak')
ks1_granges <- getGrangesFromNarrowPeaks('all_KS1_samples_T_peaks.narrowPeak')
ks2_granges <- getGrangesFromNarrowPeaks('all_KS2_samples_T_peaks.narrowPeak')

length(unique(queryHits(findOverlaps(wt_granges, combined_encode_T))))

