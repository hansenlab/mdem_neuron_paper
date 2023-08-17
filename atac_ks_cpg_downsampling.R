setwd('/dcl01/hansen/data/kabuki_mice/HiSeq325_ATAC_new/ATAC/dedup_files/merged_files')
bamFiles_filepaths <- list.files(pattern = "noMito_merged_downsampled_25_percent_with_samtools.bam")
load(file = "combined_granges_from_peaks_called_from_condition_specific_bams_with_all_fragments_B.rda")
#For T cells instead
load(file = "combined_granges_from_peaks_called_from_condition_specific_bams_with_all_fragments_T.rda")

count_reads_in_peaks <- featureCounts(bamFiles_filepaths, annot.ext = combined_peaks_ann, 
                                      isPairedEnd = TRUE, minOverlap = 3,
                                      requireBothEndsMapped = TRUE, countChimericFragments = FALSE, 
                                      nthreads = 30, countMultiMappingReads = FALSE)
save(count_reads_in_peaks, 
     file = "DOWNSAMPLED25percent_count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_kabuki_cohort.rda")


#For T cells instead
#save(count_reads_in_peaks, 
#     file = "DOWNSAMPLED25percent_count_reads_in_combined_peaks_from_condition_specific_bam_files_T_cells_only_kabuki_cohort.rda")



#####
load(file = "rt_atacseq/DOWNSAMPLED25percent_count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_kabuki_cohort.rda")
peaks_by_samples <- getPeaksBySamplesMatrix("rt_atacseq/DOWNSAMPLED25percent_count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_kabuki_cohort.rda", 
                                            "blood_all_samples")

#For T cells instead
load(file = "rt_atacseq/DOWNSAMPLED25percent_count_reads_in_combined_peaks_from_condition_specific_bam_files_T_cells_only_kabuki_cohort.rda")
peaks_by_samples <- getPeaksBySamplesMatrix("rt_atacseq/DOWNSAMPLED25percent_count_reads_in_combined_peaks_from_condition_specific_bam_files_T_cells_only_kabuki_cohort.rda", 
                                            "blood_all_samples")


###
downsampled_atac_B_KS1_granges_proms <- getPromoterResultsRanges(downsampled_atac_B_KS1_granges)
downsampled_atac_B_KS2_granges_proms <- getPromoterResultsRanges(downsampled_atac_B_KS2_granges)
downsampled_atac_T_KS1_granges_proms <- getPromoterResultsRanges(downsampled_atac_T_KS1_granges)
downsampled_atac_T_KS2_granges_proms <- getPromoterResultsRanges(downsampled_atac_T_KS2_granges)

##
downsampled_atac_B_KS1_granges_proms$overlaps_cgi <- "no"
downsampled_atac_B_KS1_granges_proms$overlaps_cgi[unique(queryHits(findOverlaps(
  downsampled_atac_B_KS1_granges_proms, cpg)))] <- "yes"

downsampled_atac_T_KS1_granges_proms$overlaps_cgi <- "no"
downsampled_atac_T_KS1_granges_proms$overlaps_cgi[unique(queryHits(findOverlaps(
  downsampled_atac_T_KS1_granges_proms, cpg)))] <- "yes"

downsampled_atac_B_KS2_granges_proms$overlaps_cgi <- "no"
downsampled_atac_B_KS2_granges_proms$overlaps_cgi[unique(queryHits(findOverlaps(
  downsampled_atac_B_KS2_granges_proms, cpg)))] <- "yes"

downsampled_atac_T_KS2_granges_proms$overlaps_cgi <- "no"
downsampled_atac_T_KS2_granges_proms$overlaps_cgi[unique(queryHits(findOverlaps(
  downsampled_atac_T_KS2_granges_proms, cpg)))] <- "yes"

##
atac_B_KS1_granges_proms$overlaps_cgi <- "no"
atac_B_KS1_granges_proms$overlaps_cgi[unique(queryHits(findOverlaps(
  atac_B_KS1_granges_proms, cpg)))] <- "yes"

atac_T_KS1_granges_proms$overlaps_cgi <- "no"
atac_T_KS1_granges_proms$overlaps_cgi[unique(queryHits(findOverlaps(
  atac_T_KS1_granges_proms, cpg)))] <- "yes"

atac_B_KS2_granges_proms$overlaps_cgi <- "no"
atac_B_KS2_granges_proms$overlaps_cgi[unique(queryHits(findOverlaps(
  atac_B_KS2_granges_proms, cpg)))] <- "yes"

atac_T_KS2_granges_proms$overlaps_cgi <- "no"
atac_T_KS2_granges_proms$overlaps_cgi[unique(queryHits(findOverlaps(
  atac_T_KS2_granges_proms, cpg)))] <- "yes"

##
res_proms$overlaps_cgi <- "no"
res_proms$overlaps_cgi[unique(queryHits(findOverlaps(
  res_proms, cpg)))] <- "yes"

res2_proms$overlaps_cgi <- "no"
res2_proms$overlaps_cgi[unique(queryHits(findOverlaps(
  res2_proms, cpg)))] <- "yes"


getDevianceToNullDevianceRatio <- function(results_granges){
  model <- glm(as.factor(results_granges$overlaps_cgi) ~ results_granges$pvalue, family = "binomial")
  1 - model$deviance/model$null.deviance
}

getDevianceToNullDevianceRatio(atac_B_KS1_granges_proms)
getDevianceToNullDevianceRatio(atac_T_KS1_granges_proms)
getDevianceToNullDevianceRatio(downsampled_atac_B_KS1_granges_proms)
getDevianceToNullDevianceRatio(downsampled_atac_T_KS1_granges_proms)

getDevianceToNullDevianceRatio(atac_B_KS2_granges_proms)
getDevianceToNullDevianceRatio(atac_T_KS2_granges_proms)
getDevianceToNullDevianceRatio(downsampled_atac_B_KS2_granges_proms)
getDevianceToNullDevianceRatio(downsampled_atac_T_KS2_granges_proms)



