###neurons
#separately get the KS1 and KS2 matrices
mat <- peaks_by_samples[, -which(colnames(peaks_by_samples) == "Cohort9_2-22_wt")]
sampinf <- sample_info[-which(colnames(peaks_by_samples) == "Cohort9_2-22_wt"), ]
dds <- getDDSObject(mat, sampinf, "KS1", "9")

dds2 <- getDDSObject(peaks_by_samples, sample_info, "KS2", c("1", "8"))

mat1 <- mat[, which(sampinf$genotype %in% c("WT", "KS1") & 
                      sampinf$cohort %in% c("9"))]
mat2 <- peaks_by_samples[, which(sample_info$genotype %in% c("WT", "KS2") & 
                                   sample_info$cohort %in% c("1", "8"))]

#combine the KS1 and KS2 matrices
geo_mat_neurons <- cbind(mat1, mat2)

geo_mat_neurons_df <- as.data.frame(geo_mat_neurons)
library(tibble)
geo_mat_neurons_df <- rownames_to_column(geo_mat_neurons_df, "mm10_coordinates")
write_csv(geo_mat_neurons_df, "GEO/neurons_paper_2023/atac_peaks_by_samples_matrix_neurons.csv")

###T cells
#get the KS1 and KS2 matrices combined
geo_mat_T <- peaks_by_samples[, grep("_T_", colnames(peaks_by_samples))]
geo_mat_T_df <- as.data.frame(geo_mat_T)
geo_mat_T_df <- rownames_to_column(geo_mat_T_df, "mm10_coordinates")
write_csv(geo_mat_T_df, "GEO/neurons_paper_2023/atac_peaks_by_samples_matrix_T.csv")

