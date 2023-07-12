###supp tables for the genes downstream of promoters with disrupted accessibility in 
###all three cell types: neurons, B and T cells
res_granges$qvalue <- qvalue(res_granges$pvalue, fdr.level = 0.1)$qvalues
res2_granges$qvalue <- qvalue(res2_granges$pvalue, fdr.level = 0.1)$qvalues

atac_B_KS1_granges$qvalue <- qvalue(atac_B_KS1_granges$pvalue, fdr.level = 0.1)$qvalues
atac_T_KS1_granges$qvalue <- qvalue(atac_T_KS1_granges$pvalue, fdr.level = 0.1)$qvalues

atac_B_KS2_granges$qvalue <- qvalue(atac_B_KS2_granges$pvalue, fdr.level = 0.1)$qvalues
atac_T_KS2_granges$qvalue <- qvalue(atac_T_KS2_granges$pvalue, fdr.level = 0.1)$qvalues



##KS1
res_granges$log2FoldChange <- -res_granges$log2FoldChange
atac_B_KS1_granges$log2FoldChange <- -atac_B_KS1_granges$log2FoldChange
atac_T_KS1_granges$log2FoldChange <- -atac_T_KS1_granges$log2FoldChange

T_shared_KS1 <- KS1_neurons_B_T[-which(is.na(KS1_neurons_B_T$prom_overlapping_id))]
T_shared_KS1$log2FoldChange <- -T_shared_KS1$log2FoldChange
shared_genes_df_KS1_T <- annoGR2DF(T_shared_KS1)
shared_genes_df_KS1_T <- shared_genes_df_KS1_T[, c("chr", "start", "end", "width", "log2FoldChange", "pvalue", "qvalue")]

overlaps <- findOverlaps(res_granges, T_shared_KS1)
neurons_shared <- res_granges[unique(queryHits(overlaps))]
shared_genes_df_KS1_neurons <- annoGR2DF(neurons_shared)
shared_genes_df_KS1_neurons <- shared_genes_df_KS1_neurons[, c("chr", "start", "end", "width", "log2FoldChange", "pvalue", "qvalue")]

overlaps <- findOverlaps(atac_B_KS1_granges, T_shared_KS1)
B_shared <- atac_B_KS1_granges[unique(queryHits(overlaps))]
shared_genes_df_KS1_B <- annoGR2DF(B_shared)
shared_genes_df_KS1_B <- shared_genes_df_KS1_B[, c("chr", "start", "end", "width", "log2FoldChange", "pvalue", "qvalue")]

KS1_names <- annoGR2DF(T_shared_KS1)
KS1_names <- KS1_names[, c("prom_overlapping_id", "prom_overlapping_name")]
KS1_names <- KS1_names[-which(duplicated(KS1_names$prom_overlapping_id)), ]
KS1_names <- data.frame(gene_id = unlist(KS1_names$prom_overlapping_id), 
                        gene_name = unlist(KS1_names$prom_overlapping_name))

#csv files
write_csv(KS1_names, file = "atac_shared_genes_KS1.csv")
write_csv(shared_genes_df_KS1_T, file = "atac_shared_locations_KS1_T_coordinates.csv")
write_csv(shared_genes_df_KS1_B, file = "atac_shared_locations_KS1_B_coordinates.csv")
write_csv(shared_genes_df_KS1_neurons, file = "atac_shared_locations_KS1_neurons_coordinates.csv")


##KS2
res2_granges$log2FoldChange <- -res2_granges$log2FoldChange
atac_B_KS2_granges$log2FoldChange <- -atac_B_KS2_granges$log2FoldChange
atac_T_KS2_granges$log2FoldChange <- -atac_T_KS2_granges$log2FoldChange

T_shared_KS2 <- KS2_neurons_B_T[-which(is.na(KS2_neurons_B_T$prom_overlapping_id))]
T_shared_KS2$log2FoldChange <- -T_shared_KS2$log2FoldChange
shared_genes_df_KS2_T <- annoGR2DF(T_shared_KS2)
shared_genes_df_KS2_T <- shared_genes_df_KS2_T[, c("chr", "start", "end", "width", "log2FoldChange", "pvalue", "qvalue")]

overlaps <- findOverlaps(res_granges, T_shared_KS2)
neurons_shared <- res_granges[unique(queryHits(overlaps))]
shared_genes_df_KS2_neurons <- annoGR2DF(neurons_shared)
shared_genes_df_KS2_neurons <- shared_genes_df_KS2_neurons[, c("chr", "start", "end", "width", "log2FoldChange", "pvalue", "qvalue")]

overlaps <- findOverlaps(atac_B_KS2_granges, T_shared_KS2)
B_shared <- atac_B_KS2_granges[unique(queryHits(overlaps))]
shared_genes_df_KS2_B <- annoGR2DF(B_shared)
shared_genes_df_KS2_B <- shared_genes_df_KS2_B[, c("chr", "start", "end", "width", "log2FoldChange", "pvalue", "qvalue")]

KS2_names <- annoGR2DF(T_shared_KS2)
KS2_names <- KS2_names[, c("prom_overlapping_id", "prom_overlapping_name")]
KS2_names <- KS2_names[-which(duplicated(KS2_names$prom_overlapping_id)), ]
KS2_names <- data.frame(gene_id = unlist(KS2_names$prom_overlapping_id), 
                        gene_name = unlist(KS2_names$prom_overlapping_name))

#csv files
write_csv(KS2_names, file = "atac_shared_genes_KS2.csv")
write_csv(shared_genes_df_KS2_T, file = "atac_shared_locations_KS2_T_coordinates.csv")
write_csv(shared_genes_df_KS2_B, file = "atac_shared_locations_KS2_B_coordinates.csv")
write_csv(shared_genes_df_KS2_neurons, file = "atac_shared_locations_KS2_neurons_coordinates.csv")

