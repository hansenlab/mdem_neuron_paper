###supp tables for the genes downstream of promoters with disrupted accessibility in 
###all three cell types: neurons, B and T cells
res_granges$qvalue <- qvalue(res_granges$pvalue, fdr.level = 0.1)$qvalues
res2_granges$qvalue <- qvalue(res2_granges$pvalue, fdr.level = 0.1)$qvalues

atac_B_KS1_granges$qvalue <- qvalue(atac_B_KS1_granges$pvalue, fdr.level = 0.1)$qvalues
atac_T_KS1_granges$qvalue <- qvalue(atac_T_KS1_granges$pvalue, fdr.level = 0.1)$qvalues

atac_B_KS2_granges$qvalue <- qvalue(atac_B_KS2_granges$pvalue, fdr.level = 0.1)$qvalues
atac_T_KS2_granges$qvalue <- qvalue(atac_T_KS2_granges$pvalue, fdr.level = 0.1)$qvalues




###shared gene ensembl IDs and gene names
write_csv(shared_genes_df[, 1:2], file = "atac_shared_genes_KS1.csv")
write_csv(shared_genes2_df[, 1:2], file = "atac_shared_genes_KS2.csv")


###
addOverlappingGenes <- function(genomic_ranges){
  genomic_ranges$gene_id <- NA
  genomic_ranges$gene_name <- NA
  overlaps <- findOverlaps(genomic_ranges, proms_mouse)
  results_proms <- genomic_ranges[queryHits(overlaps)]
  results_proms$gene_id <- proms_mouse$gene_id[subjectHits(overlaps)]
  results_proms$gene_name <- proms_mouse$gene_name[subjectHits(overlaps)]
  
  gene_level <- split(results_proms, results_proms$gene_name)
  gene_level <- endoapply(gene_level, function(xx) unique(xx))
  c(unlist(gene_level), genomic_ranges[-queryHits(overlaps)])
}

mergeDuplicates <- function(genomic_ranges){
  dup_ind <- which(duplicated(genomic_ranges))
  indices_drop <- vector()
  for (i in dup_ind){
    indices <- which(genomic_ranges == genomic_ranges[i])
    gene_ids_same_range <- paste(genomic_ranges[indices]$gene_id, collapse = ";")
    gene_names_same_range <- paste(genomic_ranges[indices]$gene_name, collapse = ";")
    genomic_ranges$gene_id[indices[1]] <- gene_ids_same_range # Assign to the first occurrence
    genomic_ranges$gene_name[indices[1]] <- gene_names_same_range # Assign to the first occurrence
    
    indices_drop <- c(indices_drop, indices[-1])
    genomic_ranges
    indices_drop
  }
  genomic_ranges[-indices_drop]
}

##KS1
T_shared_KS1 <- KS1_neurons_B_T
T_shared_KS1$log2FoldChange <- -T_shared_KS1$log2FoldChange
T_shared_KS1 <- addOverlappingGenes(T_shared_KS1)
T_shared_KS1 <- mergeDuplicates(T_shared_KS1)
shared_genes_df_KS1_T <- annoGR2DF(T_shared_KS1)
shared_genes_df_KS1_T <- shared_genes_df_KS1_T[, c("chr", "start", "end", "width", 
                                                   "log2FoldChange", "pvalue", "qvalue", "gene_id", "gene_name")]

overlaps <- findOverlaps(res_granges, T_shared_KS1)
neurons_shared_KS1 <- res_granges[unique(queryHits(overlaps))]
neurons_shared_KS1$log2FoldChange <- -neurons_shared_KS1$log2FoldChange
neurons_shared_KS1 <- addOverlappingGenes(neurons_shared_KS1)
neurons_shared_KS1 <- mergeDuplicates(neurons_shared_KS1)
shared_genes_df_KS1_neurons <- annoGR2DF(neurons_shared_KS1)
shared_genes_df_KS1_neurons <- shared_genes_df_KS1_neurons[, c("chr", "start", "end", "width", 
                                                               "log2FoldChange", "pvalue", "qvalue", "gene_id", "gene_name")]

overlaps <- findOverlaps(atac_B_KS1_granges, T_shared_KS1)
B_shared_KS1 <- atac_B_KS1_granges[unique(queryHits(overlaps))]
B_shared_KS1$log2FoldChange <- -B_shared_KS1$log2FoldChange
B_shared_KS1 <- addOverlappingGenes(B_shared_KS1)
B_shared_KS1 <- mergeDuplicates(B_shared_KS1)
shared_genes_df_KS1_B <- annoGR2DF(B_shared_KS1)
shared_genes_df_KS1_B <- shared_genes_df_KS1_B[, c("chr", "start", "end", "width", 
                                                   "log2FoldChange", "pvalue", "qvalue", "gene_id", "gene_name")]

#KS1_names <- annoGR2DF(T_shared_KS1)
#KS1_names <- KS1_names[, c("prom_overlapping_id", "prom_overlapping_name")]
#KS1_names <- KS1_names[-which(duplicated(KS1_names$prom_overlapping_id)), ]
#KS1_names <- data.frame(gene_id = unlist(KS1_names$prom_overlapping_id), 
#                        gene_name = unlist(KS1_names$prom_overlapping_name))

#csv files
write_csv(shared_genes_df_KS1_T, file = "neuron_mdem_biorxiv_1/supp_tables/SuppTable8_atac_shared_locations_KS1_T_coordinates.csv")
write_csv(shared_genes_df_KS1_B, file = "neuron_mdem_biorxiv_1/supp_tables/SuppTable7_atac_shared_locations_KS1_B_coordinates.csv")
write_csv(shared_genes_df_KS1_neurons, file = "neuron_mdem_biorxiv_1/supp_tables/SuppTable9_atac_shared_locations_KS1_neurons_coordinates.csv")


##KS2
T_shared_KS2 <- KS2_neurons_B_T
T_shared_KS2$log2FoldChange <- -T_shared_KS2$log2FoldChange
T_shared_KS2 <- addOverlappingGenes(T_shared_KS2)
T_shared_KS2 <- mergeDuplicates(T_shared_KS2)
shared_genes_df_KS2_T <- annoGR2DF(T_shared_KS2)
shared_genes_df_KS2_T <- shared_genes_df_KS2_T[, c("chr", "start", "end", "width", 
                                                   "log2FoldChange", "pvalue", "qvalue", "gene_id", "gene_name")]

overlaps <- findOverlaps(res_granges, T_shared_KS2)
neurons_shared_KS2 <- res_granges[unique(queryHits(overlaps))]
neurons_shared_KS2$log2FoldChange <- -neurons_shared_KS2$log2FoldChange
neurons_shared_KS2 <- addOverlappingGenes(neurons_shared_KS2)
neurons_shared_KS2 <- mergeDuplicates(neurons_shared_KS2)
shared_genes_df_KS2_neurons <- annoGR2DF(neurons_shared_KS2)
shared_genes_df_KS2_neurons <- shared_genes_df_KS2_neurons[, c("chr", "start", "end", "width", 
                                                               "log2FoldChange", "pvalue", "qvalue", "gene_id", "gene_name")]

overlaps <- findOverlaps(atac_B_KS2_granges, T_shared_KS2)
B_shared_KS2 <- atac_B_KS2_granges[unique(queryHits(overlaps))]
B_shared_KS2$log2FoldChange <- -B_shared_KS2$log2FoldChange
B_shared_KS2 <- addOverlappingGenes(B_shared_KS2)
B_shared_KS2 <- mergeDuplicates(B_shared_KS2)
shared_genes_df_KS2_B <- annoGR2DF(B_shared_KS2)
shared_genes_df_KS2_B <- shared_genes_df_KS2_B[, c("chr", "start", "end", "width", 
                                                   "log2FoldChange", "pvalue", "qvalue", "gene_id", "gene_name")]

#KS2_names <- annoGR2DF(T_shared_KS2)
#KS2_names <- KS2_names[, c("prom_overlapping_id", "prom_overlapping_name")]
#KS2_names <- data.frame(gene_id = unlist(KS2_names$prom_overlapping_id), 
#                        gene_name = unlist(KS2_names$prom_overlapping_name))

#csv files
write_csv(shared_genes_df_KS2_T, file = "neuron_mdem_biorxiv_1/supp_tables/SuppTable12_atac_shared_locations_KS2_T_coordinates.csv")
write_csv(shared_genes_df_KS2_B, file = "neuron_mdem_biorxiv_1/supp_tables/SuppTable11_atac_shared_locations_KS2_B_coordinates.csv")
write_csv(shared_genes_df_KS2_neurons, file = "neuron_mdem_biorxiv_1/supp_tables/SuppTable13_atac_shared_locations_KS2_neurons_coordinates.csv")

