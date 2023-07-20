load(file = "atac_differential_results_neurons_B_T_july2023.rda")

####
prom_overlaps_B <- findOverlaps(atac_B_KS1_granges, proms_mouse)
atac_B_KS1_granges_proms <- atac_B_KS1_granges[queryHits(prom_overlaps_B)]
atac_B_KS1_granges_proms$gene_id <- proms_mouse$gene_id[subjectHits(prom_overlaps_B)]
atac_B_KS1_granges_proms$gene_name <- proms_mouse$gene_name[subjectHits(prom_overlaps_B)]
atac_B_KS1_granges_proms <- atac_B_KS1_granges_proms[-which(duplicated(atac_B_KS1_granges_proms))]

prom_overlaps2_B <- findOverlaps(atac_B_KS2_granges, proms_mouse)
atac_B_KS2_granges_proms <- atac_B_KS2_granges[queryHits(prom_overlaps2_B)]
atac_B_KS2_granges_proms$gene_id <- proms_mouse$gene_id[subjectHits(prom_overlaps2_B)]
atac_B_KS2_granges_proms$gene_name <- proms_mouse$gene_name[subjectHits(prom_overlaps2_B)]
atac_B_KS2_granges_proms <- atac_B_KS2_granges_proms[-which(duplicated(atac_B_KS2_granges_proms))]

prom_overlaps_T <- findOverlaps(atac_T_KS1_granges, proms_mouse)
atac_T_KS1_granges_proms <- atac_T_KS1_granges[queryHits(prom_overlaps_T)]
atac_T_KS1_granges_proms$gene_id <- proms_mouse$gene_id[subjectHits(prom_overlaps_T)]
atac_T_KS1_granges_proms$gene_name <- proms_mouse$gene_name[subjectHits(prom_overlaps_T)]
atac_T_KS1_granges_proms <- atac_T_KS1_granges_proms[-which(duplicated(atac_T_KS1_granges_proms))]

prom_overlaps2_T <- findOverlaps(atac_T_KS2_granges, proms_mouse)
atac_T_KS2_granges_proms <- atac_T_KS2_granges[queryHits(prom_overlaps2_T)]
atac_T_KS2_granges_proms$gene_id <- proms_mouse$gene_id[subjectHits(prom_overlaps2_T)]
atac_T_KS2_granges_proms$gene_name <- proms_mouse$gene_name[subjectHits(prom_overlaps2_T)]
atac_T_KS2_granges_proms <- atac_T_KS2_granges_proms[-which(duplicated(atac_T_KS2_granges_proms))]
