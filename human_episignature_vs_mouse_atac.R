##import the episignature differentially methylated positions and covert to granges
episignature_ks1 <- read_csv('~/Downloads/human_episignature_aref_eshgi_et_al.csv') ##from Aref Eshgi et al. 2017
episignature_genes <- episignature_ks1[-which(is.na(episignature_ks1$Gene)), ]
colnames(episignature_genes)[3] <- "Start"
episignature_genes$End <- episignature_genes$Start
episignature_genes <- episignature_genes[grep("chr", episignature_genes$Chr), ]
episignature_genes_granges <- makeGRangesFromDataFrame(episignature_genes, keep.extra.columns = TRUE)

##now get human promoter coordinates (promoter = +/- 500bp from TSS)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)
edb <- EnsDb.Hsapiens.v75
proms_all <- promoters(edb, upstream = 500, downstream = 500, 
                       columns = c("gene_name", "tx_id", "tx_cds_seq_start", "tx_cds_seq_end",
                                   "tx_biotype", "gene_id"))
genome(seqinfo(proms_all)) <- "hg19"
seqlevelsStyle(proms_all) <- "ucsc"

chrs <- names(Hsapiens)[1:24]
proms_all <- proms_all[which(seqnames(proms_all) %in% chrs[1:24])] 

##get gene IDs of promoters that overlap episignature differentially methylated positions
overlaps <- findOverlaps(proms_all, episignature_genes_granges)
human_episignature_gene_ids <- unique(proms_all$gene_id[unique(queryHits(overlaps))])


##covert gene ensembl IDs to mouse ensembl IDs
library(biomaRt)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human_mouse_homologs <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                             "mmusculus_homolog_ensembl_gene", 
                                             "mmusculus_homolog_associated_gene_name", 
                                             "mmusculus_homolog_orthology_confidence", 
                                             "mmusculus_homolog_perc_id_r1"),
                              filters = "ensembl_gene_id",
                              values = human_episignature_gene_ids,
                              mart = human)
human_mouse_homologs <- human_mouse_homologs[which(human_mouse_homologs$mmusculus_homolog_orthology_confidence == 1), ]
mouse_episignature_gene_ids <- unique(human_mouse_homologs$mmusculus_homolog_ensembl_gene)

##now do the statistics (the gene level results are obtained from the atac_pathways_neurons.R script)
indices_episignature_all_dmps <- which(gene_level_results_neurons_KS1$gene_ids %in% 
                                  mouse_episignature_gene_ids)
indices_episignature_B_all_dmps <- which(gene_level_results_B_KS1$gene_ids %in% 
                                    mouse_episignature_gene_ids)
indices_episignature_T_all_dmps <- which(gene_level_results_T_KS1$gene_ids %in% 
                                    mouse_episignature_gene_ids)

#
wilcox_stat_neurons_KS1_all_dmps <- wilcox.test(gene_level_results_neurons_KS1$pval[indices_episignature_all_dmps], 
                                         gene_level_results_neurons_KS1$pval[-indices_episignature_all_dmps])$statistic
null_distribution_neurons_KS1_all_dmps <- replicate(1000, {
  random_indices <- sample(1:length(gene_level_results_neurons_KS1$pval), length(indices_episignature_all_dmps))
  random_wilcox_stat <- wilcox.test(gene_level_results_neurons_KS1$pval[random_indices], 
                                    gene_level_results_neurons_KS1$pval[-random_indices])$statistic
  random_wilcox_stat
})
#
wilcox_stat_B_KS1_all_dmps <- wilcox.test(gene_level_results_B_KS1$pval[indices_episignature_B_all_dmps], 
                                   gene_level_results_B_KS1$pval[-indices_episignature_B_all_dmps])$statistic
null_distribution_B_KS1_all_dmps <- replicate(1000, {
  random_indices <- sample(1:length(gene_level_results_B_KS1$pval), length(indices_episignature_B_all_dmps))
  random_wilcox_stat <- wilcox.test(gene_level_results_B_KS1$pval[random_indices], 
                                    gene_level_results_B_KS1$pval[-random_indices])$statistic
  random_wilcox_stat
})
#
wilcox_stat_T_KS1_all_dmps <- wilcox.test(gene_level_results_T_KS1$pval[indices_episignature_T_all_dmps], 
                                   gene_level_results_T_KS1$pval[-indices_episignature_T_all_dmps])$statistic
null_distribution_T_KS1_all_dmps <- replicate(1000, {
  random_indices <- sample(1:length(gene_level_results_T_KS1$pval), length(indices_episignature_T_all_dmps))
  random_wilcox_stat <- wilcox.test(gene_level_results_T_KS1$pval[random_indices], 
                                    gene_level_results_T_KS1$pval[-random_indices])$statistic
  random_wilcox_stat
})

quartz(file = "episignature_ranks_KS1_all.pdf", width = 6, height = 2, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
makeWilcoxPlot(null_distribution_neurons_KS1_all_dmps, wilcox_stat_neurons_KS1_all_dmps, 
               main_lab = "KS1 neurons")
legend <- legend("topright", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
                 cex = 1, lty = "solid", lwd = 2.5)

makeWilcoxPlot(null_distribution_B_KS1_all_dmps, wilcox_stat_B_KS1_all_dmps, 
               main_lab = "KS1 B cells")
makeWilcoxPlot(null_distribution_T_KS1_all_dmps, wilcox_stat_T_KS1_all_dmps, 
               main_lab = "KS1 T cells")
dev.off()














