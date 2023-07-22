###
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)
load(file = "ezh2_granges_hg19.rda") ##contains binding locations of EZH2 as identified by chip-seq experiments from the encode project

ezh2_all <- ezh2
ezh2_moderate_evicence <- ezh2[which(ezh2$sourceCount >= 2)]
ezh2_strong_evicence <- ezh2[which(ezh2$sourceCount >= 5)]

###find overlaps of ezh2 binding locations with promoters to obtain "ezh2 target genes"
edb <- EnsDb.Hsapiens.v75
proms_all <- promoters(edb, upstream = 2000, downstream = 2000, 
                       columns = c("gene_name", "tx_id", "tx_cds_seq_start", "tx_cds_seq_end",
                                   "tx_biotype", "gene_id"))
genome(seqinfo(proms_all)) <- "hg19"
seqlevelsStyle(proms_all) <- "ucsc"

chrs <- names(Hsapiens)[1:24]
proms_all <- proms_all[which(seqnames(proms_all) %in% chrs[1:24])] 

ezh2_target_ids_all <- unique(proms_all$gene_id[queryHits(findOverlaps(proms_all, ezh2_all))])
ezh2_target_ids_moderate <- unique(proms_all$gene_id[queryHits(findOverlaps(proms_all, ezh2_moderate_evicence))])
ezh2_target_ids_strong <- unique(proms_all$gene_id[queryHits(findOverlaps(proms_all, ezh2_strong_evicence))])


###get mouse homologs
library(biomaRt)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

human_mouse_homologs <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                             "mmusculus_homolog_ensembl_gene", 
                                             "mmusculus_homolog_associated_gene_name", 
                                             "mmusculus_homolog_orthology_confidence", 
                                             "mmusculus_homolog_perc_id_r1"),
                              filters = "ensembl_gene_id",
                              values = ezh2_target_ids_all,
                              mart = human)
human_mouse_homologs <- human_mouse_homologs[which(human_mouse_homologs$mmusculus_homolog_orthology_confidence == 1), ]
mouse_ezh2_targets_all <- unique(human_mouse_homologs$mmusculus_homolog_ensembl_gene)

human_mouse_homologs <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                             "mmusculus_homolog_ensembl_gene", 
                                             "mmusculus_homolog_associated_gene_name", 
                                             "mmusculus_homolog_orthology_confidence", 
                                             "mmusculus_homolog_perc_id_r1"),
                              filters = "ensembl_gene_id",
                              values = ezh2_target_ids_strong,
                              mart = human)
human_mouse_homologs <- human_mouse_homologs[which(human_mouse_homologs$mmusculus_homolog_orthology_confidence == 1), ]
mouse_ezh2_targets_strong <- unique(human_mouse_homologs$mmusculus_homolog_ensembl_gene)










