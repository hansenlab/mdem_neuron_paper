library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
edb_mouse <- EnsDb.Mmusculus.v79
proms_mouse <- promoters(edb_mouse, filter = TxBiotypeFilter("protein_coding"), upstream = 2000, downstream = 2000, columns = c("gene_name", "tx_id", "gene_id"))
#genome(seqinfo(proms_mouse)) <- "mm10"
genome(seqinfo(proms_mouse)) <- "mm10"
#seqlevelsStyle(proms_mouse) <- "ucsc"
seqlevelsStyle(proms_mouse) <- "ucsc"
#proms_mouse <- proms_mouse[-which(duplicated(proms_mouse$gene_id))]
proms_mouse <- proms_mouse[which(seqnames(proms_mouse) %in% seqnames(Mmusculus)[1:21])]