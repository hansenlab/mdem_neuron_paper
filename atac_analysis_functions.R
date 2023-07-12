####
info_df <- read_csv('/Users/leandrosboukas/Desktop/KS1_KS2_analysis/brain_samples_identifiers.csv') # this is necessary for loading the brain samples
info_df$sample_name <- gsub("[.]", "-", info_df$sample_name)
  
getPeaksBySamplesMatrix <- function(file_path, what = c("neurons", "blood_all_samples", 
                                                        "blood_first_run", "blood_second_run", "rt_blood")){
  load(file = file_path)
  peaks_by_samples <- count_reads_in_peaks$counts
  
  if (what == "neurons"){
    
    for (i in 1:length(colnames(peaks_by_samples))){
      index_to_use <- grep(paste0(info_df$sample_name[i], "_"), colnames(peaks_by_samples))
      colnames(peaks_by_samples)[index_to_use] <- paste0(info_df$sample_name[i], "_", info_df$genotype[i])
      peaks_by_samples
    }
    colnames(peaks_by_samples)[19:25] <- gsub("Cohort9.2", "Cohort9_2", colnames(peaks_by_samples)[19:25])
  } else if (what == "blood_all_samples"){
    colnames(peaks_by_samples) <- c("100_B_KDM6A", "100_T_KDM6A", "101_B_WT", "101_T_WT", 
                                    "102_B_KDM6A", "102_T_KDM6A", "103_B_KDM6A", "103_T_KDM6A",
                                    "104_B_WT", "104_T_WT", "106_B_KDM6A", "106_T_KDM6A",
                                    "21_B_WT", "21_T_WT", "23_B_KMT2D", "23_T_KMT2D", 
                                    "24_B_WT", "24_T_WT", "26_B_WT", "26_T_WT", 
                                    "27_B_WT", "27_T_WT", "27_B_KMT2D", "27_T_KMT2D",
                                    "29_B_KMT2D", "29_T_KMT2D", "31_B_KMT2D", "31_T_KMT2D", 
                                    "32_B_WT", "32_T_WT", "37_B_WT", "37_T_WT", 
                                    "40_T_WT", "42_B_KDM6A", "42_T_KDM6A", "49_B_WT", 
                                    "49_T_WT", "63_B_WT", "63_T_WT", "65_B_WT", "65_T_WT", 
                                    "6_B_KMT2D", "6_T_KMT2D", "79_B_KMT2D", "79_T_KMT2D", 
                                    "84_B_KMT2D", "84_T_KMT2D", "40_B_WT")
  } else if (what == "blood_second_run"){
    colnames(peaks_by_samples) <- c("21_B_WT", "21_T_WT",      ####this is only for the blood samples from the new seq run
                                     "24_B_KMT2D", "24_T_KMT2D", 
                                     "26_B_WT", "26_T_WT", 
                                     "27_B_WT", "27_T_WT", 
                                     "29_B_KMT2D", "29_T_KMT2D", 
                                     "31_B_KMT2D", "31_T_KMT2D", 
                                     "37_B_WT", "37_T_WT", 
                                     "49_B_KDM6A", "49_T_KDM6A",
                                     "65_B_WT", "65_T_WT")
  } else if (what == "blood_first_run"){
    colnames(peaks_by_samples) <- c("100_B_KDM6A", "100_T_KDM6A", "101_B_WT", "101_T_WT", 
                                     "102_B_KDM6A", "102_T_KDM6A", "103_B_WT", "103_T_WT",
                                     "104_B_KDM6A", "104_T_KDM6A", "106_B_KDM6A", "106_T_KDM6A",
                                     "23_B_KMT2D", "23_T_KMT2D", "27_B_KMT2D", "27_T_KMT2D",
                                     "32_B_WT", "32_T_WT", "40_T_WT", "42_B_KDM6A", "42_T_KDM6A",
                                     "63_B_WT", "63_T_WT", "6_B_KMT2D", "6_T_KMT2D", 
                                     "79_B_KMT2D", "79_T_KMT2D", "84_B_KMT2D", "84_T_KMT2D",
                                     "40_B_WT")
  } else if (what == "rt_blood"){
    cbp_atac_info <- read_csv('cbp_info_atac.csv')
    colnames(peaks_by_samples) <- paste0(cbp_atac_info$genotype, "_", cbp_atac_info$cell_type)
  }
  peaks_by_samples
}


makeReadQCDataFrame <- function(file_path, peaks_by_samples_matrix){
  load(file = file_path)
  peaks_by_samples <- peaks_by_samples_matrix
  
  read_stats <- as.data.frame(count_reads_in_peaks$stat)
  rownames(read_stats) <- read_stats$Status
  read_stats <- read_stats[, -1]
  colnames(read_stats) <- colnames(peaks_by_samples)
  read_stats <- t(read_stats)
  
  
  sample_QC_df <- data.frame(sample_name = colnames(peaks_by_samples), 
                             reads_in_promoters = read_stats[, 1], 
                             mapped_reads = sapply(1:(dim(read_stats)[1]), function(xx) 
                               sum(read_stats[xx, ])) - read_stats[, 2],
                             total_reads = sapply(1:(dim(read_stats)[1]), function(xx) 
                               sum(read_stats[xx, ])))
  
  sample_QC_df$proportion_of_reads_in_promoters <- 100*round(sample_QC_df$reads_in_promoters / 
                                                               sample_QC_df$total_reads, 2)
  
  sample_QC_df
}

getPeakBreakdownInGenomicRegions <- function(peaks_granges, plot = c(TRUE, FALSE), main_lab){
  blacklisted_region_overlap <- peaks_granges[
    unique(queryHits(findOverlaps(peaks_granges, mm10_regions_to_exclude_granges)))]
  blacklisted_region_percentage <- 100*round(length(blacklisted_region_overlap)/length(peaks_granges), 2)
  
  proximal_promoter_overlap <- peaks_granges[
    unique(queryHits(findOverlaps(peaks_granges, proms_mouse)))]
  extended_promoter_overlap <- peaks_granges[
    unique(queryHits(findOverlaps(peaks_granges, proms_mouse_2)))]
  exon_overlap <- peaks_granges[
    unique(queryHits(findOverlaps(peaks_granges, all_exons_mouse_minus_promoters)))]
  other_overlap_length <- length(peaks_granges) - length(extended_promoter_overlap) - length(exon_overlap)
  
  
  proximal_promoter_overlap_percentage <- 100*round(length(proximal_promoter_overlap)/length(peaks_granges), 2)
  extended_promoter_overlap_percentage <- 100*round(length(extended_promoter_overlap)/length(peaks_granges), 2)
  exon_overlap_percentage <- 100*round(length(exon_overlap)/length(peaks_granges), 2)
  other_percentage <- 100 - extended_promoter_overlap_percentage - 
    exon_overlap_percentage - blacklisted_region_percentage
  
  c(proximal_promoter_overlap_percentage, extended_promoter_overlap_percentage, exon_overlap_percentage, 
    blacklisted_region_percentage)
  
  if (plot == TRUE){
    barplot(c(proximal_promoter_overlap_percentage, 
              extended_promoter_overlap_percentage, 
              exon_overlap_percentage, other_percentage, blacklisted_region_percentage), 
            ylab = "% of combined peaks", 
            names.arg = c(paste0("proximal\npromoter", "\n(n=", length(proximal_promoter_overlap), ")"),
                          paste0("extended\npromoter", "\n(n=", length(extended_promoter_overlap), ")"),
                          paste0("exon", "\n(n=", length(exon_overlap), ")"),
                          paste0("other\nregion", "\n(n=", other_overlap_length, ")"),
                          paste0("blacklisted\nregion", "\n(n=", length(blacklisted_region_overlap), ")")),  
            border = "white", las = 2, ylim = c(0, 100),
            col = rep(rgb(1,0,0,0.7), 5), 
            space = c(0.1, 0.1, 0.1, 0.1, 0.1),
            main = main_lab, yaxt = 'n', xlim = c(0, 5))
    axis(2, at = c(0, 20, 40, 60, 80, 100))
  }
}


getGrangesFromPeakAnnotation <- function(file_path){
  load(file = file_path)
  makeGRangesFromDataFrame(as.data.frame(count_reads_in_peaks$annotation[, 2:5]))
}

###
getDifferentialRegionGranges <- function(peaks_annotation_df, diff_results_tab){
  rownames(peaks_annotation_df) <- peaks_annotation_df$GeneID
  merged_df <- merge(diff_results_tab, peaks_annotation_df, by = 0, sort = FALSE, all.x = TRUE)
  merged_df$Row.names <- NULL
  
  diff_granges <- makeGRangesFromDataFrame(merged_df, keep.extra.columns = TRUE)
  genome(seqinfo(diff_granges)) <- "mm10"
  seqlevelsStyle(diff_granges) <- "ucsc"
  diff_granges
}

diff_proms_gene_ids <- unique(proms_mouse$gene_id[unique(queryHits(findOverlaps(
  proms_mouse, diff_prom_granges)))])


compareDEAnalyses <- function(de_tab_1, de_tab_2, pval_cutoff, what = c("p-value_density", "p_value_scatterplot", "log2FC"),
                                main_lab, x_lab, y_lab, legend_1, legend_2, x_lim, y_lim, bwidth){
  tab1 <- de_tab_1
  tab2 <- de_tab_2
  
  if (what == "p-value_density"){
    plot(density(tab1$P.Value[which(rownames(tab1) %in% rownames(tab2)[
      which(tab2$P.Value < pval_cutoff)])], from = 0, to = 1, bw = bwidth), col = alpha("brown", 0.75), lwd = 2.1, 
      xlab = x_lab, ylab = y_lab, bty = 'l', main = main_lab, font.main = 1, cex.lab = 1.2, xaxt = 'n', yaxt = 'n', ylim = y_lim)
    lines(density(tab1$P.Value[-which(rownames(tab1) %in% rownames(tab2)[
      which(tab2$P.Value < pval_cutoff)])], from = 0, to = 1, bw = bwidth), col = "cornflowerblue", lwd = 2.1)
    
    legend <- legend("topright", legend = c(legend_1, legend_2), 
                     bty = 'n', col = c(alpha("brown", 0.75), "cornflowerblue"), lty = "solid", lwd = 3, cex = 0.8)
    axis(1, at = c(0, 0.5, 1), cex.axis = 1)
    axis(2, at = y_lim, cex.axis = 1)
    
  } else if(what == "p_value_scatterplot"){
    tab1 <- tab1[which(rownames(tab1) %in% rownames(tab2)[
      which(tab2$P.Value < pval_cutoff)]), ]
    
    tab2 <- tab2[which(tab2$P.Value < pval_cutoff), ]
    no_overlap_indices <- which(!(rownames(tab2) %in% rownames(tab1)))
    if (length(no_overlap_indices) > 0){
      tab2 <- tab2[-no_overlap_indices, ] 
    }
    
    plot(tab2$P.Value[order(match(rownames(tab2), rownames(tab1)))], tab1$P.Value, 
         col = rgb(1,0,0,0.25), pch = 19, xlab = x_lab, ylab = y_lab, 
         main = main_lab, bty = 'l')
    
  } else if (what == "log2FC"){
    tab1 <- tab1[which(rownames(tab1) %in% rownames(tab2)[
      which(tab2$P.Value < pval_cutoff)]), ]
    
    tab2 <- tab2[which(tab2$P.Value < pval_cutoff), ]
    no_overlap_indices <- which(!(rownames(tab2) %in% rownames(tab1)))
    if (length(no_overlap_indices) > 0){
      tab2 <- tab2[-no_overlap_indices, ] 
    }
    
    plot(tab2$logFC[order(match(rownames(tab2), rownames(tab1)))], tab1$logFC, 
         col = rgb(1,0,0,0.25), pch = 19, xlab = x_lab, ylab = y_lab, 
         main = main_lab, bty = 'l', xlim = x_lim, ylim = y_lim, font.main = 1, xaxt = 'n', yaxt = 'n', cex.lab = 1.2, cex = 0.8)
    axis(1, at = c(x_lim[1], 0, x_lim[2]), cex.axis = 1)
    axis(2, at = c(y_lim[1], 0, y_lim[2]), cex.axis = 1)
    
    abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
    abline(h = 0, lty = "longdash", col = rgb(0,0,0,0.7))
  }
}

#the following two functions are useful for the blood samples
getGenotypeVector <- function(sample_names_vector){
  genotype <- rep(NA, length(sample_names_vector))
  genotype[grep("WT", sample_names_vector)] <- "WT"
  genotype[grep("KDM6A", sample_names_vector)] <- "KS2"
  genotype[grep("KMT2D", sample_names_vector)] <- "KS1"
  genotype[grep("CBP", sample_names_vector)] <- "RT1"
  genotype
}

getCellTypeVector <- function(sample_names_vector){
  cell_type <- rep(NA, length(sample_names_vector))
  cell_type[grep("_B", sample_names_vector)] <- "B"
  cell_type[grep("_T", sample_names_vector)] <- "T"
  cell_type
}

###
getEmpiricalNullPValues <- function(all_results_df){
  res <- all_results_df
  colnames(res)[c(3,4,5)] <- c("stat", "pvalue", "padj")
  res <- res[ !is.na(res$padj), ]
  res <- res[ !is.na(res$pvalue), ]
  res <- res[, -which(names(res) == "padj")]
  fdr.res <- fdrtool(res$stat, statistic= "normal", plot = FALSE)
  res[,"padj"]  <- p.adjust(fdr.res$pval, method = "BH")
  res$pvalue <- fdr.res$pval
  colnames(res)[c(4, 6)] <- c("P.Value", "adj.P.Val")
  res
}

getEmpiricalNullPValues2 <- function(all_results_df){
  res <- all_results_df
  colnames(res)[c(3,4,5)] <- c("stat", "pvalue", "padj")
  res <- res[ !is.na(res$padj), ]
  res <- res[ !is.na(res$pvalue), ]
  res <- res[, -which(names(res) == "padj")]
  fdr.res <- fdrtool(res$stat, statistic= "normal", plot = FALSE)
  res[,"padj"]  <- p.adjust(fdr.res$pval, method = "BH")
  res$pvalue <- fdr.res$pval
  colnames(res)[c(4, 6)] <- c("P.Value", "adj.P.Val")
  res
}


getDifferentialAnalysisResults <- function(features_by_samples_mat, genotype, 
      celltype = c("B", "T", "neurons"), data_type = c("RNA", "ATAC"), what = c("voom_matrix", "results_tab")){
  if (celltype == "B" | celltype == "T"){
    sample_info <- data.frame(genotype = getGenotypeVector(colnames(features_by_samples_mat)), 
                              cell_type = getCellTypeVector(colnames(features_by_samples_mat)))
    
    features_by_samples_mat <- features_by_samples_mat[, which(sample_info$genotype %in% c("WT", genotype) & 
                                                                 sample_info$cell_type == celltype)]
    
    sample_info2 <- sample_info[which(sample_info$genotype %in% c("WT", genotype) & 
                                        sample_info$cell_type == celltype), , drop = FALSE]
    sample_info2$genotype <- droplevels(sample_info2$genotype)
  } else if (celltype == "neurons"){
    sample_info <- data.frame(genotype = gsub(".*[_]", "", colnames(features_by_samples_mat)))
    
    features_by_samples_mat <- features_by_samples_mat[, which(sample_info$genotype %in% c("wt", genotype))]
    
    sample_info2 <- sample_info[which(sample_info$genotype %in% c("wt", genotype)), , drop = FALSE]
    sample_info2$genotype <- droplevels(sample_info2$genotype)
  }
  
  if (data_type == "RNA"){
    dgelist <- DGEList(features_by_samples_mat)
  } else if (data_type == "ATAC"){
    dgelist <- DGEList(counts = features_by_samples_mat, 
            genes = count_reads_in_peaks$annotation[, c("GeneID","Length")]) ###I need to have loaded the count_reads_in_peaks object before
  }

  mod <- model.matrix(~as.factor(genotype), data = sample_info2)
  mod0 <- model.matrix(~1, data = sample_info2)
  
  keep <- filterByExpr(dgelist, mod)
  dgelist <- dgelist[keep, , keep.lib.size=FALSE]
  
  nf <- calcNormFactors(dgelist$counts)
  transformed_counts <- dgelist$counts
  norm_factors <- nf
  
  dgelist <- voom(transformed_counts, mod, lib.size = colSums(transformed_counts)*norm_factors)
  
  sva <- sva(dgelist$E, mod, mod0)
  if(sva$n.sv > 0){
    modsv <- cbind(mod, sva$sv)
    
    dgelist2 <- voom(transformed_counts, modsv, lib.size = colSums(transformed_counts)*norm_factors)
    
    fit <- lmFit(dgelist2, modsv)
    contr.matrix = matrix(c(0, -1, rep(0,sva$n.sv)), ncol = 1)
  } else {
    dgelist2 <- dgelist
    fit <- lmFit(dgelist2, mod)
    contr.matrix = matrix(c(0, -1), ncol = 1)
  }
  
  if (what == "voom_matrix"){
    dgelist2$E
  } else if (what == "results_tab"){
    
    fit2 <- contrasts.fit(fit, contr.matrix)                     
    fit2 <- eBayes(fit2)
    tab <- topTable(fit2, coef = 1, n = Inf, adjust="fdr")
    
    tab
  }
}


getIHWPalues <- function(de_tab_1, de_tab_2, stratification_method = c("all", "threshold"), pval_threshold, 
                         alpha_level){ 
#this is a function that applies independent hypothesis weighting 
#using the p-values of one experiment (e.g. in B cells or in KS1) 
#as a covariate to increase power for the other experiment (e.g. in T cells or in KS2) 
  
  tab1 <- de_tab_1
  tab2 <- de_tab_2
  if (length(which(is.na(tab2$P.Value))) > 0){
    tab2 <- tab2[-which(is.na(tab2$P.Value)), ]
  }
  if (length(which(is.na(tab1$P.Value))) > 0){
    tab1 <- tab1[-which(is.na(tab1$P.Value)), ]
  }
  tab1 <- tab1[which(rownames(tab1) %in% rownames(tab2)), ]
                     
  no_overlap_indices <- which(!(rownames(tab2) %in% rownames(tab1)))
  if (length(no_overlap_indices) > 0){
    tab2 <- tab2[-no_overlap_indices, ] 
  }
  
  tab2 <- tab2[order(match(rownames(tab2), rownames(tab1))), ]
  
  if (stratification_method == "threshold"){
    ihw_res <- ihw(pvalues = tab1$P.Value, covariates = as.factor(tab2$P.Value < pval_threshold), 
                   alpha = alpha_level)
  } else if (stratification_method == "all"){
    ihw_res <- ihw(pvalues = tab1$P.Value, covariates = tab2$P.Value, 
                   alpha = alpha_level)
  }
  list(tab1, tab2, ihw_res)
}

plotDEStatisticsComparison <- function(de_object1, de_object2, what = c("variance", "logFC", "p-values"), 
                                       x_lab, y_lab, legend_1, legend_2, add = c(TRUE, FALSE), color){
  if (what == "variance"){
    mat1 <- de_object1
    mat2 <- de_object2
    
    mat1 <- mat1[
      which(rownames(mat1) %in% rownames(mat2)), ]
    
    
    mat1 <- mat1[order(match(rownames(mat1), rownames(mat2))), ]
    if (add == FALSE){
      plot(rowVars(mat2), rowVars(mat1), 
           xlim = c(0, max(max(rowVars(mat1)), max(rowVars(mat2)))), ylim = c(0, max(max(rowVars(mat1)), max(rowVars(mat2)))), 
           col = alpha("cornflowerblue", 0.1), pch = 19, 
           xlab = x_lab, ylab = y_lab, 
           main = "", bty = 'l')
      abline(0, 1)
    } else if (add == TRUE){
      points(rowVars(mat2), rowVars(mat1), pch = 19, 
           col = color)
    }

  } else if (what == "logFC"){
    tab1 <- de_object1
    tab2 <- de_object2
    
    tab1 <- tab1[
      which(rownames(tab1) %in% rownames(tab2)), ]
    
    plot(tab1$logFC, 
         tab2$logFC[order(match(rownames(tab2), rownames(tab1)))], 
         xlab = x_lab, ylab = y_lab, main = paste0("cor = ", cor(tab1$logFC, tab2$logFC[order(match(rownames(tab2), rownames(tab1)))])), 
         pch = 19, col = rgb(1,0,0,0.1))
    abline(0, 1)
  } else if (what == "p-values"){
    tab1 <- de_object1
    tab2 <- de_object2
    
    tab1 <- tab1[
      which(rownames(tab1) %in% rownames(tab2)), ]
    
    
    hist(tab1$P.Value, xlab = "p-values",
         col = 2, border = FALSE, main = "neurons (promoters only)")
    hist(tab2$P.Value, col = "cornflowerblue", add = TRUE, border = FALSE)
    
    legend <- legend("topright", legend = c(legend_1, legend_2), 
                     bty = 'n', fill = c(2, "cornflowerblue"), border = "white")
  }
}

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



diff_granges_B_KS1_atac <- getDifferentialRegionGranges(rownames(tab_B_KS1_atac))
diff_granges_B_KS1_atac$pval <- tab_B_KS1_atac$P.Value

overlaps <- findOverlaps(diff_granges_B_KS1_atac, proms_mouse)
diff_granges_B_KS1_atac <- diff_granges_B_KS1_atac[-queryHits(overlaps)[which(duplicated(queryHits(overlaps)))]]
diff_granges_B_KS1_atac$gene_id <- sapply(1:length(diff_granges_B_KS1_atac), function(xx) 
  proms_mouse$gene_id[unique(queryHits(findOverlaps(proms_mouse, diff_granges_B_KS1_atac[xx])))])
diff_granges_B_KS1_atac$gene_id <- as.character(diff_granges_B_KS1_atac$gene_id)
diff_granges_B_KS1_atac <- diff_granges_B_KS1_atac[-which(diff_granges_B_KS1_atac$gene_id == "character(0)")]



####
getIntersectionCSV <- function(tab1, tab2, output_file_name){
  new_pvals_tab1 <- getIHWPalues(tab1, tab2, "all", alpha_level = 0.1)
  new_tab1 <- new_pvals_tab1[[1]]
  new_adj_pvals_tab1 <- adj_pvalues(new_pvals_tab1[[3]])
  region_intersection <- intersect(rownames(tab2)[which(tab2$adj.P.Val < 0.1)], 
                                   rownames(new_tab1)[which(new_adj_pvals_tab1 < 0.1)])
  intersection_df <- t(as.data.frame(strsplit(region_intersection, "_")))
  rownames(intersection_df) <- paste0("common_differential_region_", 1:length(region_intersection))
  colnames(intersection_df) <- c("chr", "start", "end")
  write_csv(as.data.frame(intersection_df), output_file_name)
}

getCSVs <- function(tab1, tab2, output_file_names){
  new_pvals_tab1 <- getIHWPalues(tab1, tab2, "all", alpha_level = 0.1)
  new_tab1 <- new_pvals_tab1[[1]]
  new_adj_pvals_tab1 <- adj_pvalues(new_pvals_tab1[[3]])
  
  rownames2 <- rownames(tab2)[which(tab2$adj.P.Val < 0.1)]
  rownames2_df <- t(as.data.frame(strsplit(rownames2, "_")))
  rownames(rownames2_df) <- paste0("differential_region_", 1:length(rownames2))
  colnames(rownames2_df) <- c("chr", "start", "end")
  
  rownames1 <- rownames(new_tab1)[which(new_adj_pvals_tab1 < 0.1)]
  rownames1_df <- t(as.data.frame(strsplit(rownames1, "_")))
  rownames(rownames1_df) <- paste0("differential_region_", 1:length(rownames1))
  colnames(rownames1_df) <- c("chr", "start", "end")
  
  write_csv(as.data.frame(rownames1_df), output_file_names[1])
  write_csv(as.data.frame(rownames2_df), output_file_names[2])
}


getTabForIPA <- function(de_results_tab, output_name){
  df <- data.frame(gene_id = rownames(de_results_tab), logFC = de_results_tab$logFC, pval = de_results_tab$P.Value)
  colnames(df) <- c("Gene ID", "Log2Ratio", "p-value")
  write_tsv(df, output_name)
}


###########
getFullResultsDF <- function(differential_results_tab){
  results_granges <- getDifferentialRegionGranges(count_reads_in_peaks$annotation, as.data.frame(differential_results_tab))
  
  overlaps <- findOverlaps(results_granges, proms_mouse)
  overlaps_splitted <- split(overlaps, queryHits(overlaps))
  
  results_granges$prom_overlapping_id <- NA
  results_granges$prom_overlapping_id[unique(queryHits(overlaps))] <- lapply(overlaps_splitted, function(xx) {
    df <- as.data.frame(xx)
    unique(proms_mouse$gene_id[df$subjectHits])
  })
  
  results_granges$prom_overlapping_name <- NA
  results_granges$prom_overlapping_name[unique(queryHits(overlaps))] <- lapply(overlaps_splitted, function(xx) {
    df <- as.data.frame(xx)
    unique(proms_mouse$gene_name[df$subjectHits])
  })
  
  results_df <- values(results_granges)
  results_df$chr <- seqnames(results_granges)
  results_df$start <- start(results_granges)
  results_df$end <- end(results_granges)
  
  results_df <- results_df[, c(11:15, 1:2, 5:8)]
  results_df
}



getFullResultsDF2 <- function(differential_results_tab){
  results_granges <- getDifferentialRegionGranges(count_reads_in_peaks$annotation, as.data.frame(differential_results_tab))
  
  overlaps <- findOverlaps(results_granges, proms_mouse)
  overlaps_splitted <- split(overlaps, queryHits(overlaps))
  
  results_granges$prom_overlapping_id <- NA
  results_granges$prom_overlapping_id[unique(queryHits(overlaps))] <- lapply(overlaps_splitted, function(xx) {
    df <- as.data.frame(xx)
    unique(proms_mouse$gene_id[df$subjectHits])
  })
  
  results_granges$prom_overlapping_name <- NA
  results_granges$prom_overlapping_name[unique(queryHits(overlaps))] <- lapply(overlaps_splitted, function(xx) {
    df <- as.data.frame(xx)
    unique(proms_mouse$gene_name[df$subjectHits])
  })
  
  results_granges
}

getGOEnrichedTerms <- function(de_genes, assayed_genes){
  assayed.genes <- assayed_genes
  de.genes <- de_genes
  
  gene.vector <- as.integer(assayed.genes %in% de.genes)
  names(gene.vector) <- assayed.genes
  
  pwf <- nullp(gene.vector, 'mm10', 'ensGene')
  GO.wall <- goseq(pwf, "mm10", "ensGene")
  
  enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH") < 0.1]
  GO_enriched_names <- sapply(enriched.GO, function(xx) GOTERM[[xx]])
  Terms <- sapply(1:length(GO_enriched_names), function(xx) Term(GO_enriched_names[xx][[1]]))
  Terms
}

getReactomeEnrichedPathways <- function(de_genes, assayed_genes){
  assayed.genes <- assayed_genes
  de.genes <- de_genes
  
  ids <- bitr(assayed_genes, fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = "org.Mm.eg.db")
  ids <- ids[-which(ids$ENTREZID %in% ids$ENTREZID[which(duplicated(ids$ENTREZID))]), ]
  
  reactome_list <- as.list(reactomeEXTID2PATHID)
  ids <- ids[which(ids$ENTREZID %in% names(reactome_list)), ]
  
  reactome_pathways <- lapply(ids$ENTREZID, function(xx) reactome_list[[which(names(reactome_list) == xx)]])
  names(reactome_pathways) <- sapply(ids$ENTREZID, function(xx) ids$ENSEMBL[which(ids$ENTREZID == xx)])
  
  assayed.genes <- assayed.genes[which(assayed.genes %in% ids$ENSEMBL)]
  de.genes <- de_genes[which(de_genes %in% assayed.genes)]

  gene.vector <- as.integer(assayed.genes %in% de.genes)
  names(gene.vector) <- assayed.genes
  
  pwf <- nullp(gene.vector, 'mm10', 'ensGene')
  reactome <- goseq(pwf, gene2cat = reactome_pathways)
  
  all_pathway_names <- as.list(reactomePATHID2NAME)
  top_pathway_df <- data.frame(pathway = sapply(reactome$category, 
                                                function(xx) all_pathway_names[[which(names(all_pathway_names) == xx)]]), 
                               overrepresentation_pval = reactome$over_represented_pvalue)
  top_pathway_df
}

getPathwayGenes <- function(assayed_genes, pathway_reactome_id){
  ids <- bitr(assayed_genes, fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = "org.Mm.eg.db")
  ids <- ids[-which(ids$ENTREZID %in% ids$ENTREZID[which(duplicated(ids$ENTREZID))]), ]
  reactome_list <- as.list(reactomePATHID2EXTID)
  
  entrez_ids <- reactome_list[which(names(reactome_list) == pathway_reactome_id)][[1]]
  ensembl_to_entrez_df <- ids[which(ids$ENTREZID %in% entrez_ids), ]
  unique(ensembl_to_entrez_df$ENSEMBL)
}


