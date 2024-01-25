getOverlapPi0 <- function(ranges1, ranges2){ #ranges 1 must have a column named "pvalue"
  1 - pi0est(p = ranges1$pvalue[unique(queryHits(findOverlaps(ranges1, ranges2)))], 
             pi0.method = "bootstrap")$pi0
}

getOverlapPi02 <- function(ranges1, ranges2, lambda_value){ #ranges 1 must have a column named "pvalue"
  1 - pi0est(p = ranges1$pvalue[unique(queryHits(findOverlaps(ranges1, ranges2)))], 
             lambda = lambda_value)$pi0
}

getOverlapNullDistribution1 <- function(ranges1, ranges2){ #ranges 1 must have a column named "pvalue"
  replicate(10000, {
    pi0_random <- pi0est(p = sample(ranges1$pvalue, 
                                    length(unique(queryHits(findOverlaps(ranges1, ranges2))))), 
                         pi0.method = "bootstrap")$pi0
    1- pi0_random
    #length(which(qobj_random$significant == TRUE))
  })
}

getOverlapNullDistribution2 <- function(ranges1, ranges2, pval_vector, threshold){ #ranges 1 must have a column named "pvalue"
                                                           #pval_vector can be either adjusted or unadjusted pvals
  replicate(10000, {
    pi0_random <- pi0est(p = sample(ranges1$pvalue[unique(queryHits(findOverlaps(ranges1, ranges2)))], 
                                    length(unique(queryHits(findOverlaps(ranges1, 
                                                                         ranges2[which(pval_vector < threshold)]))))), 
                         pi0.method = "bootstrap")$pi0
    1- pi0_random
    #length(which(qobj_random$significant == TRUE))
  })
}

getOverlapNullDistribution3 <- function(ranges1, ranges2, pval_vector, threshold, lambda_value){ #ranges 1 must have a column named "pvalue"
  #pval_vector can be either adjusted or unadjusted pvals
  replicate(10000, {
    pi0_random <- pi0est(p = sample(ranges1$pvalue[unique(queryHits(findOverlaps(ranges1, ranges2)))], 
                                    length(unique(queryHits(findOverlaps(ranges1, 
                                                                         ranges2[which(pval_vector < threshold)]))))), 
                         lambda = lambda_value)$pi0
    1- pi0_random
    #length(which(qobj_random$significant == TRUE))
  })
}


getWilcoxStat <- function(ranges1, ranges2){
  indices <- unique(queryHits(findOverlaps(ranges1, ranges2)))
  wilcox.test(ranges1$pvalue[indices], ranges1$pvalue[-indices])$statistic
}

getWilcoxStat2 <- function(full_df, gene_ids){
  indices <- which(full_df$gene_ids %in% gene_ids)
  wilcox.test(full_df$pval[indices], full_df$pval[-indices])$statistic
}

getWilcoxNullDistribution <- function(ranges1, ranges2){
  length_overlap <- length(unique(queryHits(findOverlaps(ranges1, ranges2))))
  permutation_rank <- replicate(1000, {
    indices <- sample(1:length(ranges1), length_overlap)
    wilcox.test(ranges1$pvalue[indices], ranges1$pvalue[-indices])$statistic
  })
}

getWilcoxNullDistribution2 <- function(full_df, gene_ids){
  length_overlap <- length(which(full_df$gene_ids %in% gene_ids))
  permutation_rank <- replicate(1000, {
    full_length <- dim(full_df)[1]
    indices <- sample(1:full_length, length_overlap)
    wilcox.test(full_df$pval[indices], full_df$pval[-indices])$statistic
  })
}

getPercentInTopPval <- function(full_df, gene_ids, threshold){
  top_threshold <- quantile(full_df$pval, threshold)
  indices <- which(full_df$gene_ids %in% gene_ids)
  length(which(full_df$pval[indices] < top_threshold))/length(indices)
}

getRandomPercentInTopPval <- function(full_df, gene_ids, threshold){
  indices <- which(full_df$gene_ids %in% gene_ids)
  random_gene_ids <- sample(full_df$gene_ids[-indices], length(indices))
  getPercentInTopPval(full_df, random_gene_ids, threshold)
}

compareLogFC <- function(de_results1, de_results2, x_lab, y_lab, point_color){ #results must be in granges form
  overlaps <- findOverlaps(de_results1, de_results2)
  logfc_products <- de_results1$log2FoldChange[queryHits(overlaps)]*de_results2$log2FoldChange[subjectHits(overlaps)]
  sign_corcondance <- prop.table(table(as.factor(logfc_products > 0)))[2]
  #plot(de_results1$log2FoldChange[queryHits(overlaps)], de_results2$log2FoldChange[subjectHits(overlaps)], 
  #     xlab = x_lab, ylab = y_lab, main = paste0("cor = ", round(cor(de_results1$log2FoldChange[queryHits(overlaps)], 
  #                                                             de_results2$log2FoldChange[subjectHits(overlaps)]), 2)), 
  #     pch = 19, col = point_color)
  plot(de_results1$log2FoldChange[queryHits(overlaps)], de_results2$log2FoldChange[subjectHits(overlaps)], 
       xlab = x_lab, ylab = y_lab, main = paste0("sign concordance = ", round(sign_corcondance, 2)), 
       pch = 19, col = point_color)
  abline(h = 0)
  abline(v = 0)
}

getRandomEffectSizeConcordance <- function(de_results1, de_results2, pval_vector, threshold){
  de_results2_below_cutoff <- de_results2[which(pval_vector < threshold)]
  
  n_pos <- length(which(de_results2_below_cutoff$log2FoldChange > 0))
  n_neg <- length(which(de_results2_below_cutoff$log2FoldChange < 0))
  
  ind_pos_all <- which(de_results2$log2FoldChange > 0)
  ind_pos_sampled <- sample(ind_pos_all, n_pos)
  
  ind_neg_all <- which(de_results2$log2FoldChange < 0)
  ind_neg_sampled <- sample(ind_neg_all, n_neg)
  
  overlaps <- findOverlaps(de_results1, de_results2[c(ind_pos_sampled, ind_neg_sampled)])
  logfc_products <- de_results1$log2FoldChange[queryHits(overlaps)]*de_results2$log2FoldChange[subjectHits(overlaps)]
  sign_concordance <- prop.table(table(as.factor(logfc_products > 0)))[2]
  sign_concordance
}

getEffectSizeConcordance <- function(de_results1, de_results2){
  overlaps <- findOverlaps(de_results1, de_results2)
  logfc_products <- de_results1$log2FoldChange[queryHits(overlaps)]*de_results2$log2FoldChange[subjectHits(overlaps)]
  sign_concordance <- prop.table(table(as.factor(logfc_products > 0)))[2]
  sign_concordance
}

getStatisticInExp2 <- function(de_results1, de_results2, what = c("mean", "variance", "logFC")){
  overlaps <- findOverlaps(de_results1, de_results2)
  if (what == "mean"){summary(de_results1$baseMean[queryHits(overlaps)])}
  else if (what == "variance"){summary(de_results1$lfcSE[queryHits(overlaps)])}
  else if (what == "logFC"){summary(de_results1$log2FoldChange[queryHits(overlaps)])}
}


