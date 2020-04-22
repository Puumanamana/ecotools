################################
########## NMDS utils ##########
################################

nmds <- function(dist_matrix, metadata=NULL, info.label='factor', info.value='', ndims=2, cores=50) {
  #' Wrapper around metaMDS to only compute scores for a subset of the distance matrix
  
  distances <- dist_matrix
  
  if(info.value != '') {
    distances <- dist_subset(dist_matrix, metadata[[info.label]]==info.value)
  }
  
  nmds.obj <- metaMDS(distances, k=ndims, parallel=cores, trymax=200, noshare=0.1)
  components <- as.data.frame(scores(nmds.obj))
  components[[info.label]] = info.value
  
  colnames(components) <- c(paste0('NMDS', 1:(ncol(components)-1)), info.label)
  
  return(list(obj=nmds.obj, data=components))
}

nmds_levels <- function(dist_matrix, ndims=2, metadata=NULL, factor=NULL) {
  #' - Get NMDS components for each level of [factor] 
  #' using the distance matrix in the outer scope
  #' - Returns matrix of both components for each level
  
  if (is.null(factor)) {
    return(nmds(dist_matrix, ndims=ndims)$data)
  }
  
  lvls <- levels(metadata[[factor]])
  
  nmds_results <- lapply(lvls, 
                         function(x) nmds(dist_matrix, metadata=metadata, ndims=ndims,
                                          info.label=factor, info.value=x)$data)
  nmds_results <- do.call(rbind, nmds_results)
  
  return(nmds_results)
}

################################
######### Adonis tests #########
################################

# make.pairwise.table <- function(phylo_data, strat_cond, factor) {
# 
#   phylo_data.subset <- subset_samples(phylo_data, strat_cond)
# 
#   model <- pairwise.adonis(t(phylo_data.subset@otu_table), 
#                            phylo_data.subset@sam_data[[factor]]
#   
#   levels_names <- unique(phylo_data.subset@sam_data[[factor]])
# 
#   pvals.df <- data.frame(matrix(ncol=length(levels_names), 
#                                 nrow=length(levels_names)))
#   colnames(pvals.df) <- levels_names
#   rownames(pvals.df) <- levels_names
  
#   stat.df <- data.frame(pvals.df)
#   r2.df <- data.frame(pvals.df)
  
#   for (col in colnames(pvals.df)) {
#     for (row in rownames(pvals.df)) {
#       if (is.na(pvals.df[row, col])) {
#       
#         factor.cond <- (
#           (phylo_data.subset@sam_data[[factor]] == row) |
#           (phylo_data.subset@sam_data[[factor]] == col)
#           )
#       
#         print(sprintf('%s vs %s: %s items', col, row, sum(factor.cond)))
#       
#         if (sum(factor.cond) == 0) {
#           next
#         }
#   
#         phylo_data.subset.lvl <- subset_samples(phylo_data.subset, factor.cond)
#         otu_table.subset.lvl <- as.data.frame(as.matrix(t(phylo_data.subset.lvl@otu_table)))
#         meta.subset.lvl <- as.data.frame(phylo_data.subset.lvl@sam_data)[[factor]]
#       
#         model <- adonis(otu_table.subset.lvl ~ meta.subset.lvl)
# 
#         result <- c(stat=model$aov.tab$F.Model[1], 
#                     R=model$aov.tab$R2[1], 
#                     p=model$aov.tab$`Pr(>F)`[1])
#       
#         ## res <- make.pairwise.test(phylo_data.subset.lvl)
#         pvals.df[row, col] = result['p']
#         pvals.df[col, row] = result['p']
#         r2.df[col, row] = result['R']
#         r2.df[row, col] = result['R']
#         stat.df[col, row] = result['stat']
#         stat.df[row, col] = result['stat']
#       }
#     }
#   }
  
#   params <- sprintf("adonis_%s-%s_by-%s", strat, strat_val, factor)
#   
#   write.table(pvals.df, file=paste0('results/pval_', params, '.tsv'), 
#               quote=FALSE, row.names=TRUE, col.names=NA, sep='\t')
#   write.table(r2.df, file=paste0('results/r2_', params, '.tsv'), 
#               quote=FALSE, row.names=TRUE, col.names=NA, sep='\t')
#   write.table(stat.df, file=paste0('results/statistic_', params, '.tsv'), 
#               quote=FALSE, row.names=TRUE, col.names=NA, sep='\t')  
# }
