library(tidyr)
library(data.table)
library(phyloseq)
library(vegan)
library(rhdf5)
library(ape)

################################
######### Data loaders #########
################################

fast_table_load <- function(filename, cores=1, drop=c(), row.names=1) {
    #' Wrapper around fread to significantly
    #' increase .csv reading speed
  
    data <- as.data.frame(fread(filename, nThread=cores, drop=drop, header=T, blank.lines.skip=T))
    
    if (row.names > 0) {
        rownames(data) <- data[[row.names]]
        data[, row.names] <- NULL
    }
    return(data)
}

load_dist <- function(name=NULL, data=NULL) {

    if (!file.exists(name)) {
        dist <- distance(data, "bray")
        saveRDS(dist, name)
    } else {
        dist <- readRDS(name)
    }
    return(dist)
}

################################
######### Main loader ##########
################################

assr <- function(x) {
    return(asin(sqrt(x/sum(x))))
}

sqrtsum <- function(x) {
    return(sqrt(x/sum(x)))
}

ratio <- function(x) {
    return(x/sum(x))
}

load_data_h5 <- function(h5_path, tree_path) {

    samples <- h5read(h5_path, 'samples')
    OTU.labels <- h5read(h5_path, 'OTUs')

    ## Phylogenetic tree
    tree <- read.tree(tree_path)
  
    ## abundance
    community_matrix <- h5read(h5_path, 'abundance')
    colnames(community_matrix) <- samples
    rownames(community_matrix) <- OTU.labels
  
    ## taxonomy
    taxonomy <- t(h5read(h5_path, 'taxa'))
    colnames(taxonomy) <- h5read(h5_path, 'ranks')
    rownames(taxonomy) <- OTU.labels
  
    ## metadata
    meta.quant <- as.data.frame(t(h5read(h5_path, 'metadata_quant')))
    colnames(meta.quant) <- h5read(h5_path, 'quant_vars')
    rownames(meta.quant) <- samples
    
    meta.qual <- as.data.frame(t(h5read(h5_path, 'metadata_qual')))
    colnames(meta.qual) <- h5read(h5_path, 'qual_vars')
    rownames(meta.qual) <- samples

    meta <- cbind(meta.quant, meta.qual)
    meta[meta == 'nan'] <- NA
    meta <- droplevels(meta)
  
    data <- phyloseq(
        otu_table(community_matrix, taxa_are_rows=T),
        sample_data(meta),
        tax_table(taxonomy),
        phy_tree(tree)
    )

    return(data)
}

filter_data <- function(ps, group='Genus', prevalenceThreshold=5) {
  ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
  
  ## Compute prevalence of each feature, store as data.frame
  info.df <- data.frame(
      abundance=taxa_sums(ps0),
      prevalence=rowSums(otu_table(ps0) > 0),
      taxonomy=tax_table(ps0)[, "Phylum"]
  )  
  ## Execute prevalence filter, using `prune_taxa()` function
  keepTaxa <- rownames(info.df)[(info.df$prevalence >= prevalenceThreshold)]
  ps2 <- prune_taxa(keepTaxa, ps0)
  
  ## How many genera would be present after filtering?
  length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
  
  if (!is.null(group)) {
      ps2 <- tax_glom(ps2, group, NArm=TRUE)  
  }

  return(ps2)
}

load_and_preproc() <- function() {
    data <- load_data_h5(h5_path, tree_path)
    data <- filter_data(data, group='Genus')
    data <- transform_sample_counts(raw, NORM_FN)

    return(data)
}
